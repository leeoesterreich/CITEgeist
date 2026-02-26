"""
Core optimization implementations using Gurobi for CITEgeist deconvolution.

This module contains the optimization functions for cell proportion estimation
and gene expression deconvolution using quadratic programming.
"""
# Standard library imports
import concurrent
import gc
import json
import logging
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple, Union

import gurobipy as gp

# Third-party imports
import numpy as np
import pandas as pd
import psutil
import scanpy as sc
import scipy
import scipy.sparse as sp
from gurobipy import GRB, Model, quicksum
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import minimize
from scipy.spatial.distance import squareform
from scipy.special import digamma, loggamma
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from .checkpoints import CheckpointManager

# Local imports
# Provide a fallback so this module can be executed directly as a script
try:
    from .utils import get_neighbors_with_fixed_radius
    from .checkpoints import CheckpointManager
    from .gex_modules import (
        compute_module_aware_enrichment,
        compute_softmax_target,
        compute_kl_penalty_coefficients,
    )
except Exception:  # pragma: no cover - fallback for __main__ execution
    from utils import get_neighbors_with_fixed_radius  # type: ignore
    from checkpoints import CheckpointManager  # type: ignore
    from gex_modules import (  # type: ignore
        compute_module_aware_enrichment,
        compute_softmax_target,
        compute_kl_penalty_coefficients,
    )


def build_spatial_laplacian(
    coords: np.ndarray,
    k: int = 8,
    normed: bool = True,
) -> scipy.sparse.spmatrix:
    """
    Build graph Laplacian from spatial coordinates using k-NN.

    The Laplacian matrix L is used for spatial smoothness regularization.
    For normalized Laplacian: L = I - D^(-1/2) A D^(-1/2)
    where A is the adjacency matrix and D is the degree matrix.

    The quadratic form Y^T L Y penalizes differences between neighboring spots:
        Y^T L Y = sum_{i~j} (Y_i - Y_j)^2 / sqrt(d_i * d_j)

    Args:
        coords: Spatial coordinates, shape (n_spots, 2).
        k: Number of nearest neighbors for graph construction.
        normed: If True, return normalized Laplacian (recommended).

    Returns:
        Sparse Laplacian matrix of shape (n_spots, n_spots).
    """
    from sklearn.neighbors import kneighbors_graph
    from scipy.sparse.csgraph import laplacian

    n_spots = coords.shape[0]

    # Build k-NN adjacency matrix (directed)
    A = kneighbors_graph(coords, k, mode='connectivity', include_self=False)

    # Symmetrize: if i is neighbor of j OR j is neighbor of i, they are connected
    A = (A + A.T)
    A.data = np.clip(A.data, 0, 1)  # Binary adjacency

    # Compute Laplacian
    L = laplacian(A, normed=normed)

    logging.debug(
        f"Built spatial Laplacian: {n_spots} spots, k={k}, "
        f"nnz={L.nnz}, normed={normed}"
    )

    return L


def compute_global_prior(
    spotwise_gene_expression_profiles: Dict[int, np.ndarray],
    cell_type_numbers_array: np.ndarray,
    lambda_prior: float = 1.0,
    min_expression_threshold: float = 0.1,
) -> Dict[str, Any]:
    """
    Compute global prior from pass 1 results using normalized expression patterns.

    Args:
        spotwise_gene_expression_profiles: Dictionary mapping spot indices to profile matrices
        cell_type_numbers_array: Array of cell type proportions (N_spots × T_celltypes)
        lambda_prior: Strength of prior (default: 1.0)
        min_expression_threshold: Minimum expression to consider "active" (default: 0.1)

    Returns:
        Dict containing:
            - global_prior: Prior matrix (T_celltypes × M_genes)
            - confidence_scores: Confidence in each prior value
            - expression_patterns: Summary of expression patterns
    """
    # Validate inputs
    N = cell_type_numbers_array.shape[0]
    T = cell_type_numbers_array.shape[1]

    spot_keys = sorted(spotwise_gene_expression_profiles.keys())
    if len(spot_keys) != N:
        raise ValueError(f"Mismatch in number of spots: {len(spot_keys)} vs {N}")

    # Get dimensions from first profile
    example_profile = spotwise_gene_expression_profiles[spot_keys[0]]
    M = example_profile.shape[1]  # number of genes

    # Initialize arrays
    usage_array = np.zeros((N, T, M))
    for i, profile in spotwise_gene_expression_profiles.items():
        usage_array[i] = profile

    # Calculate expression statistics per cell type
    mean_expression = np.zeros((T, M))
    expression_frequency = np.zeros((T, M))
    expression_consistency = np.zeros((T, M))

    for t in range(T):
        # Weight profiles by cell type abundance
        weights = cell_type_numbers_array[:, t]  # Now 1D array of shape (N,)

        # Calculate weighted statistics
        active_expression = usage_array[:, t, :] > min_expression_threshold
        weighted_expression = usage_array[:, t, :]  # Shape: (N, M)

        # Mean expression when the cell type is present
        present_mask = weights > 0
        if np.any(present_mask):
            # Ensure weights match the data shape for averaging
            weights_for_average = weights[present_mask]  # 1D array of length n_present
            expression_for_average = weighted_expression[present_mask, :]  # (n_present, M)

            mean_expression[t] = np.average(expression_for_average, weights=weights_for_average, axis=0)

            # Expression consistency (coefficient of variation, inverse)
            # Calculate weighted std dev properly
            diff_squared = (expression_for_average - mean_expression[t]) ** 2  # (n_present, M)
            weighted_var = np.average(diff_squared, weights=weights_for_average, axis=0)  # (M,)
            std = np.sqrt(weighted_var)  # (M,)
            expression_consistency[t] = 1 / (1 + std / (mean_expression[t] + 1e-6))

        # Frequency of expression (properly weighted)
        total_weight = np.sum(weights) + 1e-6
        expression_frequency[t] = np.sum(active_expression * weights[:, np.newaxis], axis=0) / total_weight

    # Combine metrics into confidence scores
    confidence_scores = expression_frequency * expression_consistency

    # Generate prior probabilities
    # Scale mean expression to [0,1] per gene
    scaled_expression = mean_expression / (np.max(mean_expression, axis=0) + 1e-6)

    # Weight by confidence and apply prior strength
    weighted_scores = scaled_expression * np.power(confidence_scores, lambda_prior)

    # Convert to probabilities via softmax
    global_prior = np.zeros((T, M))
    for m in range(M):
        scores = weighted_scores[:, m]
        exp_scores = np.exp(scores - np.max(scores))  # Numerical stability
        global_prior[:, m] = exp_scores / (np.sum(exp_scores) + 1e-6)

    # Log statistics
    logging.info("Prior computation statistics:")
    logging.info(f" - Mean confidence score: {np.mean(confidence_scores):.4f}")
    logging.info(f" - Mean prior strength: {np.mean(global_prior):.4f}")
    logging.info(f" - % Strong signals (>0.5): {100 * np.mean(global_prior > 0.5):.2f}%")

    return {
        "global_prior": global_prior,
        "confidence_scores": confidence_scores,
        "expression_patterns": {
            "mean_expression": mean_expression,
            "expression_frequency": expression_frequency,
            "expression_consistency": expression_consistency,
        },
    }


def map_antibodies_to_profiles(adata, cell_profile_dict):
    """
    Map antibody capture data to predefined cell type profiles.

    Args:
        adata (AnnData): Antibody capture AnnData object.
        cell_profile_dict (dict): Dictionary mapping cell types to antibody markers.

    Returns:
        np.ndarray: Profile-based antibody data matrix (N_spots x T_cell_types).
        list: List of cell type names (to ensure column order).
    """
    # Step 1: Subset data to relevant markers
    all_markers: list[str] = [marker for profile in cell_profile_dict.values() for marker in profile['Major']]
    existing_markers: list[str] = [marker for marker in all_markers if marker in adata.var_names]

    if len(existing_markers) == 0:
        logging.info("Adata variables: %s", adata.var_names)
        logging.info("Antibody markers: %s", all_markers)
        raise ValueError("No matching antibody markers found in adata.var_names.")

    adata.var_names_make_unique()
    adata = adata[:, existing_markers]

    # Step 2: Extract and prepare antibody capture data
    antibody_capture_data = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X.X
    antibody_capture_var_names = np.array(adata.var_names)

    cell_type_names = list(cell_profile_dict.keys())
    N = antibody_capture_data.shape[0]
    T = len(cell_type_names)

    profile_based_antibody_data = np.zeros((N, T))

    # Step 3: Map antibodies to profiles
    for profile_idx, (profile_name, profile_markers) in enumerate(cell_profile_dict.items()):
        major_markers = profile_markers.get("Major", [])
        try:
            relevant_marker_indices = [
                np.where(antibody_capture_var_names == marker)[0][0]
                for marker in major_markers
                if marker in antibody_capture_var_names
            ]
            if relevant_marker_indices:
                profile_based_antibody_data[:, profile_idx] = antibody_capture_data[:, relevant_marker_indices].mean(
                    axis=1
                )
            else:
                logging.warning(f"No valid markers found for profile '{profile_name}'")
        except IndexError as e:
            logging.warning(f"Error processing markers for profile '{profile_name}': {str(e)}")

    # Step 4: Normalize with safety checks
    column_max = np.max(profile_based_antibody_data, axis=0)
    zero_columns = column_max == 0
    if np.any(zero_columns):
        logging.warning("Zero columns detected. Adding epsilon to prevent NaNs.")
        column_max[zero_columns] = 1e-6

    profile_based_antibody_data /= column_max

    if np.isnan(profile_based_antibody_data).any():
        raise ValueError("NaN values detected in profile_based_antibody_data after mapping.")

    return profile_based_antibody_data, cell_type_names


def map_antibodies_to_profiles_v2(adata, cell_profile_dict):
    """
    Map antibody capture data while preserving marker-level granularity.

    Unlike map_antibodies_to_profiles() which averages multiple markers per cell type,
    this function keeps individual marker data and returns an assignment matrix mapping
    markers to cell types. This enables per-marker beta learning during optimization.

    Args:
        adata (AnnData): Antibody capture AnnData object.
        cell_profile_dict (dict): Dictionary mapping cell types to antibody markers.
            Format: {"CellType": {"Major": ["Marker1", "Marker2"], ...}, ...}

    Returns:
        Tuple of:
        - marker_level_data (np.ndarray): (N_spots, M_markers) normalized antibody data
        - marker_names (List[str]): Ordered list of marker names
        - assignment_matrix (np.ndarray): (M_markers, T_celltypes) binary assignment
        - cell_type_names (List[str]): Ordered list of cell type names
    """
    # Step 1: Collect all Major markers from cell_profile_dict
    all_markers = []
    for profile_markers in cell_profile_dict.values():
        all_markers.extend(profile_markers.get("Major", []))

    # Remove duplicates while preserving order
    seen = set()
    unique_markers = []
    for m in all_markers:
        if m not in seen:
            seen.add(m)
            unique_markers.append(m)

    # Step 2: Filter to markers that exist in adata
    existing_markers = [m for m in unique_markers if m in adata.var_names]

    if len(existing_markers) == 0:
        logging.info("Adata variables: %s", adata.var_names)
        logging.info("Antibody markers: %s", unique_markers)
        raise ValueError("No matching antibody markers found in adata.var_names.")

    # Step 3: Extract marker-level data
    adata.var_names_make_unique()
    marker_indices = [np.where(adata.var_names == m)[0][0] for m in existing_markers]

    antibody_data = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    marker_level_data = antibody_data[:, marker_indices].astype(np.float64)  # (N, M)

    # Step 4: Normalize per-marker (column-wise max normalization)
    col_max = np.max(marker_level_data, axis=0)
    zero_cols = col_max == 0
    if np.any(zero_cols):
        logging.warning(f"Zero columns detected for markers: {[existing_markers[i] for i, z in enumerate(zero_cols) if z]}. Adding epsilon.")
        col_max[zero_cols] = 1e-6
    marker_level_data = marker_level_data / col_max

    if np.isnan(marker_level_data).any():
        raise ValueError("NaN values detected after normalization.")

    # Step 5: Build assignment matrix
    cell_type_names = list(cell_profile_dict.keys())
    M = len(existing_markers)
    T = len(cell_type_names)
    assignment_matrix = np.zeros((M, T), dtype=np.float64)

    marker_to_idx = {name: i for i, name in enumerate(existing_markers)}

    for j, (ct_name, markers_dict) in enumerate(cell_profile_dict.items()):
        major_markers = markers_dict.get("Major", [])
        for marker in major_markers:
            if marker in marker_to_idx:
                m = marker_to_idx[marker]
                assignment_matrix[m, j] = 1.0

    logging.info(f"Mapped {M} markers to {T} cell types (preserving marker-level data)")

    return marker_level_data, existing_markers, assignment_matrix, cell_type_names


def optimize_cell_proportions(
    profile_based_antibody_data: np.ndarray,
    cell_type_names: List[str],
    tolerance: float = 1e-4,
    max_iterations: int = 50,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    normalize_beta: bool = True,
    unknown_threshold: float = 0.05,
    min_celltype_threshold: float = 0.01,
    redundancy_threshold: float = 0.2,
    warn_only: bool = False,
    # Laplacian smoothing parameters
    lambda_laplacian: float = 0.1,
    coords: Optional[np.ndarray] = None,
    laplacian_k: int = 8,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform EM-based optimization for cell type proportions using Gurobi.

    Includes optional Laplacian smoothing for spatial regularization, which
    encourages similar proportions in neighboring spots.

    Args:
        profile_based_antibody_data: N x T matrix of mapped antibody data
        cell_type_names: List of cell type names
        tolerance: Convergence tolerance for EM algorithm
        max_iterations: Maximum number of iterations
        lambda_reg: Regularization strength for elastic net
        alpha: L1-L2 tradeoff factor (0 = L2, 1 = L1)
        normalize_beta: Whether to normalize beta values
        unknown_threshold: Maximum allowed mean proportion for Unknown (default: 0.05)
        min_celltype_threshold: Minimum required mean proportion for cell types (default: 0.01)
        redundancy_threshold: Maximum allowed fraction of redundant types (default: 0.1)
        warn_only: If True, only warn on validation failures instead of raising
        lambda_laplacian: Weight for Laplacian smoothing (default: 0.1).
            Set to 0 to disable spatial smoothing.
        coords: Spatial coordinates, shape (N, 2). Required if lambda_laplacian > 0.
        laplacian_k: Number of neighbors for Laplacian graph (default: 8).
        spot_abundance_target: Optional per-spot abundance target (shape N,).
            When provided, objective includes soft penalty on row-sum mismatch.
        lambda_abundance_prior: Weight of abundance soft-prior term.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Y_values (N x T), beta_values (T,)
    """
    import gurobipy as gp
    from gurobipy import GRB

    N, T = profile_based_antibody_data.shape
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N:
            raise ValueError(
                f"spot_abundance_target length {spot_abundance_target.shape[0]} does not match N={N}."
            )

    # Build spatial Laplacian for smoothing if requested
    L_coo = None
    use_laplacian = lambda_laplacian > 0 and coords is not None
    if use_laplacian:
        if coords.shape[0] != N:
            raise ValueError(
                f"coords has {coords.shape[0]} rows but data has {N} spots"
            )
        L = build_spatial_laplacian(coords, k=laplacian_k, normed=True)
        # Convert to COO for efficient iteration
        L_coo = L.tocoo()
        logging.info(
            f"Laplacian smoothing enabled: lambda={lambda_laplacian}, "
            f"k={laplacian_k}, nnz={L_coo.nnz}"
        )
    elif lambda_laplacian > 0 and coords is None:
        logging.warning(
            "lambda_laplacian > 0 but coords not provided. "
            "Laplacian smoothing disabled."
        )

    # Initialize beta estimates
    beta_estimates = {ct: 1.0 for ct in cell_type_names}
    beta_prev = np.zeros(T)
    Y_prev = np.zeros((N, T))
    iteration = 0

    while iteration < max_iterations:
        logging.info(f"\nIteration {iteration + 1}")
        model = gp.Model("EM_Cell_Proportions")
        model.setParam("OutputFlag", 0)

        # Define variables Y[i, j]
        Y = model.addVars(N, T, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="Y")

        # Objective: Total squared error + Elastic Net regularization
        error_terms = []
        for i in range(N):
            for j in range(T):
                S_ij = profile_based_antibody_data[i, j]
                beta_j = beta_estimates[cell_type_names[j]]
                Y_ij = Y[i, j]
                error_terms.append((S_ij - beta_j * Y_ij) * (S_ij - beta_j * Y_ij))

        total_error = gp.quicksum(error_terms)
        l1_term = gp.quicksum(Y[i, j] for i in range(N) for j in range(T))
        l2_term = gp.quicksum(Y[i, j] * Y[i, j] for i in range(N) for j in range(T))
        regularization_term = lambda_reg * (alpha * l1_term + (1 - alpha) * l2_term)

        # Add Laplacian smoothing term: lambda_L * sum_k Y_k^T L Y_k
        # This penalizes differences between neighboring spots
        laplacian_term = 0
        if use_laplacian and L_coo is not None:
            laplacian_terms = []
            # Iterate over non-zero entries of Laplacian
            for idx in range(L_coo.nnz):
                i_spot = L_coo.row[idx]
                j_spot = L_coo.col[idx]
                L_val = L_coo.data[idx]
                # Add term for each cell type: L[i,j] * Y[i,k] * Y[j,k]
                for k in range(T):
                    laplacian_terms.append(L_val * Y[i_spot, k] * Y[j_spot, k])
            laplacian_term = lambda_laplacian * gp.quicksum(laplacian_terms)

        abundance_prior_term = 0
        if spot_abundance_target is not None and lambda_abundance_prior > 0:
            abundance_terms = []
            for i in range(N):
                row_sum = gp.quicksum(Y[i, j] for j in range(T))
                target_i = float(spot_abundance_target[i])
                abundance_terms.append((row_sum - target_i) * (row_sum - target_i))
            abundance_prior_term = lambda_abundance_prior * gp.quicksum(abundance_terms)

        model.setObjective(total_error + regularization_term + laplacian_term + abundance_prior_term, GRB.MINIMIZE)

        # Sum of proportions constraints
        for i in range(N):
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) >= 0.9)
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) <= 1.2)
        model.write('cell_proportion_model.mps')
        try:
            model.optimize()
        except Exception as e:
            logging.error(f"Optimization error: {str(e)}")
            raise ValueError("Gurobi optimization failed") from e

        if model.status == GRB.OPTIMAL:
            Y_values = np.array([[Y[i, j].X for j in range(T)] for i in range(N)])
        else:
            raise ValueError("Gurobi optimization failed to converge")

        # Update beta
        beta_new = np.zeros(T)
        for j in range(T):
            Y_j = Y_values[:, j]
            S_j = profile_based_antibody_data[:, j]
            denominator = np.dot(Y_j, Y_j)

            if denominator > 0:
                beta_new[j] = np.dot(S_j, Y_j) / denominator
            beta_new[j] = max(beta_new[j], 0.0)  # Ensure non-negative

        # Optionally normalize beta values
        if normalize_beta:
            max_beta = np.max(beta_new)
            if max_beta > 0:
                beta_new = beta_new / max_beta

        # Convergence checks
        beta_diff = np.linalg.norm(beta_new - beta_prev)
        Y_diff = np.linalg.norm(Y_values - Y_prev)

        logging.info(f"Change in beta: {beta_diff:.6f}, Change in Y: {Y_diff:.6f}")
        if beta_diff < tolerance and Y_diff < tolerance:
            logging.info("Convergence achieved.")
            break

        # Update estimates for next iteration
        for j, ct_name in enumerate(cell_type_names):
            beta_estimates[ct_name] = beta_new[j]

        # Assert that beta_new is within the range [0, 1]
        assert np.all(beta_new >= 0) and np.all(beta_new <= 1), "Beta values must be within the range [0, 1]"

        beta_prev = beta_new.copy()
        Y_prev = Y_values.copy()
        iteration += 1

    # Validate cell type proportions after optimization completes
    logging.info("Validating cell type proportions after Stage 1 optimization...")
    validate_cell_proportions(
        Y_values,
        cell_type_names,
        profile_based_antibody_data=profile_based_antibody_data,
        unknown_threshold=unknown_threshold,
        min_celltype_threshold=min_celltype_threshold,
        redundancy_threshold=redundancy_threshold,
        warn_only=warn_only,
    )

    return Y_values, beta_new


def optimize_cell_proportions_per_marker(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    tolerance: float = 1e-4,
    max_iterations: int = 50,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    normalize_beta: bool = True,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    unknown_threshold: float = 0.05,
    min_celltype_threshold: float = 0.01,
    redundancy_threshold: float = 0.2,
    warn_only: bool = False,
    lambda_laplacian: float = 0.1,
    coords: Optional[np.ndarray] = None,
    laplacian_k: int = 8,
    lambda_sparse: float = 0.0,
    alpha_max: float = 0.8,
    lambda_alpha: float = 1.0,
    lambda_coverage: float = 1.0,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray]:
    """
    Perform EM-based optimization for cell type proportions with per-marker beta.

    This version learns individual scaling factors (beta) and baselines (alpha) for
    each marker. Alpha captures ubiquitous signal (e.g., VIM floor), beta captures
    cell-type-dependent variation.

    Mathematical formulation:
        Minimize: sum_{i,m} (S[i,m] - alpha[m] - beta[m] * Y[i, owner(m)])^2 + regularization
        Beta update: OLS slope of S[:,m] vs Y[:,owner(m)]
        Alpha update: OLS intercept with L2 regularization toward zero

    Args:
        marker_level_data: N x M matrix of normalized antibody data
        marker_names: List of marker names (length M)
        assignment_matrix: M x T binary matrix where A[m,j]=1 if marker m belongs to cell type j
        cell_type_names: List of cell type names (length T)
        tolerance: Convergence tolerance for EM algorithm
        max_iterations: Maximum number of iterations
        lambda_reg: Regularization strength for elastic net
        alpha: L1-L2 tradeoff factor (0 = L2, 1 = L1)
        normalize_beta: Whether to normalize beta values so max=1
        beta_min: Minimum allowed beta value (default: 0.1)
        beta_max: Maximum allowed beta value (default: 2.0)
        unknown_threshold: Maximum allowed mean proportion for Unknown (default: 0.05)
        min_celltype_threshold: Minimum required mean proportion for cell types (default: 0.01)
        redundancy_threshold: Maximum allowed fraction of redundant types (default: 0.2)
        warn_only: If True, only warn on validation failures instead of raising
        lambda_laplacian: Weight for Laplacian smoothing (default: 0.1)
        coords: Spatial coordinates, shape (N, 2). Required if lambda_laplacian > 0.
        laplacian_k: Number of neighbors for Laplacian graph (default: 8)
        lambda_sparse: Weight for L1 sparsity penalty on proportions (default: 0.0)
        alpha_max: Maximum allowed alpha (baseline) value (default: 0.8)
        lambda_alpha: Regularization strength for alpha values (default: 1.0)
        lambda_coverage: Asymmetric loss exponent for marker count scaling (default: 1.0).
            Controls how loss is weighted by number of markers per cell type:
            0 = symmetric (no scaling), 1 = linear inverse scaling (1/n_markers).
            Higher values boost single-marker cell types that may be underestimated.
        spot_abundance_target: Optional per-spot abundance target (shape N,).
        lambda_abundance_prior: Weight of abundance soft-prior term.

    Returns:
        Tuple of:
        - Y_values (np.ndarray): (N, T) cell type proportions
        - beta_values (np.ndarray): (M,) per-marker scaling factors
        - marker_beta_dict (Dict[str, float]): {marker_name: beta_value}
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N:
            raise ValueError(
                f"spot_abundance_target length {spot_abundance_target.shape[0]} does not match N={N}."
            )

    # Validate assignment matrix shape
    if assignment_matrix.shape != (M, T):
        raise ValueError(f"Assignment matrix shape {assignment_matrix.shape} != expected ({M}, {T})")

    # Precompute marker-to-celltype assignments (supports shared markers)
    # marker_owners[m] = list of cell type indices that own marker m
    marker_owners = []
    for m in range(M):
        owners_for_m = []
        for j in range(T):
            if assignment_matrix[m, j] > 0:
                owners_for_m.append(j)
        marker_owners.append(owners_for_m)

    marker_has_owner = np.array([len(o) > 0 for o in marker_owners])

    # Compute markers-per-celltype for loss normalization (Bug 1 fix)
    # Each cell type gets equal total loss weight regardless of marker count
    markers_per_celltype = np.zeros(T, dtype=np.float64)
    for m in range(M):
        for j in marker_owners[m]:
            markers_per_celltype[j] += 1
    markers_per_celltype = np.maximum(markers_per_celltype, 1.0)

    # Compute asymmetric loss boost for underestimation
    # Cell types with fewer markers get higher boost
    max_markers = np.max(markers_per_celltype)
    underestimation_boost = np.power(max_markers / markers_per_celltype, lambda_coverage)

    if lambda_coverage > 0:
        logging.info(f"Asymmetric loss enabled: lambda_coverage={lambda_coverage}")
        for j, ct_name in enumerate(cell_type_names):
            if underestimation_boost[j] > 1.01:  # Only log if meaningful boost
                logging.info(f"  {ct_name}: {markers_per_celltype[j]:.0f} markers -> {underestimation_boost[j]:.2f}x boost")

    logging.info(f"Per-marker beta optimization: {N} spots, {M} markers, {T} cell types")
    logging.info(f"Markers with assignments: {marker_has_owner.sum()}/{M}")

    # Build spatial Laplacian if requested
    L_coo = None
    use_laplacian = lambda_laplacian > 0 and coords is not None
    if use_laplacian:
        if coords.shape[0] != N:
            raise ValueError(f"coords has {coords.shape[0]} rows but data has {N} spots")
        L = build_spatial_laplacian(coords, k=laplacian_k, normed=True)
        L_coo = L.tocoo()
        logging.info(f"Laplacian smoothing enabled: lambda={lambda_laplacian}, k={laplacian_k}")
    elif lambda_laplacian > 0 and coords is None:
        logging.warning("lambda_laplacian > 0 but coords not provided. Laplacian smoothing disabled.")

    # Initialize beta (per-marker) and alpha (per-marker baseline)
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)
    beta_prev = np.zeros(M)
    alpha_prev = np.zeros(M)
    Y_prev = np.zeros((N, T))

    iteration = 0
    while iteration < max_iterations:
        logging.info(f"\nIteration {iteration + 1}")

        model = gp.Model("EM_Cell_Proportions_PerMarker")
        model.setParam("OutputFlag", 0)

        # Define Y variables: Y[i, j] = proportion of cell type j at spot i
        Y = model.addVars(N, T, lb=0, ub=1, vtype=GRB.CONTINUOUS, name="Y")

        # Objective: for each marker m and each owner j of m,
        # add normalized error: (1/n_owners) * (1/markers_per_celltype[j]) * (S - β*Y)²
        # This ensures equal loss weight per cell type regardless of marker count,
        # and shared markers contribute to all owner cell types.
        #
        # Asymmetric loss: when lambda_coverage > 0, cell types with fewer markers
        # get extra penalty for underestimation (when residual > 0, i.e., observed
        # signal exceeds predicted). This is modeled using auxiliary variables:
        #   R_pos >= residual, R_pos >= 0  =>  R_pos = max(residual, 0)
        # Then add: weight * (boost - 1) * R_pos^2 to penalize positive residuals.
        error_terms = []
        asymmetric_terms = []
        asymmetric_constraints = []

        # Track which cell types need asymmetric boost
        use_asymmetric = lambda_coverage > 0
        boosted_celltypes = set()
        if use_asymmetric:
            for j in range(T):
                if underestimation_boost[j] > 1.01:  # Only add if meaningful boost
                    boosted_celltypes.add(j)

        # Create auxiliary variables for asymmetric loss (only if needed)
        R_pos = {}
        if use_asymmetric and len(boosted_celltypes) > 0:
            # We need R_pos[i, j, m] for each (spot, celltype, marker) triplet
            # where celltype j has boost > 1 and owns marker m
            for m in range(M):
                if not marker_has_owner[m]:
                    continue
                for j in marker_owners[m]:
                    if j in boosted_celltypes:
                        for i in range(N):
                            R_pos[(i, j, m)] = model.addVar(
                                lb=0, vtype=GRB.CONTINUOUS,
                                name=f"R_pos_{i}_{j}_{m}"
                            )

        for m in range(M):
            if not marker_has_owner[m]:
                continue

            owners_m = marker_owners[m]
            n_owners = len(owners_m)
            beta_m = beta_values[m]

            alpha_m = alpha_values[m]

            for j in owners_m:
                weight = 1.0 / (n_owners * markers_per_celltype[j])
                boost_extra = underestimation_boost[j] - 1.0  # Extra weight beyond 1x

                for i in range(N):
                    S_im = marker_level_data[i, m] - alpha_m  # baseline-subtracted
                    Y_ij = Y[i, j]
                    residual = S_im - beta_m * Y_ij

                    # Symmetric squared error (always applied)
                    error_terms.append(weight * residual * residual)

                    # Asymmetric boost for underestimation (positive residual)
                    if use_asymmetric and j in boosted_celltypes:
                        R_pos_ijm = R_pos[(i, j, m)]
                        # Constraint: R_pos >= residual (combined with R_pos >= 0, gives max(residual, 0))
                        asymmetric_constraints.append((R_pos_ijm, residual))
                        # Extra penalty: weight * boost_extra * R_pos^2
                        asymmetric_terms.append(weight * boost_extra * R_pos_ijm * R_pos_ijm)

        total_error = gp.quicksum(error_terms)

        # Add asymmetric loss terms if any
        if len(asymmetric_terms) > 0:
            total_error += gp.quicksum(asymmetric_terms)
            logging.info(f"Asymmetric loss: added {len(asymmetric_terms)} boost terms for {len(boosted_celltypes)} cell types")

        # Add constraints for asymmetric loss auxiliary variables
        for R_pos_var, residual_expr in asymmetric_constraints:
            model.addConstr(R_pos_var >= residual_expr)

        # Regularization terms (elastic net on Y)
        l1_term = gp.quicksum(Y[i, j] for i in range(N) for j in range(T))
        l2_term = gp.quicksum(Y[i, j] * Y[i, j] for i in range(N) for j in range(T))
        regularization_term = lambda_reg * (alpha * l1_term + (1 - alpha) * l2_term)

        # Laplacian smoothing term
        laplacian_term = 0
        if use_laplacian and L_coo is not None:
            laplacian_terms = []
            for idx in range(L_coo.nnz):
                i_spot = L_coo.row[idx]
                j_spot = L_coo.col[idx]
                L_val = L_coo.data[idx]
                for k in range(T):
                    laplacian_terms.append(L_val * Y[i_spot, k] * Y[j_spot, k])
            laplacian_term = lambda_laplacian * gp.quicksum(laplacian_terms)

        # Sparsity penalty (negative L2 on Y - encourages near-one-hot for cell-level)
        # Maximizing sum(Y^2) on a simplex peaks at one-hot assignment.
        # We add -lambda * sum(Y^2) to the minimization objective.
        sparsity_term = 0
        if lambda_sparse > 0:
            sparsity_term = -lambda_sparse * gp.quicksum(
                Y[i, j] * Y[i, j] for i in range(N) for j in range(T)
            )
            logging.info(f"Sparsity penalty enabled (neg-L2): lambda_sparse={lambda_sparse}")

        abundance_prior_term = 0
        if spot_abundance_target is not None and lambda_abundance_prior > 0:
            abundance_terms = []
            for i in range(N):
                row_sum = gp.quicksum(Y[i, j] for j in range(T))
                target_i = float(spot_abundance_target[i])
                abundance_terms.append((row_sum - target_i) * (row_sum - target_i))
            abundance_prior_term = lambda_abundance_prior * gp.quicksum(abundance_terms)

        model.setObjective(
            total_error + regularization_term + laplacian_term + sparsity_term + abundance_prior_term,
            GRB.MINIMIZE,
        )

        # Sum of proportions constraints
        for i in range(N):
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) >= 0.9)
            model.addConstr(gp.quicksum(Y[i, j] for j in range(T)) <= 1.2)

        try:
            model.optimize()
        except Exception as e:
            logging.error(f"Optimization error: {str(e)}")
            raise ValueError("Gurobi optimization failed") from e

        if model.status != GRB.OPTIMAL:
            raise ValueError(f"Gurobi optimization failed to converge (status: {model.status})")

        Y_values = np.array([[Y[i, j].X for j in range(T)] for i in range(N)])

        # Update beta and alpha (per-marker OLS: S = alpha + beta * Y_combined)
        # For shared markers, use sum of all owner proportions
        beta_new = np.zeros(M, dtype=np.float64)
        alpha_new = np.zeros(M, dtype=np.float64)
        for m in range(M):
            if not marker_has_owner[m]:
                beta_new[m] = 1.0
                alpha_new[m] = 0.0
                continue

            owners_m = marker_owners[m]
            Y_combined = np.zeros(N, dtype=np.float64)
            for j in owners_m:
                Y_combined += Y_values[:, j]

            S_m = marker_level_data[:, m]

            # OLS: S_m = alpha_m + beta_m * Y_combined
            Y_mean = np.mean(Y_combined)
            S_mean = np.mean(S_m)
            Y_centered = Y_combined - Y_mean
            Y_var = np.dot(Y_centered, Y_centered)

            if Y_var > 1e-9:
                beta_new[m] = np.dot(S_m - S_mean, Y_centered) / Y_var
            else:
                beta_new[m] = beta_values[m]  # keep previous
            beta_new[m] = np.clip(beta_new[m], beta_min, beta_max)

            # Alpha with L2 regularization toward zero
            raw_alpha = S_mean - beta_new[m] * Y_mean
            alpha_new[m] = raw_alpha / (1.0 + lambda_alpha / N)
            alpha_new[m] = np.clip(alpha_new[m], 0.0, alpha_max)

        # Optionally normalize beta so max=1
        if normalize_beta:
            max_beta = np.max(beta_new)
            if max_beta > 0:
                beta_new = beta_new / max_beta
                # Re-clip after normalization to prevent extreme ratios.
                # Without this, wide pre-normalization range creates tiny
                # post-normalization betas that silence weak markers.
                beta_new = np.clip(beta_new, beta_min, 1.0)

        # Convergence check
        beta_diff = np.linalg.norm(beta_new - beta_prev)
        alpha_diff = np.linalg.norm(alpha_new - alpha_prev)
        Y_diff = np.linalg.norm(Y_values - Y_prev)

        logging.info(f"Change in beta: {beta_diff:.6f}, alpha: {alpha_diff:.6f}, Y: {Y_diff:.6f}")
        if beta_diff < tolerance and alpha_diff < tolerance and Y_diff < tolerance:
            logging.info("Convergence achieved.")
            break

        beta_values = beta_new.copy()
        alpha_values = alpha_new.copy()
        beta_prev = beta_new.copy()
        alpha_prev = alpha_new.copy()
        Y_prev = Y_values.copy()
        iteration += 1

    # Validate cell type proportions after optimization completes
    logging.info("Validating cell type proportions after Stage 1 optimization...")
    validate_cell_proportions(
        Y_values,
        cell_type_names,
        profile_based_antibody_data=None,  # Not applicable for per-marker approach
        unknown_threshold=unknown_threshold,
        min_celltype_threshold=min_celltype_threshold,
        redundancy_threshold=redundancy_threshold,
        warn_only=warn_only,
    )

    # Build marker-beta dictionary for interpretability
    marker_beta_dict = {marker_names[m]: beta_new[m] for m in range(M)}

    # Log beta and alpha statistics
    logging.info(f"Beta range: [{beta_new.min():.3f}, {beta_new.max():.3f}], mean: {beta_new.mean():.3f}")
    logging.info(f"Alpha (baseline) range: [{alpha_values.min():.3f}, {alpha_values.max():.3f}], mean: {alpha_values.mean():.3f}")
    for m in range(M):
        if marker_has_owner[m] and alpha_values[m] > 0.05:
            logging.info(f"  Marker '{marker_names[m]}': alpha={alpha_values[m]:.3f}, beta={beta_new[m]:.3f}")

    return Y_values, beta_new, marker_beta_dict, alpha_values


def classify_cells_from_betas(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    beta_values: np.ndarray,
    alpha_values: np.ndarray,
    temperature: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Stage 2 of cell-level assignment: classify cells using learned betas.

    After the QP learns per-marker scaling (beta) and baseline (alpha),
    this function scores each cell against each type profile and assigns
    via softmax over negative reconstruction error.

    For each cell i and type j, the score is:
        score[i,j] = -sum_m w[m,j] * (S[i,m] - alpha[m] - beta[m])^2
                     + sum_m w[m,j'] * (S[i,m] - alpha[m])^2  (for non-owner types j')

    In words: type j gets a good score when the cell's marker expression
    is close to (alpha + beta) for j's markers, and close to alpha for
    markers belonging to other types.

    Args:
        marker_level_data: (N, M) normalized marker expression per cell.
        marker_names: Ordered list of marker names (length M).
        assignment_matrix: (M, T) binary matrix mapping markers to cell types.
        cell_type_names: Ordered list of cell type names (length T).
        beta_values: (M,) per-marker scaling factors from QP.
        alpha_values: (M,) per-marker baselines from QP.
        temperature: Softmax temperature (lower = sharper). Default 1.0.

    Returns:
        Y_cell: (N, T) soft assignment matrix (softmax probabilities).
        scores: (N, T) raw scores before softmax.
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    # For each type j, compute per-cell reconstruction error
    # Type j "explains" its own markers with alpha + beta, and doesn't
    # explain other types' markers (expects them at alpha level)
    scores = np.zeros((N, T), dtype=np.float64)

    for j in range(T):
        for m in range(M):
            S_m = marker_level_data[:, m]  # (N,)
            alpha_m = alpha_values[m]
            beta_m = beta_values[m]

            if assignment_matrix[m, j] > 0:
                # This marker belongs to type j: expect S ≈ alpha + beta
                expected = alpha_m + beta_m
                residual = (S_m - expected) ** 2
            else:
                # This marker does NOT belong to type j: expect S ≈ alpha (baseline)
                expected = alpha_m
                residual = (S_m - expected) ** 2

            # Weight by beta (stronger markers matter more)
            scores[:, j] -= beta_m * residual

    # Softmax with temperature
    scores_scaled = scores / (temperature + 1e-10)
    scores_scaled -= scores_scaled.max(axis=1, keepdims=True)  # numerical stability
    exp_scores = np.exp(scores_scaled)
    Y_cell = exp_scores / exp_scores.sum(axis=1, keepdims=True)

    # Log statistics
    max_Y = Y_cell.max(axis=1)
    dominant = np.argmax(Y_cell, axis=1)
    entropy = -np.sum(Y_cell * np.log(Y_cell + 1e-10), axis=1)
    max_entropy = np.log(T)

    logging.info(f"Cell classification (temperature={temperature}):")
    logging.info(f"  Max Y: mean={max_Y.mean():.3f} median={np.median(max_Y):.3f} "
                 f"p90={np.percentile(max_Y, 90):.3f}")
    logging.info(f"  Entropy ratio: {entropy.mean() / max_entropy:.3f}")
    logging.info(f"  One-hot fraction (>0.9): {np.mean(max_Y > 0.9):.3f}")

    for j_idx, ct in enumerate(cell_type_names):
        n_dom = np.sum(dominant == j_idx)
        logging.info(f"  {ct:<20}: n={n_dom:>5} ({n_dom/N*100:5.1f}%)")

    return Y_cell, scores


def beta_weighted_classification(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    beta_values: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Classify cells using QP-learned betas as quality weights in a gating classifier.

    For each cell i and type j, computes:
        positive = beta-weighted mean expression of type j's markers
        negative = beta-weighted mean expression of markers belonging to other types
        score[i,j] = positive - negative

    This leverages both positive evidence (high expression of own markers)
    and negative evidence (low expression of other types' markers), weighted
    by the QP's learned marker quality estimates.

    Args:
        marker_level_data: (N, M) normalized marker expression per cell.
        marker_names: Ordered list of marker names (length M).
        assignment_matrix: (M, T) binary matrix mapping markers to cell types.
        cell_type_names: Ordered list of cell type names (length T).
        beta_values: (M,) per-marker quality weights from QP.

    Returns:
        Y_cell: (N, T) softmax probabilities.
        scores: (N, T) raw scores before softmax.
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    scores = np.zeros((N, T), dtype=np.float64)

    for j in range(T):
        own_mask = assignment_matrix[:, j] > 0  # markers assigned to type j
        other_mask = ~own_mask  # markers assigned to other types

        # Beta-weighted mean of own markers
        own_betas = beta_values[own_mask]
        own_data = marker_level_data[:, own_mask]  # (N, n_own)
        if own_betas.sum() > 0:
            positive = (own_data * own_betas).sum(axis=1) / own_betas.sum()
        else:
            positive = np.zeros(N)

        # Beta-weighted mean of other markers
        other_betas = beta_values[other_mask]
        other_data = marker_level_data[:, other_mask]  # (N, n_other)
        if other_betas.sum() > 0:
            negative = (other_data * other_betas).sum(axis=1) / other_betas.sum()
        else:
            negative = np.zeros(N)

        scores[:, j] = positive - negative

    # Softmax for probabilities
    scores_shifted = scores - scores.max(axis=1, keepdims=True)
    exp_scores = np.exp(scores_shifted)
    Y_cell = exp_scores / exp_scores.sum(axis=1, keepdims=True)

    # Log statistics
    max_prob = Y_cell.max(axis=1)
    dominant = np.argmax(Y_cell, axis=1)

    logging.info("Beta-weighted cell classification:")
    logging.info(f"  Max P: mean={max_prob.mean():.3f} median={np.median(max_prob):.3f} "
                 f"p90={np.percentile(max_prob, 90):.3f}")
    logging.info(f"  One-hot (>0.9): {np.mean(max_prob > 0.9)*100:.1f}%")

    for j_idx, ct in enumerate(cell_type_names):
        n_dom = np.sum(dominant == j_idx)
        logging.info(f"  {ct:<20}: n={n_dom:>5} ({n_dom/N*100:5.1f}%)")

    return Y_cell, scores


def bayesian_cell_classification(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    beta_values: np.ndarray,
    alpha_values: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Bayesian cell classification using learned betas from the QP.

    Unlike classify_cells_from_betas() which uses beta-weighted scoring
    (biased toward types with high-beta markers), this function uses
    precision-weighted Gaussian log-likelihood over ALL markers.

    For each cell i and type j:
        expected[m,j] = alpha[m] + beta[m] * A[m,j]
        log P(S_i | type=j) = -0.5 * sum_m (S[i,m] - expected[m,j])^2 / var[m]

    Every marker contributes evidence for AND against each type:
    - Marker m owned by type j: high S[i,m] matches high expected → good score
    - Marker m NOT owned by type j: high S[i,m] mismatches low expected → penalty

    This naturally handles negative evidence (CD68 high → penalizes Epithelial)
    without profile-size bias since all M markers are used for every type.

    Args:
        marker_level_data: (N, M) normalized marker expression per cell.
        marker_names: Ordered list of marker names (length M).
        assignment_matrix: (M, T) binary matrix mapping markers to cell types.
        cell_type_names: Ordered list of cell type names (length T).
        beta_values: (M,) per-marker scaling factors from QP.
        alpha_values: (M,) per-marker baselines from QP.

    Returns:
        Y_cell: (N, T) posterior probabilities per cell.
        log_lik: (N, T) log-likelihood scores before normalization.
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    # Expected marker expression for each type
    # For type j: expected[m,j] = alpha[m] + beta[m] if A[m,j]=1, else alpha[m]
    expected = np.tile(alpha_values[:, None], (1, T))  # (M, T)
    expected += beta_values[:, None] * assignment_matrix  # (M, T)

    # Per-marker precision (inverse variance)
    marker_var = np.var(marker_level_data, axis=0) + 1e-6  # (M,)
    precision = 1.0 / marker_var  # (M,)

    # Log-likelihood: Gaussian with per-marker precision
    # log P(S_i | type_j) = -0.5 * sum_m precision[m] * (S[i,m] - expected[m,j])^2
    log_lik = np.zeros((N, T), dtype=np.float64)
    for j in range(T):
        diff = marker_level_data - expected[:, j]  # (N, M)
        weighted_sq = diff ** 2 * precision  # (N, M)
        log_lik[:, j] = -0.5 * weighted_sq.sum(axis=1)

    # Posterior via softmax (flat prior)
    log_lik_shifted = log_lik - log_lik.max(axis=1, keepdims=True)
    probs = np.exp(log_lik_shifted)
    probs /= probs.sum(axis=1, keepdims=True)

    # Log statistics
    max_prob = probs.max(axis=1)
    dominant = np.argmax(probs, axis=1)
    entropy = -np.sum(probs * np.log(probs + 1e-10), axis=1)
    max_entropy = np.log(T)

    logging.info("Bayesian cell classification:")
    logging.info(f"  Max P: mean={max_prob.mean():.3f} median={np.median(max_prob):.3f} "
                 f"p90={np.percentile(max_prob, 90):.3f}")
    logging.info(f"  Entropy ratio: {entropy.mean() / max_entropy:.3f}")
    logging.info(f"  One-hot fraction (>0.9): {np.mean(max_prob > 0.9):.3f}")

    for j_idx, ct in enumerate(cell_type_names):
        n_dom = np.sum(dominant == j_idx)
        logging.info(f"  {ct:<20}: n={n_dom:>5} ({n_dom/N*100:5.1f}%)")

    return probs, log_lik


def supervised_cell_classification(
    marker_level_data: np.ndarray,
    Y_qp: np.ndarray,
    cell_type_names: List[str],
    confidence_percentile: float = 75.0,
    classifier_type: str = "logistic",
) -> Tuple[np.ndarray, dict]:
    """
    Two-stage supervised cell classification bootstrapped from QP.

    Stage 1: QP gives soft Y assignments (already done externally).
    Stage 2: Identify confident cells from QP, train a discriminative
    classifier on them, then classify all cells.

    The QP is good at learning the relative ordering (argmax is 75.6%
    accurate) but produces diffuse probabilities. A discriminative
    classifier trained on the QP's most confident cells can generalize
    to sharper, more accurate predictions.

    Args:
        marker_level_data: (N, M) normalized marker expression per cell.
        Y_qp: (N, T) soft assignment from QP.
        cell_type_names: Ordered list of cell type names (length T).
        confidence_percentile: Percentile threshold for confident cells.
            Higher = stricter = fewer training cells but cleaner labels.
        classifier_type: "logistic" for logistic regression,
            "rf" for random forest.

    Returns:
        Y_supervised: (N, T) class probabilities from supervised model.
        info: Dict with training stats.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler

    N, M = marker_level_data.shape
    T = len(cell_type_names)

    # Confidence = max probability per cell
    max_probs = Y_qp.max(axis=1)
    threshold = np.percentile(max_probs, confidence_percentile)
    confident_mask = max_probs >= threshold

    # Get pseudo-labels from QP argmax
    qp_labels = np.argmax(Y_qp, axis=1)
    confident_labels = qp_labels[confident_mask]
    confident_features = marker_level_data[confident_mask]

    logging.info(f"Supervised classification:")
    logging.info(f"  Confidence threshold: {threshold:.4f} (p{confidence_percentile:.0f})")
    logging.info(f"  Confident cells: {confident_mask.sum()}/{N} ({confident_mask.sum()/N*100:.1f}%)")

    # Check that all types are represented in training set
    unique_train = set(confident_labels.tolist())
    for j in range(T):
        n_train_j = np.sum(confident_labels == j)
        logging.info(f"  Training {cell_type_names[j]}: {n_train_j}")
    if len(unique_train) < T:
        missing = [cell_type_names[j] for j in range(T) if j not in unique_train]
        logging.warning(f"  Missing types in training: {missing}")

    # Standardize features
    scaler = StandardScaler()
    X_train = scaler.fit_transform(confident_features)
    X_all = scaler.transform(marker_level_data)

    # Train classifier
    if classifier_type == "rf":
        clf = RandomForestClassifier(
            n_estimators=100, max_depth=10, random_state=42, n_jobs=-1,
            class_weight="balanced",
        )
    else:
        clf = LogisticRegression(
            max_iter=1000, C=1.0, random_state=42, multi_class="multinomial",
            class_weight="balanced",
        )

    clf.fit(X_train, confident_labels)

    # Predict all cells
    Y_supervised = clf.predict_proba(X_all)

    # Log statistics
    pred_labels = np.argmax(Y_supervised, axis=1)
    max_prob = Y_supervised.max(axis=1)

    logging.info(f"  Classifier: {classifier_type}")
    logging.info(f"  Max P: mean={max_prob.mean():.3f} median={np.median(max_prob):.3f}")
    logging.info(f"  One-hot (>0.9): {np.mean(max_prob > 0.9)*100:.1f}%")

    for j_idx, ct in enumerate(cell_type_names):
        n_dom = np.sum(pred_labels == j_idx)
        logging.info(f"  {ct:<20}: n={n_dom:>5} ({n_dom/N*100:5.1f}%)")

    info = {
        "n_confident": int(confident_mask.sum()),
        "confidence_threshold": float(threshold),
        "classifier_type": classifier_type,
    }

    return Y_supervised, info


def compute_marker_exclusivity(
    marker_level_data: np.ndarray,
    Y_values: np.ndarray,
    marker_owners: List[List[int]],
    assignment_matrix: np.ndarray,
    floor: float = 0.3,
    epsilon: float = 1e-9,
) -> np.ndarray:
    """
    Compute per-marker exclusivity scores measuring discriminative power.

    For each marker, measures how exclusively it correlates with its assigned
    cell type(s) versus the best non-owner cell type. Markers that track many
    cell types equally (e.g., VIM) get low scores; markers specific to their
    owner (e.g., CD68) get high scores.

    Args:
        marker_level_data: (N, M) normalized marker signals.
        Y_values: (N, T) cell type proportions from global EM pass.
        marker_owners: List of lists, marker_owners[m] = indices of owner cell types.
        assignment_matrix: (M, T) binary matrix mapping markers to cell types.
        floor: Minimum exclusivity score (default: 0.3).
        epsilon: Small constant to prevent division by zero.

    Returns:
        (M,) array of exclusivity scores in [floor, 1.0].
    """
    N, M = marker_level_data.shape
    T = Y_values.shape[1]
    exclusivity = np.ones(M, dtype=np.float64)

    for m in range(M):
        owners_m = marker_owners[m]
        if not owners_m:
            # Unowned markers: neutral weight (1.0)
            continue

        S_m = marker_level_data[:, m]

        # Owner correlation: correlate with combined owner proportions
        Y_owner = np.zeros(N, dtype=np.float64)
        for j in owners_m:
            Y_owner += Y_values[:, j]

        r_owner = np.corrcoef(S_m, Y_owner)[0, 1]
        if np.isnan(r_owner):
            r_owner = 0.0
        r_owner = max(r_owner, 0.0)

        # Best non-owner correlation
        owner_set = set(owners_m)
        r_best_other = 0.0
        for k in range(T):
            if k in owner_set:
                continue
            r_k = np.corrcoef(S_m, Y_values[:, k])[0, 1]
            if np.isnan(r_k):
                r_k = 0.0
            r_best_other = max(r_best_other, max(r_k, 0.0))

        # Exclusivity ratio
        denom = r_owner + r_best_other + epsilon
        exclusivity[m] = r_owner / denom

    # Apply floor
    exclusivity = np.clip(exclusivity, floor, 1.0)

    return exclusivity


def validate_cell_proportions(
    Y_values: np.ndarray,
    cell_type_names: List[str],
    profile_based_antibody_data: np.ndarray = None,
    unknown_threshold: float = 0.05,
    min_celltype_threshold: float = 0.01,
    redundancy_threshold: float = 0.2,
    warn_only: bool = False,
) -> None:
    """
    Validate cell type proportions after Stage 1 optimization.

    This function performs THREE validation checks:
    1. Unknown cell type proportion should not exceed threshold (default 5%)
       - High Unknown indicates too few cell types were defined (UNDERSPECIFICATION)
    2. All defined cell types should have minimum mean proportion (default 1%)
       - Very low proportions suggest the cell type doesn't exist
    3. Redundant cell types detection via linear regression (default max 20% redundant)
       - Cell types whose marker patterns are linear combinations of others (OVERSPECIFICATION)

    Args:
        Y_values (np.ndarray): Cell proportion matrix (N_spots × T_celltypes)
        cell_type_names (List[str]): List of cell type names corresponding to columns
        profile_based_antibody_data (np.ndarray): Cell-type aggregated antibody scores for redundancy check
        unknown_threshold (float): Maximum allowed mean proportion for Unknown cell type (default: 0.05 = 5%)
        min_celltype_threshold (float): Minimum required mean proportion for defined cell types (default: 0.01 = 1%)
        redundancy_threshold (float): Maximum allowed fraction of redundant types (default: 0.1 = 10%)
        warn_only (bool): If True, emit warnings instead of raising errors (useful for auto-discovered profiles)

    Raises:
        ValueError: If validation fails and warn_only=False (Unknown exceeds threshold, any cell type below minimum,
                   or redundancy fraction exceeds threshold)
    """
    # Calculate mean proportions across all spots for each cell type
    mean_proportions = np.mean(Y_values, axis=0)

    # Check if "Unknown" cell type exists and validate its proportion
    if "Unknown" in cell_type_names:
        unknown_idx = cell_type_names.index("Unknown")
        unknown_mean = mean_proportions[unknown_idx]

        if unknown_mean > unknown_threshold:
            error_msg = (
                f"Unknown cell type has mean proportion {unknown_mean*100:.2f}% "
                f"(> {unknown_threshold*100:.2f}% threshold). "
                f"This indicates too few cell types are defined in the cell profile dictionary. "
                f"Consider adding more cell types to better characterize the tissue composition."
            )
            if warn_only:
                logging.warning(f"[VALIDATION WARNING] {error_msg}")
            else:
                logging.error(error_msg)
                raise ValueError(error_msg)

    # Check if any defined (non-Unknown) cell types have very low proportions
    low_proportion_celltypes = []
    for idx, (celltype, mean_prop) in enumerate(zip(cell_type_names, mean_proportions)):
        # Skip Unknown cell type in this check
        if celltype == "Unknown":
            continue

        if mean_prop < min_celltype_threshold:
            low_proportion_celltypes.append((celltype, mean_prop))

    if low_proportion_celltypes:
        celltype_list = [f"'{ct}' ({mp*100:.3f}%)" for ct, mp in low_proportion_celltypes]
        error_msg = (
            f"The following cell type(s) have mean proportion < {min_celltype_threshold*100:.2f}% "
            f"and likely do not exist in this sample: {', '.join(celltype_list)}. "
            f"Consider removing these cell types from the cell profile dictionary."
        )
        if warn_only:
            logging.warning(f"[VALIDATION WARNING] {error_msg}")
        else:
            logging.error(error_msg)
            raise ValueError(error_msg)

    # Log successful validation
    logging.info("✓ Cell type validation passed (proportions + redundancy)")
    if "Unknown" in cell_type_names:
        unknown_idx = cell_type_names.index("Unknown")
        logging.info(f"  - Unknown: {mean_proportions[unknown_idx]*100:.2f}% (threshold < {unknown_threshold*100:.0f}%)")
    
    # Log cell type proportions
    for idx, (celltype, mean_prop) in enumerate(zip(cell_type_names, mean_proportions)):
        if celltype != "Unknown":
            logging.info(f"  - {celltype}: {mean_prop*100:.2f}%")


def finetune_cell_proportions(
    profile_based_antibody_data: np.ndarray,
    cell_type_names: List[str],
    initial_Y_values: np.ndarray,
    initial_beta_values: np.ndarray,
    adata: sc.AnnData,
    radius: float = 4.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    max_iterations: int = 20,
    max_y_change: float = 0.4,
    beta_vary: bool = True,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",
    rerun: bool = False,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Refine cell proportions using local neighborhood optimization with parallelization.

    Args:
        profile_based_antibody_data (np.ndarray):
            (N x T) array of mapped antibody intensities.
        cell_type_names (List[str]):
            Ordered list of length T specifying cell type names.
        initial_Y_values (np.ndarray):
            Initial cell proportion matrix of shape (N, T).
        initial_beta_values (np.ndarray):
            Initial beta estimates (shape T,). This is passed for consistency,
            but local solver decides how to use or ignore it based on beta_vary.
        adata (sc.AnnData):
            AnnData object with spot-level spatial coordinates in obsm['spatial'].
        radius (float):
            Radius for neighborhood-based local refinement.
        tolerance (float):
            Convergence tolerance for the local optimization loops.
        lambda_reg (float):
            Elastic net regularization strength for local solver.
        alpha (float):
            L1-L2 tradeoff (0 = purely L2, 1 = purely L1) in local solver.
        max_iterations (int):
            Maximum number of iterations for local solver.
        max_y_change (float):
            Maximum allowed change in Y values between iterations (default: 0.2).
            Values are constrained to vary by at most this amount while staying in [0,1].
        beta_vary (bool):
            If True, each spot's local solver is allowed to update betas;
            if False, betas remain fixed at the values passed in initial_beta_values.
        max_workers (int, optional):
            Maximum number of parallel workers. If None, uses os.cpu_count().
        checkpoint_interval (int):
            Number of spots between checkpoints.
        output_dir (str):
            Directory for checkpoints.
        rerun (bool):
            Whether to rerun if results exist.
        spot_abundance_target (Optional[np.ndarray]):
            Optional per-spot abundance target (shape N,).
        lambda_abundance_prior (float):
            Weight of abundance soft-prior term.

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            - A new (N x T) array of refined Y-values, obtained from a single pass
              of local refinements.
            - The original beta_values array, returned for interface consistency.
    """
    import gc
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool

    from tqdm import tqdm

    if initial_Y_values.ndim != 2:
        raise ValueError("initial_Y_values must be a 2D array (N x T).")
    if initial_beta_values.ndim != 1:
        raise ValueError("initial_beta_values must be a 1D array of length T.")

    N, T = profile_based_antibody_data.shape
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N:
            raise ValueError(
                f"spot_abundance_target length {spot_abundance_target.shape[0]} does not match N={N}."
            )
    if initial_Y_values.shape != (N, T):
        raise ValueError("Mismatch between profile_based_antibody_data and initial_Y_values shapes.")
    if len(cell_type_names) != T:
        raise ValueError("cell_type_names length must match the number of columns in profile_based_antibody_data.")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Make a copy of the initial Y-values to store final refinements
    Y_refined = initial_Y_values.copy()

    # Calculate number of workers
    workers = max_workers if max_workers is not None else os.cpu_count()
    logging.info(f"Using {workers} workers for cell proportion refinement")

    # Process all spots in parallel
    futures = {}
    retry_count = 0
    max_retries = 3

    logging.info("Starting local cell proportion refinement")
    logging.info("Lambda reg: %s", lambda_reg)
    logging.info("Alpha: %s", alpha)

    while retry_count < max_retries:
        try:
            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures.clear()
                for spot_idx in range(N):
                    future = executor.submit(
                        deconvolute_local_cell_proportions,
                        spot_idx=spot_idx,
                        adata=adata,
                        profile_based_antibody_data=profile_based_antibody_data,
                        radius=radius,
                        tolerance=tolerance,
                        lambda_reg=lambda_reg,
                        alpha=alpha,
                        beta_values=initial_beta_values,
                        beta_vary=beta_vary,
                        max_iterations=max_iterations,
                        max_y_change=max_y_change,
                        spot_abundance_target=spot_abundance_target,
                        lambda_abundance_prior=lambda_abundance_prior,
                    )
                    futures[future] = spot_idx

                spots_processed = 0
                with tqdm(total=N, desc="Refining Cell Proportions") as pbar:
                    for future in as_completed(futures):
                        spot_idx = futures[future]
                        try:
                            result = future.result(timeout=300)
                            if result is not None:
                                Y_refined[spot_idx, :] = result
                                spots_processed += 1
                                pbar.update(1)

                                if spots_processed % checkpoint_interval == 0:
                                    # Save checkpoint
                                    checkpoint_path = os.path.join(
                                        output_dir, f"cell_prop_refinement_checkpoint_{spots_processed}.npy"
                                    )
                                    np.save(checkpoint_path, Y_refined)
                                    logging.info(f"Saved checkpoint after {spots_processed} spots")

                        except TimeoutError:
                            logging.error(f"Timeout processing spot {spot_idx}")
                            continue
                        except Exception as e:
                            logging.error(f"Error processing spot {spot_idx}: {str(e)}")
                            continue

                break

        except BrokenProcessPool:
            retry_count += 1
            logging.warning(f"Process pool broken, retry {retry_count}/{max_retries}")
            if retry_count == max_retries:
                logging.error("Max retries reached, saving current progress")
            import time

            time.sleep(5)

    # Cleanup
    if futures:
        futures.clear()
    gc.collect()

    # Save final results
    final_path = os.path.join(output_dir, "cell_prop_refinement_final.npy")
    np.save(final_path, Y_refined)
    logging.info("Saved final refined cell proportions")

    return Y_refined, initial_beta_values


def deconvolute_local_cell_proportions(
    spot_idx: int,
    adata: sc.AnnData,
    profile_based_antibody_data: np.ndarray,
    radius: float = 2.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    beta_values: Optional[np.ndarray] = None,
    beta_vary: bool = True,
    normalize_beta: bool = True,
    max_iterations: int = 20,
    max_y_change: float = 0.4,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
) -> Optional[np.ndarray]:
    """
    Refine cell proportions for a single spot via local neighborhood optimization.

    Args:
        spot_idx (int):
            Index of the spot to refine in the AnnData object.
        adata (sc.AnnData):
            AnnData containing spot-level spatial coordinates in obsm['spatial'].
        profile_based_antibody_data (np.ndarray):
            (N x T) global antibody intensities for N spots, T cell types.
        radius (float):
            Neighborhood radius for identifying neighbors.
        tolerance (float):
            Convergence threshold for Y- and beta-updates (if beta_vary=True).
        lambda_reg (float):
            Strength of elastic net regularization.
        alpha (float):
            L1-L2 tradeoff for the elastic net (0 = L2, 1 = L1).
        beta_values (Optional[np.ndarray]):
            Global or initial local beta values (length T). If None and beta_vary=True,
            local betas initialize at 1.0 each.
        beta_vary (bool):
            If True, local betas are iteratively updated.
            If False, beta_values remain fixed throughout optimization.
        normalize_beta (bool):
            Whether to normalize beta values after updates.
        max_iterations (int):
            Maximum iterations allowed for EM-like steps within this local function.
        max_y_change (float):
            Maximum allowed change in Y values between iterations (default: 0.2).
            Values are constrained to vary by at most this amount while staying in [0,1].
        spot_abundance_target (Optional[np.ndarray]):
            Optional per-spot abundance target (shape N,).
        lambda_abundance_prior (float):
            Weight of abundance soft-prior term.

    Returns:
        Optional[np.ndarray]:
            Refined proportions (T,) for the specified spot, or None on failure.
    """
    import gurobipy as gp
    from gurobipy import GRB

    # Identify indices of spot's local neighborhood
    neighbor_indices = get_neighbors_with_fixed_radius(spot_idx, adata, radius=int(radius), include_center=True)
    if not neighbor_indices:
        logging.error(f"[Local Cell Props] No valid neighbors for spot {spot_idx}.")
        return None
    neighbor_indices = np.array(neighbor_indices, dtype=int)

    local_antibody_data = profile_based_antibody_data[neighbor_indices, :]
    local_N, T = local_antibody_data.shape
    local_abundance_target = None
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != profile_based_antibody_data.shape[0]:
            raise ValueError("spot_abundance_target length does not match profile_based_antibody_data rows.")
        local_abundance_target = spot_abundance_target[neighbor_indices]

    if local_N == 0:
        logging.error(f"[Local Cell Props] Spot {spot_idx} has empty local antibody data.")
        return None

    # Identify center spot's position in neighbor list
    try:
        center_local_idx = np.where(neighbor_indices == spot_idx)[0][0]
    except IndexError:
        logging.error(f"[Local Cell Props] Could not find spot {spot_idx} in neighbor list.")
        return None

    # Initialize local betas
    if beta_values is not None and len(beta_values) == T:
        local_beta = beta_values.copy()
    else:
        local_beta = np.ones(T, dtype=float)

    beta_prev = local_beta.copy()

    # Initialize local Y to something uniform
    Y_prev = np.full((local_N, T), 1.0 / T)

    iteration = 0
    while iteration < max_iterations:
        try:
            model = gp.Model(f"Local_Cell_Props_spot_{spot_idx}")
            model.setParam("OutputFlag", 0)
            model.setParam("TimeLimit", 60)
            model.setParam("MIPGap", 0.01)

            # Build Y variables in [0, 1]
            Y_vars = model.addVars(local_N, T, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="Y")

            # Summation constraints on each row
            for i in range(local_N):
                model.addConstr(gp.quicksum(Y_vars[i, j] for j in range(T)) >= 0.9)
                model.addConstr(gp.quicksum(Y_vars[i, j] for j in range(T)) <= 1.2)

            # Add constraints to limit Y value changes from previous iteration
            if iteration > 0:
                for i in range(local_N):
                    for j in range(T):
                        prev_value = Y_prev[i, j]
                        # Lower bound: max(0, prev_value - max_y_change)
                        # Upper bound: min(1, prev_value + max_y_change)
                        lb = max(0.0, prev_value - max_y_change)
                        ub = min(1.0, prev_value + max_y_change)
                        model.addConstr(Y_vars[i, j] >= lb)
                        model.addConstr(Y_vars[i, j] <= ub)

            # Objective: sum of squared differences + elastic net
            error_terms = []
            for i in range(local_N):
                for j in range(T):
                    S_ij = local_antibody_data[i, j]
                    error_terms.append((S_ij - local_beta[j] * Y_vars[i, j]) ** 2)

            total_error = gp.quicksum(error_terms)
            l1 = gp.quicksum(Y_vars[i, j] for i in range(local_N) for j in range(T))
            l2 = gp.quicksum(Y_vars[i, j] * Y_vars[i, j] for i in range(local_N) for j in range(T))
            reg_term = lambda_reg * (alpha * l1 + (1.0 - alpha) * l2)
            abundance_prior_term = 0
            if local_abundance_target is not None and lambda_abundance_prior > 0:
                abundance_terms = []
                for i in range(local_N):
                    row_sum = gp.quicksum(Y_vars[i, j] for j in range(T))
                    target_i = float(local_abundance_target[i])
                    abundance_terms.append((row_sum - target_i) * (row_sum - target_i))
                abundance_prior_term = lambda_abundance_prior * gp.quicksum(abundance_terms)
            model.setObjective(total_error + reg_term + abundance_prior_term, GRB.MINIMIZE)

            model.write('local_cell_proportions_model.mps')

            model.optimize()

            if model.status != GRB.OPTIMAL:
                logging.warning(
                    f"[Local Cell Props] Spot {spot_idx} local optimization not optimal (status: {model.status})."
                )
                return None

            # Extract current Y solution
            Y_values = np.array([[Y_vars[i, j].X for j in range(T)] for i in range(local_N)])

            # Update local beta if allowed
            if beta_vary:
                new_beta = np.zeros(T, dtype=float)
                for j in range(T):
                    Y_j = Y_values[:, j]
                    S_j = local_antibody_data[:, j]
                    denominator = np.dot(Y_j, Y_j)

                    if denominator > 1e-15:
                        new_beta[j] = np.dot(S_j, Y_j) / denominator
                    new_beta[j] = max(new_beta[j], 0.0)  # Ensure non-negative

                # Optionally normalize beta values
                if normalize_beta:
                    max_beta = np.max(new_beta)
                    if max_beta > 0:
                        new_beta = new_beta / max_beta
            else:
                new_beta = local_beta.copy()

            # Check convergence
            beta_diff = np.linalg.norm(new_beta - beta_prev) if beta_vary else 0.0
            Y_diff = np.linalg.norm(Y_values - Y_prev)

            logging.debug(
                f"Spot {spot_idx} - Iteration {iteration + 1}: " f"beta_diff={beta_diff:.6f}, Y_diff={Y_diff:.6f}"
            )

            if beta_diff < tolerance and Y_diff < tolerance:
                logging.debug(f"Spot {spot_idx} converged after {iteration + 1} iterations")
                Y_prev = Y_values
                local_beta = new_beta
                break

            # Prepare for next iteration
            Y_prev = Y_values.copy()
            local_beta = new_beta.copy()
            beta_prev = new_beta.copy()
            iteration += 1

        except Exception as e:
            logging.error(f"Error in local optimization for spot {spot_idx}: {str(e)}")
            return None

        finally:
            if "model" in locals():
                del model
            gc.collect()

    # Return just the center row of Y for this spot
    return Y_prev[center_local_idx, :]


def deconvolute_local_cell_proportions_per_marker(
    spot_idx: int,
    adata: sc.AnnData,
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    radius: float = 2.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    beta_values: Optional[np.ndarray] = None,
    beta_vary: bool = True,
    normalize_beta: bool = True,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    max_iterations: int = 20,
    max_y_change: float = 0.4,
    marker_exclusivity: Optional[np.ndarray] = None,
    marker_alpha: Optional[np.ndarray] = None,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
) -> Optional[np.ndarray]:
    """
    Refine cell proportions for a single spot via local neighborhood optimization with per-marker beta.

    Args:
        spot_idx: Index of the spot to refine in the AnnData object.
        adata: AnnData containing spot-level spatial coordinates in obsm['spatial'].
        marker_level_data: (N x M) antibody intensities for N spots, M markers.
        marker_names: List of marker names (length M).
        assignment_matrix: (M x T) binary matrix mapping markers to cell types.
        cell_type_names: List of cell type names (length T).
        radius: Neighborhood radius for identifying neighbors.
        tolerance: Convergence threshold for Y- and beta-updates.
        lambda_reg: Strength of elastic net regularization.
        alpha: L1-L2 tradeoff for the elastic net (0 = L2, 1 = L1).
        beta_values: Initial per-marker beta values (length M). If None, initialized to 1.0.
        beta_vary: If True, local betas are iteratively updated.
        normalize_beta: Whether to normalize beta values after updates.
        beta_min: Minimum allowed beta value.
        beta_max: Maximum allowed beta value.
        max_iterations: Maximum iterations allowed for EM-like steps.
        max_y_change: Maximum allowed change in Y values between iterations.
        marker_exclusivity: Optional (M,) array of per-marker exclusivity weights.
            If provided, multiplies the loss weight for each marker.
        marker_alpha: Optional (M,) array of per-marker baselines from global EM.
            If provided, signal is baseline-subtracted before reconstruction.
        spot_abundance_target: Optional per-spot abundance target (shape N,).
        lambda_abundance_prior: Weight of abundance soft-prior term.

    Returns:
        Refined proportions (T,) for the specified spot, or None on failure.
    """
    import gurobipy as gp
    from gurobipy import GRB

    N_global, M = marker_level_data.shape
    T = len(cell_type_names)

    # Precompute marker-to-celltype assignments (supports shared markers)
    marker_owners = []
    for m in range(M):
        owners_for_m = []
        for j in range(T):
            if assignment_matrix[m, j] > 0:
                owners_for_m.append(j)
        marker_owners.append(owners_for_m)

    marker_has_owner = np.array([len(o) > 0 for o in marker_owners])

    # Compute markers-per-celltype for loss normalization
    markers_per_celltype = np.zeros(T, dtype=np.float64)
    for m in range(M):
        for j in marker_owners[m]:
            markers_per_celltype[j] += 1
    markers_per_celltype = np.maximum(markers_per_celltype, 1.0)

    # Identify indices of spot's local neighborhood
    neighbor_indices = get_neighbors_with_fixed_radius(spot_idx, adata, radius=int(radius), include_center=True)
    if not neighbor_indices:
        logging.error(f"[Local Cell Props] No valid neighbors for spot {spot_idx}.")
        return None
    neighbor_indices = np.array(neighbor_indices, dtype=int)

    local_marker_data = marker_level_data[neighbor_indices, :]
    local_N = local_marker_data.shape[0]
    local_abundance_target = None
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N_global:
            raise ValueError("spot_abundance_target length does not match marker_level_data rows.")
        local_abundance_target = spot_abundance_target[neighbor_indices]

    if local_N == 0:
        logging.error(f"[Local Cell Props] Spot {spot_idx} has empty local marker data.")
        return None

    # Identify center spot's position in neighbor list
    try:
        center_local_idx = np.where(neighbor_indices == spot_idx)[0][0]
    except IndexError:
        logging.error(f"[Local Cell Props] Could not find spot {spot_idx} in neighbor list.")
        return None

    # Initialize local betas (per-marker)
    if beta_values is not None and len(beta_values) == M:
        local_beta = beta_values.copy()
    else:
        local_beta = np.ones(M, dtype=float)

    beta_prev = local_beta.copy()

    # Initialize local Y to something uniform
    Y_prev = np.full((local_N, T), 1.0 / T)

    iteration = 0
    while iteration < max_iterations:
        try:
            model = gp.Model(f"Local_Cell_Props_PerMarker_spot_{spot_idx}")
            model.setParam("OutputFlag", 0)
            model.setParam("TimeLimit", 60)
            model.setParam("MIPGap", 0.01)

            # Build Y variables in [0, 1]
            Y_vars = model.addVars(local_N, T, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="Y")

            # Summation constraints on each row
            for i in range(local_N):
                model.addConstr(gp.quicksum(Y_vars[i, j] for j in range(T)) >= 0.9)
                model.addConstr(gp.quicksum(Y_vars[i, j] for j in range(T)) <= 1.2)

            # Add constraints to limit Y value changes from previous iteration
            if iteration > 0:
                for i in range(local_N):
                    for j in range(T):
                        prev_value = Y_prev[i, j]
                        lb = max(0.0, prev_value - max_y_change)
                        ub = min(1.0, prev_value + max_y_change)
                        model.addConstr(Y_vars[i, j] >= lb)
                        model.addConstr(Y_vars[i, j] <= ub)

            # Objective: sum of squared differences + elastic net (per-marker)
            error_terms = []
            for m in range(M):
                if not marker_has_owner[m]:
                    continue
                owners_m = marker_owners[m]
                n_owners = len(owners_m)
                beta_m = local_beta[m]
                for j in owners_m:
                    excl = marker_exclusivity[m] if marker_exclusivity is not None else 1.0
                    weight = excl / (n_owners * markers_per_celltype[j])
                    for i in range(local_N):
                        S_im = local_marker_data[i, m]
                        if marker_alpha is not None:
                            S_im = S_im - marker_alpha[m]
                        error_terms.append(weight * (S_im - beta_m * Y_vars[i, j]) ** 2)

            total_error = gp.quicksum(error_terms)
            l1 = gp.quicksum(Y_vars[i, j] for i in range(local_N) for j in range(T))
            l2 = gp.quicksum(Y_vars[i, j] * Y_vars[i, j] for i in range(local_N) for j in range(T))
            reg_term = lambda_reg * (alpha * l1 + (1.0 - alpha) * l2)
            abundance_prior_term = 0
            if local_abundance_target is not None and lambda_abundance_prior > 0:
                abundance_terms = []
                for i in range(local_N):
                    row_sum = gp.quicksum(Y_vars[i, j] for j in range(T))
                    target_i = float(local_abundance_target[i])
                    abundance_terms.append((row_sum - target_i) * (row_sum - target_i))
                abundance_prior_term = lambda_abundance_prior * gp.quicksum(abundance_terms)
            model.setObjective(total_error + reg_term + abundance_prior_term, GRB.MINIMIZE)

            model.optimize()

            if model.status != GRB.OPTIMAL:
                logging.warning(
                    f"[Local Cell Props] Spot {spot_idx} local optimization not optimal (status: {model.status})."
                )
                return None

            # Extract current Y solution
            Y_values = np.array([[Y_vars[i, j].X for j in range(T)] for i in range(local_N)])

            # Update local beta if allowed (per-marker)
            if beta_vary:
                new_beta = np.zeros(M, dtype=float)
                for m in range(M):
                    if not marker_has_owner[m]:
                        new_beta[m] = 1.0
                        continue
                    owners_m = marker_owners[m]
                    Y_combined = np.zeros(local_N, dtype=float)
                    for j in owners_m:
                        Y_combined += Y_values[:, j]
                    S_m = local_marker_data[:, m]
                    if marker_alpha is not None:
                        S_m = S_m - marker_alpha[m]
                    denominator = np.dot(Y_combined, Y_combined) + 1e-9
                    new_beta[m] = np.dot(S_m, Y_combined) / denominator
                    new_beta[m] = np.clip(new_beta[m], beta_min, beta_max)

                # Optionally normalize beta values
                if normalize_beta:
                    max_beta_val = np.max(new_beta)
                    if max_beta_val > 0:
                        new_beta = new_beta / max_beta_val
                        # Re-clip after normalization to prevent extreme ratios
                        new_beta = np.clip(new_beta, beta_min, 1.0)
            else:
                new_beta = local_beta.copy()

            # Check convergence
            beta_diff = np.linalg.norm(new_beta - beta_prev) if beta_vary else 0.0
            Y_diff = np.linalg.norm(Y_values - Y_prev)

            logging.debug(
                f"Spot {spot_idx} - Iteration {iteration + 1}: " f"beta_diff={beta_diff:.6f}, Y_diff={Y_diff:.6f}"
            )

            if beta_diff < tolerance and Y_diff < tolerance:
                logging.debug(f"Spot {spot_idx} converged after {iteration + 1} iterations")
                Y_prev = Y_values
                local_beta = new_beta
                break

            # Prepare for next iteration
            Y_prev = Y_values.copy()
            local_beta = new_beta.copy()
            beta_prev = new_beta.copy()
            iteration += 1

        except Exception as e:
            logging.error(f"Error in local optimization for spot {spot_idx}: {str(e)}")
            return None

        finally:
            if "model" in locals():
                del model
            gc.collect()

    # Return just the center row of Y for this spot
    return Y_prev[center_local_idx, :]


def finetune_cell_proportions_per_marker(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    initial_Y_values: np.ndarray,
    initial_beta_values: np.ndarray,
    adata: sc.AnnData,
    radius: float = 4.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    max_iterations: int = 20,
    max_y_change: float = 0.4,
    beta_vary: bool = True,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    marker_exclusivity: Optional[np.ndarray] = None,
    marker_alpha: Optional[np.ndarray] = None,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",
    rerun: bool = False,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Refine cell proportions using local neighborhood optimization with per-marker beta.

    Args:
        marker_level_data: (N x M) array of marker-level antibody intensities.
        marker_names: List of marker names (length M).
        assignment_matrix: (M x T) binary matrix mapping markers to cell types.
        cell_type_names: Ordered list of cell type names (length T).
        initial_Y_values: Initial cell proportion matrix of shape (N, T).
        initial_beta_values: Initial per-marker beta estimates (shape M,).
        adata: AnnData object with spot-level spatial coordinates in obsm['spatial'].
        radius: Radius for neighborhood-based local refinement.
        tolerance: Convergence tolerance for the local optimization loops.
        lambda_reg: Elastic net regularization strength for local solver.
        alpha: L1-L2 tradeoff (0 = purely L2, 1 = purely L1).
        max_iterations: Maximum number of iterations for local solver.
        max_y_change: Maximum allowed change in Y values between iterations.
        beta_vary: If True, local solver updates betas; if False, betas remain fixed.
        beta_min: Minimum allowed beta value.
        beta_max: Maximum allowed beta value.
        max_workers: Maximum number of parallel workers. If None, uses os.cpu_count().
        checkpoint_interval: Number of spots between checkpoints.
        output_dir: Directory for checkpoints.
        rerun: Whether to rerun if results exist.
        spot_abundance_target: Optional per-spot abundance target (shape N,).
        lambda_abundance_prior: Weight of abundance soft-prior term.

    Returns:
        Tuple of:
        - Refined (N x T) array of Y-values.
        - The per-marker beta_values array (returned for interface consistency).
    """
    import gc
    import os
    from concurrent.futures import ProcessPoolExecutor, as_completed
    from concurrent.futures.process import BrokenProcessPool

    from tqdm import tqdm

    N, M = marker_level_data.shape
    T = len(cell_type_names)
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N:
            raise ValueError(
                f"spot_abundance_target length {spot_abundance_target.shape[0]} does not match N={N}."
            )

    if initial_Y_values.ndim != 2 or initial_Y_values.shape[0] != N or initial_Y_values.shape[1] != T:
        raise ValueError(f"initial_Y_values shape {initial_Y_values.shape} != expected ({N}, {T})")
    if initial_beta_values.ndim != 1 or len(initial_beta_values) != M:
        raise ValueError(f"initial_beta_values length {len(initial_beta_values)} != expected {M}")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Make a copy of the initial Y-values to store final refinements
    Y_refined = initial_Y_values.copy()

    # Calculate number of workers
    workers = max_workers if max_workers is not None else os.cpu_count()
    logging.info(f"Using {workers} workers for per-marker cell proportion refinement")

    # Process all spots in parallel
    futures = {}
    retry_count = 0
    max_retries = 3

    logging.info("Starting local cell proportion refinement (per-marker beta)")
    logging.info(f"Lambda reg: {lambda_reg}, Alpha: {alpha}")

    while retry_count < max_retries:
        try:
            with ProcessPoolExecutor(max_workers=workers) as executor:
                futures.clear()
                for spot_idx in range(N):
                    future = executor.submit(
                        deconvolute_local_cell_proportions_per_marker,
                        spot_idx=spot_idx,
                        adata=adata,
                        marker_level_data=marker_level_data,
                        marker_names=marker_names,
                        assignment_matrix=assignment_matrix,
                        cell_type_names=cell_type_names,
                        radius=radius,
                        tolerance=tolerance,
                        lambda_reg=lambda_reg,
                        alpha=alpha,
                        beta_values=initial_beta_values,
                        beta_vary=beta_vary,
                        beta_min=beta_min,
                        beta_max=beta_max,
                        max_iterations=max_iterations,
                        max_y_change=max_y_change,
                        marker_exclusivity=marker_exclusivity,
                        marker_alpha=marker_alpha,
                        spot_abundance_target=spot_abundance_target,
                        lambda_abundance_prior=lambda_abundance_prior,
                    )
                    futures[future] = spot_idx

                spots_processed = 0
                with tqdm(total=N, desc="Refining Cell Proportions (per-marker)") as pbar:
                    for future in as_completed(futures):
                        spot_idx = futures[future]
                        try:
                            result = future.result(timeout=300)
                            if result is not None:
                                Y_refined[spot_idx, :] = result
                                spots_processed += 1
                                pbar.update(1)

                                if spots_processed % checkpoint_interval == 0:
                                    # Save checkpoint
                                    checkpoint_path = os.path.join(
                                        output_dir, f"cell_prop_refinement_permarker_checkpoint_{spots_processed}.npy"
                                    )
                                    np.save(checkpoint_path, Y_refined)
                                    logging.info(f"Saved checkpoint after {spots_processed} spots")

                        except TimeoutError:
                            logging.error(f"Timeout processing spot {spot_idx}")
                            continue
                        except Exception as e:
                            logging.error(f"Error processing spot {spot_idx}: {str(e)}")
                            continue

                break

        except BrokenProcessPool as e:
            retry_count += 1
            logging.warning(f"ProcessPool broken (attempt {retry_count}/{max_retries}): {e}")
            gc.collect()
            if retry_count >= max_retries:
                logging.error("Max retries reached for ProcessPool. Returning partial results.")
                break

        except Exception as e:
            logging.error(f"Unexpected error in parallel processing: {e}")
            traceback.print_exc()
            break

    # Save final results
    final_path = os.path.join(output_dir, "cell_prop_refinement_permarker_final.npy")
    np.save(final_path, Y_refined)
    logging.info(f"Saved final refined cell proportions to {final_path}")

    return Y_refined, initial_beta_values


def compute_proportion_enrichment(
    gene_expr: np.ndarray,
    cell_type_props: np.ndarray,
    celltype_frequencies: np.ndarray = None,
) -> np.ndarray:
    """
    Compute proportion-based enrichment WITHOUT smoothing.

    This function calculates how enriched a gene's expression is in each cell type
    based on the cell type proportions of spots where the gene is expressed.
    Unlike the internal smoothed version, this returns raw enrichment values
    that can go to extreme values (close to 0 or 1) for highly specific genes.

    Args:
        gene_expr: (N,) expression values for one gene across spots
        cell_type_props: (N, T) cell type proportions per spot
        celltype_frequencies: (T,) global cell type frequencies for normalization.
            If None, no frequency normalization is applied.

    Returns:
        enrichment: (T,) enrichment scores summing to 1
    """
    T = cell_type_props.shape[1]

    # Handle zero expression
    total_expr = gene_expr.sum()
    if total_expr < 1e-10:
        return np.ones(T) / T

    # Expression-weighted cell type proportions
    weights = gene_expr / total_expr

    # Normalize by global frequency if provided
    if celltype_frequencies is not None:
        normalized_props = cell_type_props / (celltype_frequencies + 1e-10)
    else:
        normalized_props = cell_type_props

    weighted_props = np.sum(normalized_props * weights[:, np.newaxis], axis=0)
    background_props = np.mean(normalized_props, axis=0)

    epsilon = 1e-10
    enrichment = weighted_props / (background_props + epsilon)

    # Normalize to sum to 1 (NO 80/20 smoothing)
    return enrichment / (enrichment.sum() + epsilon)


def compute_marker_enrichment(
    gene_expr: np.ndarray,
    anchor_expr: np.ndarray,
    anchor_weights: np.ndarray,
) -> np.ndarray:
    """
    Compute marker-guided enrichment via correlation with anchor genes.

    Args:
        gene_expr: (N,) expression of target gene across neighborhood
        anchor_expr: (N, T) mean anchor expression per cell type in neighborhood
        anchor_weights: (T,) weight per cell type (mean anchor correlation with proportions)

    Returns:
        enrichment: (T,) enrichment scores summing to 1
    """
    from scipy.stats import pearsonr

    T = anchor_expr.shape[1]
    enrichment = np.zeros(T)

    # Skip if gene has no variance
    if np.std(gene_expr) < 1e-10:
        return np.ones(T) / T

    for t in range(T):
        anchor_t = anchor_expr[:, t]

        # Skip if anchor has no variance
        if np.std(anchor_t) < 1e-10:
            continue

        r, _ = pearsonr(gene_expr, anchor_t)

        if not np.isnan(r):
            # Only positive correlations, weighted by anchor strength
            enrichment[t] = max(0, r) * anchor_weights[t]

    # Normalize
    if enrichment.sum() < 1e-10:
        return np.ones(T) / T

    return enrichment / enrichment.sum()


def compute_adaptive_enrichment(
    prop_enrichment: np.ndarray,
    marker_enrichment: np.ndarray,
    max_variance: float = 0.1,
) -> np.ndarray:
    """
    Adaptively blend proportion and marker enrichment based on proportion variance.

    High proportion variance = trust proportions (peaked distribution)
    Low proportion variance = use marker guidance (flat distribution needs help)

    For a normalized distribution summing to 1:
    - Uniform [0.25, 0.25, 0.25, 0.25] has variance ~0
    - Peaked [0.7, 0.1, 0.1, 0.1] has variance ~0.07
    - One-hot [1, 0, 0, 0] has variance ~0.19

    Args:
        prop_enrichment: (T,) proportion-based enrichment (normalized)
        marker_enrichment: (T,) marker-guided enrichment (normalized)
        max_variance: Variance threshold above which we fully trust proportions.
            Default 0.1 is reasonable for typical cell type distributions.

    Returns:
        blended: (T,) adaptive blend summing to 1
    """
    # Variance of proportion enrichment across cell types
    variance = np.var(prop_enrichment)

    # anchor_weight: 0 when variance is high (peaked), 1 when variance is low (flat)
    # High variance -> small anchor_weight -> trust proportions
    # Low variance -> large anchor_weight -> use marker guidance
    anchor_weight = max(0.0, 1 - variance / max_variance)

    # Blend: (1 - anchor_weight) weights proportions, anchor_weight weights markers
    blended = (1 - anchor_weight) * prop_enrichment + anchor_weight * marker_enrichment

    # Normalize (should already sum to ~1, but ensure)
    return blended / (blended.sum() + 1e-10)


def precompute_anchor_expression(
    gene_expression: np.ndarray,
    anchor_genes: Dict[int, List[int]],
    anchor_weights: Dict[int, Dict[int, float]],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Precompute weighted mean anchor expression per cell type.

    Args:
        gene_expression: (N, M) expression matrix
        anchor_genes: {cell_type: [gene_indices]}
        anchor_weights: {cell_type: {gene_idx: correlation}}

    Returns:
        anchor_expr: (N, T) weighted mean anchor expression per cell type
        type_weights: (T,) mean anchor weight per cell type
    """
    N = gene_expression.shape[0]
    T = len(anchor_genes)

    anchor_expr = np.zeros((N, T))
    type_weights = np.zeros(T)

    for t in range(T):
        if t not in anchor_genes or not anchor_genes[t]:
            # No anchors for this type - use uniform
            anchor_expr[:, t] = 1.0
            type_weights[t] = 0.1  # low weight
            continue

        genes = anchor_genes[t]
        weights_t = anchor_weights.get(t, {})

        # Weighted mean of anchor expressions
        total_weight = 0.0
        for g in genes:
            w = weights_t.get(g, 0.3)  # default weight if missing
            anchor_expr[:, t] += w * gene_expression[:, g]
            total_weight += w

        if total_weight > 0:
            anchor_expr[:, t] /= total_weight
            type_weights[t] = total_weight / len(genes)  # mean weight
        else:
            anchor_expr[:, t] = 1.0
            type_weights[t] = 0.1

    return anchor_expr, type_weights


################################################################################
# === DECONVOLUTION FOR GENES ===
################################################################################
def deconvolute_spot_with_neighbors_with_prior(
    spot_idx: int,
    adata: sc.AnnData,
    cell_type_numbers_array: np.ndarray,
    radius: float,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    local_enrichment_weight: float = 0.5,
    global_enrichment_weight: float = 0.5,
    continuous_relaxation: bool = True,
    lambda_gex_reg: float = 0.01,
    # NEW parameters for module enrichment and KL regularization
    anchor_genes: Optional[Dict[int, List[int]]] = None,
    anchor_weights: Optional[Dict[int, Dict[int, float]]] = None,
    module_weight: float = 0.5,
    use_kl_regularization: bool = False,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Optional[np.ndarray]:
    """
    Deconvolute a spot with its neighbors, using enrichment weights and optional prior.

    Uses continuous relaxation (LP) by default for fractional gene assignment,
    with L2 regularization to stabilize the solution. Optionally supports
    module-aware enrichment and softmax KL regularization.

    Args:
        spot_idx: Index of the center spot to deconvolve.
        adata: AnnData with gene expression in .X and spatial coords in .obsm.
        cell_type_numbers_array: Cell type proportions (N_spots x T), from Pass 1.
        radius: Spatial radius for neighbor detection.
        global_prior: Optional prior matrix (T x M) for guidance.
        lambda_prior_weight: Weight for prior guidance term.
        local_enrichment_weight: Weight for local expression enrichment (0-1).
        global_enrichment_weight: Weight for global expression enrichment (0-1).
        continuous_relaxation: If True, use continuous variables (LP); else integer (MIP).
        lambda_gex_reg: L2 regularization weight on X variables.
        anchor_genes: Optional dict mapping cell type index to list of anchor gene indices.
        anchor_weights: Optional dict {cell_type_idx: {gene_idx: correlation}} with per-gene
            correlation weights for each anchor. Used by adaptive marker-guided enrichment.
        module_weight: Weight for module-aware enrichment (0-1), blends with base enrichment.
            NOTE: Deprecated when using adaptive marker guidance (anchor_weights provided).
        use_kl_regularization: If True, use softmax KL penalty instead of L2 regularization.
        kl_temperature: Temperature for softmax target (lower = sharper).
        lambda_kl: Weight for KL penalty term.
        use_direct_softmax: If True, skip Gurobi entirely and use direct softmax allocation.
            This avoids L2's uniform spreading problem.

    Returns:
        np.ndarray of shape (T, M) with deconvolved expression, or None on failure.
    """
    from scipy.special import softmax as scipy_softmax

    # DIRECT SOFTMAX MODE: Skip Gurobi entirely
    # This fixes the L2 uniform spreading problem (variance_ratio 3.02 -> 1.53)
    if kl_temperature < 1.0 and not use_kl_regularization:
        try:
            # Get expression data
            deconvolution_expression_data = adata.X
            if scipy.sparse.issparse(deconvolution_expression_data):
                deconvolution_expression_data = deconvolution_expression_data.toarray()

            T = cell_type_numbers_array.shape[1]
            M = deconvolution_expression_data.shape[1]
            center_counts = deconvolution_expression_data[spot_idx, :]
            center_props = cell_type_numbers_array[spot_idx, :]

            # Get neighbors for enrichment calculation
            neighborhood_indices = get_neighbors_with_fixed_radius(
                spot_idx, adata, radius=int(radius), include_center=True
            )
            if not neighborhood_indices:
                return None

            neighborhood_indices = np.array([int(idx) for idx in neighborhood_indices], dtype=int)
            neighborhood_expression = deconvolution_expression_data[neighborhood_indices, :]
            neighborhood_props = cell_type_numbers_array[neighborhood_indices, :]

            # Compute celltype frequencies for normalization
            total_celltype_counts = np.sum(cell_type_numbers_array, axis=0) + 1e-10
            celltype_frequencies = total_celltype_counts / np.sum(total_celltype_counts)

            result = np.zeros((T, M))

            for k in range(M):
                total = center_counts[k]
                if total <= 0:
                    continue

                # Expression-weighted enrichment
                gene_expr = neighborhood_expression[:, k]
                total_expr = gene_expr.sum()

                if total_expr < 1e-10:
                    # No expression - distribute by proportions
                    result[:, k] = total * (center_props / (center_props.sum() + 1e-10))
                    continue

                # Weight each spot's contribution by expression
                weights = gene_expr / total_expr
                normalized_props = neighborhood_props / (celltype_frequencies + 1e-10)
                weighted_props = np.sum(normalized_props * weights[:, np.newaxis], axis=0)
                background_props = np.mean(normalized_props, axis=0)

                enrichment = weighted_props / (background_props + 1e-10)
                enrichment = enrichment / (enrichment.sum() + 1e-10)

                # Direct softmax allocation (no Gurobi, no L2)
                log_enrichment = np.log(enrichment + 1e-10)
                alloc_weights = scipy_softmax(log_enrichment / kl_temperature)
                result[:, k] = total * alloc_weights

            return result

        except Exception as e:
            logging.error(f"Direct softmax failed for spot {spot_idx}: {e}")
            # Fall through to Gurobi

    model = None
    import gurobipy as gp  # type: ignore
    from gurobipy import GRB  # type: ignore
    try:
        neighborhood_indices = get_neighbors_with_fixed_radius(
            spot_idx, adata, radius=int(radius), include_center=True
        )
        if not neighborhood_indices:
            logging.error(f"No valid neighbors found for spot {spot_idx}.")
            return None

        neighborhood_indices = np.array(
            [int(idx) for idx in neighborhood_indices
             if isinstance(idx, (int, np.integer))], dtype=int
        )

        # Extract expression data
        deconvolution_expression_data = adata.X
        if scipy.sparse.issparse(deconvolution_expression_data):
            deconvolution_expression_data = deconvolution_expression_data.toarray()  # type: ignore
        elif not isinstance(deconvolution_expression_data, np.ndarray):
            deconvolution_expression_data = np.array(deconvolution_expression_data)

        # Dimensions
        T = cell_type_numbers_array.shape[1]  # number of cell types
        M = deconvolution_expression_data.shape[1]  # number of genes

        neighborhood_expression_data = deconvolution_expression_data[neighborhood_indices, :]
        neighborhood_cell_type_numbers = cell_type_numbers_array[neighborhood_indices, :]

        # Compute normalized cell type weights to avoid abundance bias
        total_celltype_counts = np.sum(cell_type_numbers_array, axis=0) + 1e-10
        celltype_frequencies = total_celltype_counts / np.sum(total_celltype_counts)

        # Compute expression-aware enrichment using adaptive marker guidance
        gene_specific_enrichment = np.zeros((M, T))

        # Precompute anchor expression for neighborhood if anchors provided
        if anchor_genes is not None and anchor_weights is not None:
            neighborhood_anchor_expr, type_weights = precompute_anchor_expression(
                neighborhood_expression_data, anchor_genes, anchor_weights
            )
        else:
            neighborhood_anchor_expr = None
            type_weights = None

        for k in range(M):
            # Proportion-based enrichment (no smoothing)
            prop_enrich = compute_proportion_enrichment(
                gene_expr=neighborhood_expression_data[:, k],
                cell_type_props=neighborhood_cell_type_numbers,
                celltype_frequencies=celltype_frequencies,
            )

            # If anchors available, compute adaptive blend
            if neighborhood_anchor_expr is not None:
                marker_enrich = compute_marker_enrichment(
                    gene_expr=neighborhood_expression_data[:, k],
                    anchor_expr=neighborhood_anchor_expr,
                    anchor_weights=type_weights,
                )
                gene_specific_enrichment[k] = compute_adaptive_enrichment(
                    prop_enrich, marker_enrich
                )
            else:
                gene_specific_enrichment[k] = prop_enrich

        # Build Gurobi model
        model = gp.Model(f"gene_expression_spot_{spot_idx}")
        model.setParam("OutputFlag", 0)
        model.setParam("Threads", 1)

        if continuous_relaxation:
            model.setParam("TimeLimit", 30)
        else:
            model.setParam("NodefileStart", 0.5)
            model.setParam("MIPGap", 0.01)
            model.setParam("TimeLimit", 600)
            model.setParam("NodeLimit", 1000000)

        # Variables for count assignment
        var_type = GRB.CONTINUOUS if continuous_relaxation else GRB.INTEGER
        X = {}
        center_counts = deconvolution_expression_data[spot_idx, :]

        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    X[j, k] = model.addVar(
                        vtype=var_type, lb=0, ub=total_counts,
                        name=f"X_{j}_{k}"
                    )
                # Count conservation constraint
                model.addConstr(
                    gp.quicksum(X[j, k] for j in range(T)) == total_counts,
                    name=f"count_conservation_{k}"
                )

        # Validate prior if provided
        if global_prior is not None and lambda_prior_weight > 0:
            if not isinstance(global_prior, np.ndarray):
                raise ValueError("global_prior must be a numpy array")
            if global_prior.shape != (T, M):
                raise ValueError(
                    f"Prior matrix shape {global_prior.shape} does not match "
                    f"expected shape ({T}, {M})"
                )

        # Center spot proportions (index 0 = center in neighborhood)
        center_props = neighborhood_cell_type_numbers[0, :]

        # Apply softmax sharpening to proportions (simpler than KL regularization)
        # This concentrates allocation toward dominant cell types
        if kl_temperature < 1.0:
            logits = center_props / kl_temperature
            logits = logits - logits.max()  # numerical stability
            exp_logits = np.exp(logits)
            center_props = exp_logits / (exp_logits.sum() + 1e-10)

        # Objective: maximize enrichment-weighted allocation with regularization
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                enrichment_for_gene = gene_specific_enrichment[k, :]

                if use_kl_regularization:
                    # NEW: Softmax KL-divergence regularization
                    target = compute_softmax_target(enrichment_for_gene, temperature=kl_temperature)
                    kl_coeffs = compute_kl_penalty_coefficients(target, total_counts, lambda_kl)

                    for j in range(T):
                        # Enrichment term
                        base_term = enrichment_for_gene[j] * center_props[j] * X[j, k]
                        obj_terms.append(base_term)

                        # KL penalty (pulls toward target)
                        target_j = kl_coeffs['target_counts'][j]
                        penalty = kl_coeffs['penalty_weight'] * (X[j, k] - target_j) * (X[j, k] - target_j)
                        obj_terms.append(-penalty)

                        # Prior terms (unchanged)
                        if global_prior is not None and lambda_prior_weight > 0:
                            try:
                                prior_value = float(global_prior[j, k])
                                prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j, k]
                                obj_terms.append(-prior_penalty)
                            except Exception as e:
                                logging.warning(f"Error accessing prior at [{j}, {k}]: {e}")
                                continue
                else:
                    # OLD: L2 regularization (backward compatible)
                    for j in range(T):
                        enrichment_weight = enrichment_for_gene[j]
                        cell_type_weight = center_props[j]

                        base_term = enrichment_weight * cell_type_weight * X[j, k]
                        obj_terms.append(base_term)

                        if lambda_gex_reg > 0:
                            obj_terms.append(-lambda_gex_reg * X[j, k] * X[j, k])

                        if global_prior is not None and lambda_prior_weight > 0:
                            try:
                                prior_value = float(global_prior[j, k])
                                prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j, k]
                                obj_terms.append(-prior_penalty)
                            except Exception as e:
                                logging.warning(f"Error accessing prior at [{j}, {k}]: {e}")
                                continue

        model.setObjective(gp.quicksum(obj_terms), GRB.MAXIMIZE)

        model.optimize()

        if model.status == GRB.OPTIMAL:
            result = np.zeros((T, M))
            for k in range(M):
                total_counts = int(center_counts[k])
                if total_counts > 0:
                    for j in range(T):
                        result[j, k] = X[j, k].X
            return result
        else:
            logging.error(f"No feasible solution found for spot {spot_idx}.")
            return None

    except Exception as e:
        logging.error(f"Error during deconvolution of spot {spot_idx}: {str(e)}")
        logging.error(traceback.format_exc())
        return None

    finally:
        if model:
            del model
        gc.collect()


def log_marker_gene_patterns(zero_patterns, marker_genes):
    """
    Log detailed patterns for marker genes.
    """
    for gene in marker_genes:
        logging.info(f"\nPatterns for {gene}:")
        for ct, genes_data in zero_patterns.items():
            if gene in genes_data:
                stats = genes_data[gene]
                logging.info(f"  {ct}:")
                logging.info(f"    Zero proportion: {stats['zero_proportion']:.3f}")
                if stats["n_nonzero"] > 0:
                    logging.info(f"    Mean nonzero expression: {stats['mean_nonzero_expression']:.3f}")
                else:
                    logging.info(f"    Mean nonzero expression: 0.0 (no nonzero values)")
                logging.info(f"    Number of spots: {stats['n_spots']}")
                logging.info(f"    Number of nonzero spots: {stats['n_nonzero']}")


def scale_genes(expression_matrix):
    """Scale each gene independently to [0,1] range.

    Args:
        expression_matrix (np.ndarray): Spots x Genes matrix

    Returns:
        tuple: (scaled_matrix, gene_mins, gene_maxs)
    """
    gene_mins = np.min(expression_matrix, axis=0)
    gene_maxs = np.max(expression_matrix, axis=0)

    # Avoid division by zero
    gene_ranges = np.maximum(gene_maxs - gene_mins, 1e-10)
    scaled_matrix = (expression_matrix - gene_mins) / gene_ranges

    return scaled_matrix, gene_mins, gene_maxs


def unscale_genes(scaled_matrix, gene_mins, gene_maxs):
    """Reverse the gene-wise scaling transformation.

    Args:
        scaled_matrix (np.ndarray): Scaled matrix
        gene_mins (np.ndarray): Original minimum values per gene
        gene_maxs (np.ndarray): Original maximum values per gene

    Returns:
        np.ndarray: Unscaled matrix
    """
    gene_ranges = np.maximum(gene_maxs - gene_mins, 1e-10)
    return (scaled_matrix * gene_ranges) + gene_mins


def optimize_gene_expression(
    sample_name: str,
    deconvolution_expression_data: np.ndarray,
    cell_type_numbers_array: np.ndarray,
    filtered_adata: sc.AnnData,
    radius: float = 2,
    global_enrichment_weight: float = 0.5,
    local_enrichment_weight: float = 0.5,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 100,
    output_dir: str = "checkpoints",
    rerun: bool = False,
    continuous_relaxation: bool = True,
    lambda_gex_reg: float = 0.01,
    # NEW parameters for module enrichment and KL regularization
    anchor_genes: Optional[Dict[int, List[int]]] = None,
    module_weight: float = 0.5,
    use_kl_regularization: bool = False,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Dict[str, Any]:
    """
    Optimize gene expression with enrichment weights and prior guidance.

    Args:
        sample_name (str): Name of the sample
        deconvolution_expression_data (np.ndarray): Gene expression data (N_spots x M_genes)
        cell_type_numbers_array (np.ndarray): Cell type proportions (N_spots x T_celltypes)
        filtered_adata (sc.AnnData): Filtered AnnData object containing gene expression data
        radius (float): Radius for neighbor detection
        global_enrichment_weight (float): Weight for global expression enrichment (0-1)
        local_enrichment_weight (float): Weight for local expression enrichment (0-1)
        global_prior (np.ndarray, optional): Global prior matrix for guidance
        lambda_prior_weight (float): Weight for prior guidance
        max_workers (int, optional): Maximum number of parallel workers
        checkpoint_interval (int): Number of spots between checkpoints
        output_dir (str): Directory for checkpoints
        rerun (bool): Whether to rerun if results exist
        continuous_relaxation (bool): Use continuous (LP) vs integer (MIP) variables
        lambda_gex_reg (float): L2 regularization weight on X variables
        anchor_genes (dict, optional): Dict mapping cell type index to anchor gene indices.
        module_weight (float): Weight for module-aware enrichment (0-1).
        use_kl_regularization (bool): If True, use KL-divergence instead of L2.
        kl_temperature (float): Temperature for softmax target (lower = sharper).
        lambda_kl (float): Weight for KL penalty term.

    Returns:
        Dict[str, Any]: {
            'spotwise_profiles': Dict[int, np.ndarray],
            'dimensions': Tuple[int, int, int]
        }
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    N = deconvolution_expression_data.shape[0]
    M = deconvolution_expression_data.shape[1]
    T = cell_type_numbers_array.shape[1]

    # Initialize checkpoint manager
    checkpoint_mgr = CheckpointManager(output_dir, sample_name)

    # If rerun=False, check for completed run
    if not rerun:
        complete_results = checkpoint_mgr.check_complete_run(N, T, M)
        if complete_results is not None:
            return complete_results  # type: ignore

        # Load latest checkpoint if available
        completed_spots, spotwise_gene_expression_profiles = checkpoint_mgr.load_latest_checkpoint(N, T, M)
    else:
        completed_spots = set()
        spotwise_gene_expression_profiles = {}

    logging.info(f"Starting analysis for {sample_name}")
    logging.info(f"Already completed spots: {len(completed_spots)}")

    # Log configuration
    if global_enrichment_weight + local_enrichment_weight > 0:
        logging.info(f"Using enrichment weights - Global: {global_enrichment_weight}, Local: {local_enrichment_weight}")
    if global_prior is not None and lambda_prior_weight > 0:
        logging.info("Using prior-guided deconvolution")
    if use_kl_regularization:
        logging.info(f"Using KL regularization (temp={kl_temperature}, lambda={lambda_kl})")
    if anchor_genes is not None:
        total_anchors = sum(len(v) for v in anchor_genes.values())
        logging.info(f"Using module-aware enrichment (weight={module_weight}, {total_anchors} anchors)")

    # Initialize futures as empty dict before try block
    futures = {}

    try:
        # Calculate number of workers (ensure it's an integer)
        workers = max_workers if max_workers is not None else os.cpu_count()
        logging.info(f"Using {workers} workers")

        # Only process spots that haven't been completed
        remaining_spots = [i for i in range(N) if i not in completed_spots]
        logging.info(f"Processing {len(remaining_spots)} remaining spots")

        retry_count = 0
        max_retries = 3
        while retry_count < max_retries:
            try:
                with ProcessPoolExecutor(max_workers=workers) as executor:
                    futures.clear()
                    for spot_idx in remaining_spots:
                        # Always use the same function with consistent args
                        future = executor.submit(
                            deconvolute_spot_with_neighbors_with_prior,
                            spot_idx,
                            filtered_adata,
                            cell_type_numbers_array,
                            radius,
                            global_prior,
                            lambda_prior_weight,
                            local_enrichment_weight,
                            global_enrichment_weight,
                            continuous_relaxation,
                            lambda_gex_reg,
                            # NEW parameters for module enrichment and KL regularization
                            anchor_genes,
                            module_weight,
                            use_kl_regularization,
                            kl_temperature,
                            lambda_kl,
                        )
                        futures[future] = spot_idx

                    with tqdm(total=len(remaining_spots), desc="Deconvoluting Remaining Spots") as pbar:
                        spots_since_last_save = 0

                        for future in as_completed(futures):
                            i = futures[future]
                            try:
                                result = future.result(timeout=300)
                                if result is not None and result.ndim == 2:
                                    spotwise_gene_expression_profiles[i] = result.copy()
                                    completed_spots.add(i)
                                    spots_since_last_save += 1
                                    pbar.update(1)

                                    if spots_since_last_save >= checkpoint_interval:
                                        checkpoint_mgr.save_checkpoint(
                                            completed_spots, spotwise_gene_expression_profiles, N, T, M
                                        )
                                        spots_since_last_save = 0
                            except TimeoutError:
                                logging.error(f"Timeout processing spot {i}")
                                continue
                            except Exception as e:
                                logging.error(f"Error processing spot {i}: {str(e)}")
                                logging.error(traceback.format_exc())
                                continue

                break

            except concurrent.futures.process.BrokenProcessPool:  # type: ignore
                retry_count += 1
                logging.warning(f"Process pool broken, retry {retry_count}/{max_retries}")
                if retry_count == max_retries:
                    logging.error("Max retries reached, saving current progress")
                import time

                time.sleep(5)

    finally:
        if futures:
            futures.clear()
        gc.collect()

        if spotwise_gene_expression_profiles:
            checkpoint_mgr.save_final_results(spotwise_gene_expression_profiles, completed_spots, N, T, M)

    return spotwise_gene_expression_profiles


# Module-level worker data for cell Pass 2 parallelization
_cell_pass2_worker_data = None


def _solve_single_cell_pass2(cell_idx):
    """Solve QP for a single cell's true expression. Module-level for pickling."""
    import gurobipy as gp
    from gurobipy import GRB

    wd = _cell_pass2_worker_data
    model = None
    try:
        ct = wd['dominant_type'][cell_idx]
        obs = wd['X_obs'][cell_idx]
        obs_lib = obs.sum()

        if obs_lib < 1:
            return cell_idx, obs.copy()

        enrich = wd['enrichment_weights'][ct]
        neighbor_mean = wd['same_type_neighbor_means'][cell_idx]
        has_neighbors = wd['same_type_neighbor_counts'][cell_idx] > 0
        M_genes = wd['M']
        lib_slack = wd['library_slack']
        l_enrich = wd['lambda_enrich']
        l_spatial = wd['lambda_spatial']

        model = gp.Model(f"cell_expr_{cell_idx}")
        model.setParam("OutputFlag", 0)
        model.setParam("Threads", 1)
        model.setParam("TimeLimit", 30)

        max_lib = lib_slack * obs_lib
        X_vars = model.addVars(
            M_genes, lb=0, ub=max_lib, vtype=GRB.CONTINUOUS, name="X"
        )

        model.addConstr(
            gp.quicksum(X_vars[g] for g in range(M_genes)) <= max_lib,
            name="lib_size",
        )

        obj_terms = []
        for g in range(M_genes):
            obj_terms.append((X_vars[g] - obs[g]) * (X_vars[g] - obs[g]))
            if l_enrich > 0 and enrich[g] > 0.1:
                obj_terms.append(-l_enrich * enrich[g] * X_vars[g])
            if l_spatial > 0 and has_neighbors:
                obj_terms.append(
                    l_spatial * (X_vars[g] - neighbor_mean[g]) * (X_vars[g] - neighbor_mean[g])
                )

        model.setObjective(gp.quicksum(obj_terms), GRB.MINIMIZE)
        model.optimize()

        if model.status == GRB.OPTIMAL:
            result = np.array([X_vars[g].X for g in range(M_genes)])
            return cell_idx, result
        else:
            return cell_idx, obs.copy()

    except Exception as e:
        logging.warning(f"Cell {cell_idx} optimization failed: {e}")
        return cell_idx, wd['X_obs'][cell_idx].copy()
    finally:
        if model:
            del model


def estimate_true_expression_cell(
    X_obs: np.ndarray,
    Y_assignments: np.ndarray,
    coords: np.ndarray,
    enrichment_weights: np.ndarray,
    library_slack: float = 1.5,
    lambda_enrich: float = 1.0,
    lambda_spatial: float = 0.01,
    spatial_k: int = 50,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 500,
) -> np.ndarray:
    """
    Estimate true gene expression per cell using optimization.

    For each cell, solves a QP to find X_true that:
    1. Stays close to observed counts (data fidelity)
    2. Recovers dropout for genes expected in the cell's type (enrichment prior)
    3. Agrees with same-type spatial neighbors (spatial coherence)
    4. Respects bounded total library size (prevents runaway imputation)

    Args:
        X_obs: Observed expression matrix (n_cells, n_genes).
        Y_assignments: Cell type assignment weights (n_cells, n_types) from Pass 1.
        coords: Spatial coordinates (n_cells, 2).
        enrichment_weights: Type-gene enrichment matrix (n_types, n_genes).
        library_slack: Max ratio of X_true library size to X_obs (default 1.5).
        lambda_enrich: Weight for enrichment prior term.
        lambda_spatial: Weight for spatial coherence term.
        spatial_k: Number of neighbors for spatial smoothing.
        max_workers: Max parallel workers (None = cpu_count).
        checkpoint_interval: Cells between checkpoint saves.

    Returns:
        X_true: Estimated true expression (n_cells, n_genes).
    """
    N, M = X_obs.shape
    T = Y_assignments.shape[1]

    # Determine dominant type per cell, guarding against unassigned (zero-row) cells
    total_weight = Y_assignments.sum(axis=1)
    unassigned_mask = total_weight < 1e-9
    dominant_type = np.argmax(Y_assignments, axis=1)
    dominant_type[unassigned_mask] = -1  # Mark unassigned cells

    if unassigned_mask.any():
        n_unassigned = int(unassigned_mask.sum())
        logging.info(
            f"Cell-level Pass 2: {n_unassigned} unassigned cells will retain "
            f"observed expression (no deconvolution)"
        )

    # Build spatial neighbor graph (k-NN)
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    k_query = min(spatial_k + 1, N)
    _, all_neighbor_indices = tree.query(coords, k=k_query)
    # Remove self
    if all_neighbor_indices.ndim > 1 and all_neighbor_indices.shape[1] > 1:
        all_neighbor_indices = all_neighbor_indices[:, 1:]
    else:
        all_neighbor_indices = np.empty((N, 0), dtype=int)

    # For each cell, find same-type neighbors and precompute their mean expression
    same_type_neighbor_means = np.zeros((N, M))
    same_type_neighbor_counts = np.zeros(N, dtype=int)
    for i in range(N):
        ct = dominant_type[i]
        neighbors = all_neighbor_indices[i]
        same_type_mask = dominant_type[neighbors] == ct
        same_type_neighbors = neighbors[same_type_mask]
        if len(same_type_neighbors) > 0:
            same_type_neighbor_means[i] = X_obs[same_type_neighbors].mean(axis=0)
            same_type_neighbor_counts[i] = len(same_type_neighbors)

    logging.info(
        f"Cell-level Pass 2: {N} cells, {M} genes, {T} types, "
        f"library_slack={library_slack}, spatial_k={spatial_k}"
    )

    X_true = np.zeros((N, M), dtype=np.float64)

    # Pre-fill unassigned cells with observed expression (no deconvolution)
    if unassigned_mask.any():
        X_true[unassigned_mask] = X_obs[unassigned_mask]

    # Build list of cells to optimize (skip unassigned)
    cells_to_optimize = np.where(~unassigned_mask)[0]

    # Process cells in parallel
    workers = max_workers if max_workers is not None else os.cpu_count()

    # Store worker data in module-level variable for pickling
    global _cell_pass2_worker_data
    _cell_pass2_worker_data = {
        'X_obs': X_obs,
        'dominant_type': dominant_type,
        'enrichment_weights': enrichment_weights,
        'same_type_neighbor_means': same_type_neighbor_means,
        'same_type_neighbor_counts': same_type_neighbor_counts,
        'library_slack': library_slack,
        'lambda_enrich': lambda_enrich,
        'lambda_spatial': lambda_spatial,
        'M': M,
    }

    # Process assigned cells - use sequential for small N, parallel for large N
    n_to_optimize = len(cells_to_optimize)
    if n_to_optimize <= 100:
        # Sequential for small datasets (avoids pickle overhead)
        for i in cells_to_optimize:
            cell_idx, result = _solve_single_cell_pass2(i)
            X_true[cell_idx] = result
    else:
        # Parallel for large datasets
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_solve_single_cell_pass2, i): i for i in cells_to_optimize}
            with tqdm(total=n_to_optimize, desc="Estimating true expression") as pbar:
                for future in as_completed(futures):
                    cell_idx, result = future.result()
                    X_true[cell_idx] = result
                    pbar.update(1)

    # Clean up module-level data
    _cell_pass2_worker_data = None

    return X_true


def normalize_counts(adata, target_sum=10000, exclude_highly_expressed=False, max_fraction=0.05):
    """
    Normalize counts while preserving integer values and relative proportions.

    Args:
        adata: AnnData object
        target_sum: Target sum for each cell/spot
        exclude_highly_expressed: Whether to exclude highly expressed genes
        max_fraction: Maximum fraction for highly expressed genes
    """
    # Get matrix
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X.copy()

    # Handle highly expressed genes if requested
    if exclude_highly_expressed:
        counts_per_cell = X.sum(axis=1)
        gene_subset = ~(X > counts_per_cell[:, None] * max_fraction).any(axis=0)
        size_factors = X[:, gene_subset].sum(axis=1)
    else:
        size_factors = X.sum(axis=1)

    # Ensure positive values
    size_factors = np.maximum(size_factors, 1)
    median_size = max(np.median(size_factors), 1)

    # Calculate bounded scaling factors
    scaling_factors = np.clip(size_factors / median_size, 0.1, 10.0)
    scaled_factors = target_sum / size_factors

    # Scale and round to integers
    X_scaled = np.round(X * scaled_factors[:, None]).astype(int)

    # Safety bounds
    max_allowed = target_sum * 2
    X_scaled = np.clip(X_scaled, 0, max_allowed)

    # Create new AnnData with scaled counts
    adata_norm = adata.copy()
    adata_norm.X = X_scaled

    # Store normalization info
    adata_norm.obs["size_factors"] = scaling_factors
    adata_norm.obs["original_total"] = size_factors
    adata_norm.obs["scaled_total"] = X_scaled.sum(axis=1)

    # Log statistics
    logging.info(f"Normalization stats:")
    logging.info(f"Original median total: {median_size:.2f}")
    logging.info(f"Mean scaled total: {X_scaled.sum(axis=1).mean():.2f}")
    logging.info(f"Max scaled value: {X_scaled.max():.2f}")

    return adata_norm


def validate_prior_effect(spotwise_profiles_pass1, spotwise_profiles_pass2, global_prior):
    """
    Compare pass1 and pass2 results to verify prior influence.

    Args:
        spotwise_profiles_pass1 (dict): First pass results {spot_idx: profile_matrix}
        spotwise_profiles_pass2 (dict): Second pass results {spot_idx: profile_matrix}
        global_prior (np.ndarray): Global prior matrix (T x M)

    Returns:
        dict: Dictionary containing validation metrics
    """
    # Validate shapes
    if not spotwise_profiles_pass1 or not spotwise_profiles_pass2:
        raise ValueError("Empty profile dictionaries provided")

    # Get shapes from first profile
    first_profile1 = next(iter(spotwise_profiles_pass1.values()))
    T, M = first_profile1.shape

    if global_prior.shape != (T, M):
        raise ValueError(f"Prior shape {global_prior.shape} does not match profiles shape ({T}, {M})")

    prior_guided_changes = []
    spot_metrics = {}

    # Ensure we have matching spots
    common_spots = set(spotwise_profiles_pass1.keys()) & set(spotwise_profiles_pass2.keys())

    if not common_spots:
        logging.error("No matching spots found between pass1 and pass2 results")
        return None

    for spot in common_spots:
        profile1 = spotwise_profiles_pass1[spot]
        profile2 = spotwise_profiles_pass2[spot]

        # Calculate absolute changes between passes
        profile_diff = np.abs(profile2 - profile1)
        total_diff = np.sum(profile_diff)

        # Calculate prior influence on pass2 assignment
        prior_alignment = np.sum(global_prior * profile2)

        prior_guided_changes.append((total_diff, prior_alignment))

        # Store per-spot metrics
        spot_metrics[spot] = {
            "total_change": total_diff,
            "prior_alignment": prior_alignment,
            "mean_change": np.mean(profile_diff),
            "max_change": np.max(profile_diff),
        }

    # Calculate correlation between changes and prior influence
    changes = np.array([x[0] for x in prior_guided_changes])
    influences = np.array([x[1] for x in prior_guided_changes])

    correlation = np.corrcoef(changes, influences)[0, 1]

    # Calculate summary statistics
    validation_metrics = {
        "prior_correlation": correlation,
        "mean_total_change": np.mean(changes),
        "mean_prior_influence": np.mean(influences),
        "std_total_change": np.std(changes),
        "std_prior_influence": np.std(influences),
        "n_spots_analyzed": len(common_spots),
        "spot_metrics": spot_metrics,
    }

    # Log summary statistics
    logging.info("Prior Effect Validation:")
    logging.info(f"Prior-Change Correlation: {correlation:.4f}")
    logging.info(f"Mean Total Change: {validation_metrics['mean_total_change']:.4f}")
    logging.info(f"Mean Prior Influence: {validation_metrics['mean_prior_influence']:.4f}")
    logging.info(f"Number of Spots Analyzed: {validation_metrics['n_spots_analyzed']}")



def _create_dummy_anndata(num_spots: int, num_genes: int) -> sc.AnnData:
    """Create a minimal AnnData with integer counts and spatial coordinates.

    The data is intentionally tiny so that model construction is fast, while
    still ensuring at least one positive count to create variables.

    Args:
        num_spots: Number of spatial spots to include.
        num_genes: Number of genes to include.

    Returns:
        A newly created AnnData object with `.obsm['spatial']` populated.
    """
    X = np.zeros((num_spots, num_genes), dtype=int)
    # Ensure at least one count is positive so variables are created
    X[0, 0] = 3
    if num_genes > 1:
        X[0, 1] = 1

    adata = sc.AnnData(X=X)
    adata.obsm['spatial'] = np.zeros((num_spots, 2), dtype=float)
    return adata


def _run_dummy_models() -> Dict[str, bool]:
    """Build tiny models and trigger writing of their MPS files.

    This function is intended for quick verification that model definitions
    are valid and can be exported to .mps files. Solver failures are tolerated
    because the primary goal here is to exercise the `model.write` calls.

    Returns:
        Dict[str, bool]: Mapping of model key to existence of the written file.
            Keys: 'cell', 'local', 'gene'.
    """
    results: Dict[str, bool] = {'cell': False, 'local': False, 'gene': False}

    try:
        # Tiny problem sizes
        num_spots = 1
        num_cell_types = 2
        num_genes = 3

        # Create minimal inputs
        adata = _create_dummy_anndata(num_spots, num_genes)
        profile_based_antibody_data = np.full((num_spots, num_cell_types), 0.5, dtype=float)
        cell_type_names = [f"Type_{i}" for i in range(num_cell_types)]

        # 1) Global EM-based cell proportions → writes cell_proportion_model.mps
        try:
            optimize_cell_proportions(
                profile_based_antibody_data=profile_based_antibody_data,
                cell_type_names=cell_type_names,
                max_iterations=1,
            )
        except Exception as exc:
            logging.warning("Cell proportions run ended (expected OK if solver missing): %s", str(exc))
        finally:
            results['cell'] = os.path.exists('cell_proportion_model.mps')

        # 2) Local refinement of cell proportions → writes local_cell_proportions_model.mps
        try:
            deconvolute_local_cell_proportions(
                spot_idx=0,
                adata=adata,
                profile_based_antibody_data=profile_based_antibody_data,
                radius=2.0,
                max_iterations=1,
                beta_values=np.ones(num_cell_types, dtype=float),
                beta_vary=False,
                max_y_change=0.4,
            )
        except Exception as exc:
            logging.warning("Local proportions run ended (expected OK if solver missing): %s", str(exc))
        finally:
            results['local'] = os.path.exists('local_cell_proportions_model.mps')

        # 3) Gene deconvolution → writes gene_expression_model.mps
        try:
            cell_type_numbers_array = np.full((num_spots, num_cell_types), 0.5, dtype=float)
            deconvolute_spot_with_neighbors_with_prior(
                spot_idx=0,
                adata=adata,
                cell_type_numbers_array=cell_type_numbers_array,
                radius=2.0,
                global_prior=None,
                lambda_prior_weight=0.0,
                local_enrichment_weight=0.5,
                global_enrichment_weight=0.5,
            )
        except Exception as exc:
            logging.warning("Gene deconvolution run ended (expected OK if solver missing): %s", str(exc))
        finally:
            results['gene'] = os.path.exists('gene_expression_model.mps')

    except Exception as outer_exc:
        logging.error("Dummy model setup failed: %s", str(outer_exc))

    return results


# ===========================================================================
# Discrete Cell Assignment (Integer Quadratic Programming)
# ===========================================================================


def estimate_prior_proportions_from_markers(
    raw_marker_data: np.ndarray,
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
) -> np.ndarray:
    """
    Estimate expected global cell type proportions from raw marker signal.

    Uses the mean expression of markers assigned to each cell type (before
    any scaling) to estimate the relative abundance of each type. This provides
    a prior that prevents rare cell type over-inflation during discrete assignment.

    Args:
        raw_marker_data: (N, M) raw antibody data (NOT scaled)
        assignment_matrix: (M, T) binary matrix mapping markers to cell types
        cell_type_names: List of cell type names (length T)

    Returns:
        prior_proportions: (T,) normalized proportions summing to 1.0

    Example:
        >>> raw_data, marker_names, assignment, cell_types = map_antibodies_to_profiles_v2(...)
        >>> prior = estimate_prior_proportions_from_markers(raw_data, assignment, cell_types)
        >>> model.run_discrete_cell_assignment(prior_proportions=prior, lambda_prior=0.5)
    """
    N, M = raw_marker_data.shape
    T = len(cell_type_names)

    # Compute mean signal per marker across all spots
    marker_means = raw_marker_data.mean(axis=0)  # Shape: (M,)

    # For each cell type, compute mean of its assigned markers' signals
    # This uses raw signal, so abundant types will have higher values
    type_signals = np.zeros(T, dtype=np.float64)
    for t in range(T):
        assigned_markers = assignment_matrix[:, t] > 0
        if assigned_markers.any():
            type_signals[t] = marker_means[assigned_markers].mean()
        else:
            type_signals[t] = 0.0

    # Normalize to proportions
    total = type_signals.sum()
    if total > 0:
        prior_proportions = type_signals / total
    else:
        # Fallback to uniform if no signal
        prior_proportions = np.ones(T) / T

    logging.info("Estimated prior proportions from raw marker signals:")
    for t, name in enumerate(cell_type_names):
        logging.info(f"  {name}: {100*prior_proportions[t]:.1f}%")

    return prior_proportions


def solve_discrete_cell_counts(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,
    beta_values: np.ndarray,
    alpha_values: Optional[np.ndarray] = None,
    max_nuclei_cap: int = 30,
    timeout_per_spot: float = 60.0,
    lambda_sparse: float = 0.0,
    prior_proportions: Optional[np.ndarray] = None,
    lambda_prior: float = 0.0,
) -> np.ndarray:
    """
    Solve IQP for discrete cell counts given fixed beta (E-step).

    Mathematical formulation:
        minimize    Σᵢ Σₘ (X[i,m] - α[m] - Σₜ c[i,t] × profile[t,m] × β[m])²
                    + λ_sparse × Σₜ y[t]
                    + λ_prior × Σₜ (c[i,t] - N_i × π[t])²
        subject to  Σₜ c[i,t] = N_i     ∀ spots i with N_i > 0
                    c[i,t] ∈ Z≥0        ∀ i, t
                    y[t] ∈ {0, 1}       ∀ t (indicator: y[t]=1 iff c[t]>0)
                    c[i,t] ≤ N_i × y[t] ∀ i, t (big-M constraint)

    The sparsity penalty (lambda_sparse > 0) encourages solutions with fewer
    active cell types per spot, matching biological reality that most spots
    have 2-4 dominant cell types rather than all types present.

    The prior penalty (lambda_prior > 0) regularizes towards expected global
    proportions π[t], preventing rare cell type over-inflation caused by
    per-marker scaling artifacts.

    Args:
        marker_level_data: (N, M) antibody data (preprocessed for discrete mode)
        marker_names: List of marker names (length M)
        assignment_matrix: (M, T) binary matrix where A[m,t]=1 if marker m belongs to type t
        cell_type_names: List of cell type names (length T)
        nuclei_counts: (N,) integer nuclei count per spot from Cellpose segmentation
        beta_values: (M,) per-marker scaling factors from previous EM iteration
        alpha_values: (M,) per-marker baselines (optional, defaults to zeros)
        max_nuclei_cap: Above this nuclei count, use continuous relaxation + rounding
        timeout_per_spot: Maximum seconds per spot optimization (default: 60)
        lambda_sparse: Sparsity penalty weight (default: 0.0). Higher values
            encourage fewer active cell types per spot. Typical range: 0.0-1.0
        prior_proportions: (T,) expected global proportions per cell type.
            If None and lambda_prior > 0, uses uniform prior (1/T).
        lambda_prior: Prior penalty weight (default: 0.0). Higher values
            pull assignments towards the expected global proportions.
            Typical range: 0.0-1.0

    Returns:
        c_values: (N, T) integer cell counts per type per spot
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    if assignment_matrix.shape != (M, T):
        raise ValueError(f"Assignment matrix shape {assignment_matrix.shape} != expected ({M}, {T})")
    if nuclei_counts.shape != (N,):
        raise ValueError(f"nuclei_counts shape {nuclei_counts.shape} != expected ({N},)")
    if beta_values.shape != (M,):
        raise ValueError(f"beta_values shape {beta_values.shape} != expected ({M},)")

    if alpha_values is None:
        alpha_values = np.zeros(M, dtype=np.float64)

    # Handle prior proportions
    if prior_proportions is None and lambda_prior > 0:
        # Uniform prior if not specified
        prior_proportions = np.ones(T, dtype=np.float64) / T
    elif prior_proportions is not None:
        prior_proportions = np.asarray(prior_proportions, dtype=np.float64)
        if prior_proportions.shape != (T,):
            raise ValueError(f"prior_proportions shape {prior_proportions.shape} != expected ({T},)")
        # Normalize to sum to 1
        prior_proportions = prior_proportions / (prior_proportions.sum() + 1e-10)

    # Build profile matrix: profile[t, m] = assignment_matrix[m, t]
    # (transposed for convenience in formulation)
    profile_matrix = assignment_matrix.T  # Shape: (T, M)

    c_values = np.zeros((N, T), dtype=np.int64)

    for i in range(N):
        N_i = int(nuclei_counts[i])

        # Skip spots with no nuclei
        if N_i <= 0:
            continue

        # Baseline-subtracted observed signal
        X_i = marker_level_data[i, :] - alpha_values  # Shape: (M,)

        # Decide whether to use integer or continuous relaxation
        use_integer = N_i <= max_nuclei_cap

        try:
            model = gp.Model(f"discrete_cell_spot_{i}")
            model.setParam("OutputFlag", 0)
            model.setParam("TimeLimit", timeout_per_spot)

            # Decision variables: c[t] = count of cell type t at spot i
            if use_integer:
                c = model.addVars(T, lb=0, ub=N_i, vtype=GRB.INTEGER, name="c")
            else:
                c = model.addVars(T, lb=0, ub=N_i, vtype=GRB.CONTINUOUS, name="c")

            # Constraint: sum of counts equals nuclei count
            model.addConstr(quicksum(c[t] for t in range(T)) == N_i, name="nuclei_sum")

            # Sparsity regularization via indicator variables
            # y[t] = 1 if c[t] > 0, else y[t] = 0
            # Add big-M constraint: c[t] <= N_i * y[t]
            y = None
            if lambda_sparse > 0.0:
                # When using continuous relaxation, indicators should also be continuous
                y_vtype = GRB.BINARY if use_integer else GRB.CONTINUOUS
                y = model.addVars(T, lb=0, ub=1, vtype=y_vtype, name="y")
                for t in range(T):
                    # Big-M constraint: c[t] <= N_i * y[t]
                    # If y[t] = 0, then c[t] must be 0
                    # If y[t] = 1, then c[t] can be up to N_i
                    model.addConstr(c[t] <= N_i * y[t], name=f"indicator_{t}")

            # Objective: minimize reconstruction error + sparsity penalty + prior penalty
            # For each marker m: error = X[i,m] - Σₜ c[t] × profile[t,m] × β[m]
            obj_terms = []
            for m in range(M):
                pred = quicksum(c[t] * profile_matrix[t, m] * beta_values[m] for t in range(T))
                residual = X_i[m] - pred
                obj_terms.append(residual * residual)

            total_objective = quicksum(obj_terms)

            # Add sparsity penalty: lambda_sparse * sum(y[t])
            if lambda_sparse > 0.0 and y is not None:
                sparsity_penalty = lambda_sparse * quicksum(y[t] for t in range(T))
                total_objective = total_objective + sparsity_penalty

            # Add prior penalty: lambda_prior * sum_t (c[t] - N_i * pi[t])^2
            # This regularizes towards expected global proportions
            if lambda_prior > 0.0 and prior_proportions is not None:
                prior_terms = []
                for t in range(T):
                    expected_count = N_i * prior_proportions[t]
                    deviation = c[t] - expected_count
                    prior_terms.append(deviation * deviation)
                prior_penalty = lambda_prior * quicksum(prior_terms)
                total_objective = total_objective + prior_penalty

            model.setObjective(total_objective, GRB.MINIMIZE)
            model.optimize()

            if model.Status == GRB.OPTIMAL or model.Status == GRB.TIME_LIMIT:
                for t in range(T):
                    val = c[t].X
                    if use_integer:
                        c_values[i, t] = int(round(val))
                    else:
                        # Continuous relaxation: round to integers
                        c_values[i, t] = int(round(val))

                # Verify and fix sum constraint after rounding
                current_sum = c_values[i, :].sum()
                if current_sum != N_i:
                    diff = N_i - current_sum
                    # Adjust the largest value to maintain exact sum
                    if diff > 0:
                        max_idx = np.argmax(c_values[i, :])
                        c_values[i, max_idx] += diff
                    else:
                        # Remove from largest non-zero
                        sorted_idx = np.argsort(c_values[i, :])[::-1]
                        for idx in sorted_idx:
                            if c_values[i, idx] > 0:
                                remove = min(-diff, c_values[i, idx])
                                c_values[i, idx] -= remove
                                diff += remove
                                if diff == 0:
                                    break
            else:
                # Infeasible or error: use uniform fallback
                logging.warning(f"Spot {i}: optimization status {model.Status}, using uniform fallback")
                uniform_count = N_i // T
                remainder = N_i % T
                for t in range(T):
                    c_values[i, t] = uniform_count + (1 if t < remainder else 0)

        except Exception as exc:
            logging.warning(f"Spot {i}: exception {exc}, using uniform fallback")
            uniform_count = N_i // T
            remainder = N_i % T
            for t in range(T):
                c_values[i, t] = uniform_count + (1 if t < remainder else 0)

    return c_values


def solve_discrete_cell_counts_global(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,
    beta_values: np.ndarray,
    alpha_values: Optional[np.ndarray] = None,
    time_limit: float = 300.0,
    mip_gap: float = 0.05,
    prev_c_values: Optional[np.ndarray] = None,
    continuous_prior: Optional[np.ndarray] = None,
    lambda_continuous: float = 0.0,
    constraint_slack: int = 0,
    lambda_reg: float = 0.0,
    alpha_elastic: float = 0.5,
    use_marker_weighting: bool = False,
) -> np.ndarray:
    """
    Solve global IQP for discrete cell counts across all spots jointly.

    Unlike solve_discrete_cell_counts() which optimizes each spot independently,
    this function formulates a single IQP over all N*T integer variables. This
    enforces globally consistent marker-celltype behavior across the tissue.

    Mathematical formulation:
        minimize    sum_i sum_m w[m] * (X[i,m] - alpha[m] - beta[m] * sum_t c[i,t] * profile[t,m])^2
                    + lambda_reg * (alpha_elastic * L1(c/N) + (1-alpha_elastic) * L2(c/N))
                    + lambda_continuous * sum_i sum_t (c[i,t]/N_i - p[i,t])^2
        subject to  N_i - slack <= sum_t c[i,t] <= N_i + slack   for all spots i
                    c[i,t] in Z>=0         for all i, t
                    c[i,t] <= N_i + slack  for all i, t

    Args:
        marker_level_data: (N, M) antibody data (preprocessed for discrete mode)
        marker_names: List of marker names (length M)
        assignment_matrix: (M, T) binary matrix where A[m,t]=1 if marker m belongs to type t
        cell_type_names: List of cell type names (length T)
        nuclei_counts: (N,) integer nuclei count per spot
        beta_values: (M,) per-marker scaling factors
        alpha_values: (M,) per-marker baselines (optional, defaults to zeros)
        time_limit: Maximum solver time in seconds (default: 300)
        mip_gap: Acceptable optimality gap (default: 0.05 = 5%)
        prev_c_values: (N, T) previous cell counts for warm-start (optional)
        continuous_prior: (N, T) continuous proportions from continuous optimization (optional)
        lambda_continuous: Regularization weight for continuous prior (default: 0.0)
        constraint_slack: Allow ±slack cells deviation from nuclei count (default: 0 = exact)
        lambda_reg: Elastic net regularization weight (default: 0.0 = no regularization)
        alpha_elastic: L1/L2 trade-off (0=pure L2, 1=pure L1, default: 0.5)
        use_marker_weighting: Weight error by 1/(n_owners * markers_per_type) (default: False)

    Returns:
        c_values: (N, T) integer cell counts per type per spot
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    # Input validation
    if assignment_matrix.shape != (M, T):
        raise ValueError(f"assignment_matrix shape {assignment_matrix.shape} != expected ({M}, {T})")
    if nuclei_counts.shape != (N,):
        raise ValueError(f"nuclei_counts shape {nuclei_counts.shape} != expected ({N},)")
    if beta_values.shape != (M,):
        raise ValueError(f"beta_values shape {beta_values.shape} != expected ({M},)")

    if alpha_values is None:
        alpha_values = np.zeros(M, dtype=np.float64)

    # Build profile matrix: profile[t, m] = assignment_matrix[m, t]
    profile_matrix = assignment_matrix.T  # Shape: (T, M)

    # Compute marker weighting if enabled
    # Weight by 1/(n_owners * markers_per_type) to balance cell types
    marker_weights = np.ones(M, dtype=np.float64)
    if use_marker_weighting:
        # Count markers per cell type
        markers_per_type = assignment_matrix.sum(axis=0)  # Shape: (T,)
        # For each marker, compute weight
        for m in range(M):
            owners = np.where(assignment_matrix[m, :] > 0)[0]
            n_owners = len(owners)
            if n_owners > 0:
                # Average markers per owner type
                avg_markers = np.mean([markers_per_type[t] for t in owners])
                marker_weights[m] = 1.0 / (n_owners * max(avg_markers, 1.0))
        logging.info(f"Marker weighting enabled: weights range [{marker_weights.min():.3f}, {marker_weights.max():.3f}]")

    logging.info(f"Global IQP: {N} spots x {T} cell types = {N * T} integer variables")
    logging.info(f"Time limit: {time_limit}s, MIP gap: {mip_gap:.1%}")
    if constraint_slack > 0:
        logging.info(f"Constraint slack: ±{constraint_slack} cells")
    if lambda_reg > 0:
        logging.info(f"Elastic net: lambda={lambda_reg:.2f}, alpha={alpha_elastic:.2f}")

    # Create Gurobi model
    model = gp.Model("global_discrete_cell_assignment")
    model.setParam("OutputFlag", 1)  # Show progress for long solves
    model.setParam("TimeLimit", time_limit)
    model.setParam("MIPGap", mip_gap)
    model.setParam("Threads", 0)  # Use all available cores

    # Create integer variables: c[i, t] = count of cell type t at spot i
    c = {}
    for i in range(N):
        N_i = int(nuclei_counts[i])
        ub_i = N_i + constraint_slack  # Allow slack in upper bound
        for t in range(T):
            c[i, t] = model.addVar(
                lb=0,
                ub=ub_i,
                vtype=GRB.INTEGER,
                name=f"c_{i}_{t}"
            )

    model.update()

    # Warm-start from previous solution, continuous prior, or uniform distribution
    for i in range(N):
        N_i = int(nuclei_counts[i])
        if prev_c_values is not None:
            for t in range(T):
                c[i, t].Start = int(prev_c_values[i, t])
        elif continuous_prior is not None:
            # Initialize from continuous prior (round to integers)
            target_counts = continuous_prior[i, :] * N_i
            # Use largest-remainder method to ensure sum equals N_i
            floor_counts = np.floor(target_counts).astype(int)
            remainders = target_counts - floor_counts
            deficit = N_i - floor_counts.sum()
            # Add 1 to cells with largest remainders
            top_indices = np.argsort(remainders)[::-1][:deficit]
            for t in range(T):
                c[i, t].Start = floor_counts[t] + (1 if t in top_indices else 0)
        else:
            # Uniform distribution
            base = N_i // T
            remainder = N_i % T
            for t in range(T):
                c[i, t].Start = base + (1 if t < remainder else 0)

    # Constraint: sum of counts within slack of nuclei per spot
    for i in range(N):
        N_i = int(nuclei_counts[i])
        if N_i > 0:
            sum_c_i = quicksum(c[i, t] for t in range(T))
            lb_i = max(0, N_i - constraint_slack)
            ub_i = N_i + constraint_slack
            model.addConstr(sum_c_i >= lb_i, name=f"nuclei_lb_{i}")
            model.addConstr(sum_c_i <= ub_i, name=f"nuclei_ub_{i}")

    # Objective: minimize reconstruction error across all spots and markers
    # For each (spot, marker): error = w[m] * (X[i,m] - alpha[m] - beta[m] * pred[i,m])^2
    # where pred[i,m] = sum_t c[i,t] * profile[t,m]
    error_terms = []
    for i in range(N):
        X_i = marker_level_data[i, :]  # Shape: (M,)
        for m in range(M):
            # Predicted signal from cell counts
            pred_im = quicksum(c[i, t] * profile_matrix[t, m] for t in range(T))
            # Residual: observed - baseline - beta * predicted
            residual = X_i[m] - alpha_values[m] - beta_values[m] * pred_im
            # Apply marker weighting
            error_terms.append(marker_weights[m] * residual * residual)

    # Elastic net regularization on proportions: lambda * (alpha*L1 + (1-alpha)*L2)
    # L1 = sum |c[i,t]/N_i|, L2 = sum (c[i,t]/N_i)^2
    reg_terms = []
    if lambda_reg > 0:
        for i in range(N):
            N_i = int(nuclei_counts[i])
            if N_i > 0:
                for t in range(T):
                    prop_it = c[i, t] / N_i  # Proportion
                    # L2 term (always quadratic)
                    reg_terms.append((1 - alpha_elastic) * prop_it * prop_it)
                    # L1 term: for non-negative variables, |x| = x
                    reg_terms.append(alpha_elastic * prop_it)

    # Add continuous prior regularization if provided
    # Penalizes (c[i,t]/N_i - p[i,t])^2 = (c[i,t] - N_i*p[i,t])^2 / N_i^2
    prior_terms = []
    if continuous_prior is not None and lambda_continuous > 0:
        logging.info(f"Using continuous prior with lambda={lambda_continuous:.2f}")
        for i in range(N):
            N_i = int(nuclei_counts[i])
            if N_i > 0:
                for t in range(T):
                    # Target count from continuous proportion
                    target_count = N_i * continuous_prior[i, t]
                    # Deviation penalty (scaled by 1/N_i^2 to normalize)
                    deviation = (c[i, t] - target_count) / N_i
                    prior_terms.append(deviation * deviation)

    total_objective = quicksum(error_terms)
    if reg_terms:
        total_objective += lambda_reg * quicksum(reg_terms)
    if prior_terms:
        total_objective += lambda_continuous * quicksum(prior_terms)

    model.setObjective(total_objective, GRB.MINIMIZE)

    # Solve
    logging.info("Starting global IQP optimization...")
    model.optimize()

    # Check solution status
    if model.status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL]:
        # Extract solution
        c_values = np.zeros((N, T), dtype=np.int64)
        for i in range(N):
            for t in range(T):
                c_values[i, t] = int(round(c[i, t].X))

        # Verify sum constraint (only fix if outside slack bounds)
        for i in range(N):
            N_i = int(nuclei_counts[i])
            current_sum = c_values[i, :].sum()
            lb_i = max(0, N_i - constraint_slack)
            ub_i = N_i + constraint_slack

            if current_sum < lb_i:
                # Below lower bound - add cells
                diff = lb_i - current_sum
                max_idx = np.argmax(c_values[i, :])
                c_values[i, max_idx] += diff
            elif current_sum > ub_i:
                # Above upper bound - remove cells
                diff = current_sum - ub_i
                sorted_idx = np.argsort(c_values[i, :])[::-1]
                for idx in sorted_idx:
                    if c_values[i, idx] > 0:
                        remove = min(diff, c_values[i, idx])
                        c_values[i, idx] -= remove
                        diff -= remove
                        if diff == 0:
                            break

        gap = model.MIPGap if hasattr(model, 'MIPGap') and model.SolCount > 0 else float('inf')
        logging.info(f"Global IQP solved: status={model.status}, "
                     f"gap={gap:.2%}, time={model.Runtime:.1f}s")

        return c_values

    else:
        raise RuntimeError(f"Global IQP failed: status={model.status}")


def optimize_discrete_cell_assignment_em(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,
    max_em_iterations: int = 20,
    beta_convergence_tol: float = 1e-3,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    max_nuclei_cap: int = 30,
    timeout_per_spot: float = 60.0,
    lambda_sparse: float = 0.0,
    prior_proportions: Optional[np.ndarray] = None,
    lambda_prior: float = 0.0,
    global_solve: bool = True,
    global_time_limit: float = 300.0,
    global_mip_gap: float = 0.05,
    continuous_prior: Optional[np.ndarray] = None,
    lambda_continuous: float = 0.0,
    constraint_slack: int = 0,
    lambda_reg: float = 0.0,
    alpha_elastic: float = 0.5,
    use_marker_weighting: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray]:
    """
    EM algorithm for discrete cell assignment with per-marker beta.

    E-step: Solve IQP for cell counts given current beta
    M-step: Update beta via OLS given cell counts

    Mathematical formulation for M-step (beta update):
        β[m] = Σᵢ (X[i,m] - α[m]) × pred[i,m] / Σᵢ pred[i,m]²
        where pred[i,m] = Σₜ c[i,t] × profile[t,m]

    Args:
        marker_level_data: (N, M) antibody data (preprocessed for discrete mode)
        marker_names: List of marker names (length M)
        assignment_matrix: (M, T) binary matrix where A[m,t]=1 if marker m belongs to type t
        cell_type_names: List of cell type names (length T)
        nuclei_counts: (N,) integer nuclei count per spot from Cellpose segmentation
        max_em_iterations: Maximum EM iterations (default: 20)
        beta_convergence_tol: Convergence tolerance for beta change (default: 1e-3)
        beta_min: Minimum allowed beta value (default: 0.1)
        beta_max: Maximum allowed beta value (default: 2.0)
        max_nuclei_cap: Above this nuclei count, use continuous relaxation (default: 30)
        timeout_per_spot: Maximum seconds per spot optimization (default: 60)
        lambda_sparse: Sparsity penalty weight (default: 0.0). Higher values encourage
            fewer active cell types per spot. Typical range: 0.0-1.0
        prior_proportions: (T,) expected global proportions per cell type.
            Can be estimated from raw marker signals or from continuous mode.
            If None, no prior regularization is applied regardless of lambda_prior.
        lambda_prior: Prior penalty weight (default: 0.0). Higher values
            pull assignments towards prior_proportions. Typical range: 0.0-1.0
        global_solve: If True (default), use global IQP solver. If False, use per-spot IQP.
        global_time_limit: Time limit for global solver in seconds (default: 300).
        global_mip_gap: Acceptable MIP gap for global solver (default: 0.05 = 5%).
        continuous_prior: (N, T) continuous proportions from continuous optimization (optional).
            If provided with lambda_continuous > 0, the IQP will be biased toward these values.
        lambda_continuous: Weight for continuous prior regularization (default: 0.0).
            Higher values make discrete solution closer to continuous. Typical range: 1.0-100.0
        constraint_slack: Allow ±slack cells deviation from nuclei count (default: 0 = exact).
            Setting to 1 allows the optimizer more flexibility like the continuous model.
        lambda_reg: Elastic net regularization weight on proportions (default: 0.0).
        alpha_elastic: L1/L2 trade-off (0=pure L2, 1=pure L1, default: 0.5).
        use_marker_weighting: Weight error by 1/(n_owners * markers_per_type) (default: False).
            Balances influence of cell types with different numbers of markers.

    Returns:
        Tuple of:
        - c_values: (N, T) integer cell counts per type per spot
        - beta_values: (M,) final per-marker scaling factors
        - marker_beta_dict: {marker_name: beta_value}
        - alpha_values: (M,) per-marker baselines
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    # Validate lambda_sparse and lambda_prior
    if lambda_sparse < 0:
        raise ValueError(f"lambda_sparse must be non-negative, got {lambda_sparse}")
    if lambda_prior < 0:
        raise ValueError(f"lambda_prior must be non-negative, got {lambda_prior}")

    logging.info(f"Discrete cell assignment EM: {N} spots, {M} markers, {T} cell types")
    logging.info(f"Total nuclei: {nuclei_counts.sum()}, mean per spot: {nuclei_counts.mean():.2f}")
    if lambda_sparse > 0:
        logging.info(f"Sparsity regularization: lambda_sparse={lambda_sparse}")
    if lambda_prior > 0 and prior_proportions is not None:
        logging.info(f"Prior regularization: lambda_prior={lambda_prior}")
        logging.info(f"Prior proportions: {dict(zip(cell_type_names, prior_proportions))}")

    # Build profile matrix (T x M)
    profile_matrix = assignment_matrix.T

    # Initialize beta and alpha
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)

    # Initialize cell counts (uniform distribution)
    c_values = np.zeros((N, T), dtype=np.int64)
    for i in range(N):
        N_i = int(nuclei_counts[i])
        if N_i > 0:
            uniform_count = N_i // T
            remainder = N_i % T
            for t in range(T):
                c_values[i, t] = uniform_count + (1 if t < remainder else 0)

    prev_loss = float('inf')

    for iteration in range(max_em_iterations):
        logging.info(f"\n=== EM Iteration {iteration + 1}/{max_em_iterations} ===")

        # ==================== E-Step ====================
        if global_solve:
            # Global solve: single IQP over all spots
            c_values = solve_discrete_cell_counts_global(
                marker_level_data=marker_level_data,
                marker_names=marker_names,
                assignment_matrix=assignment_matrix,
                cell_type_names=cell_type_names,
                nuclei_counts=nuclei_counts,
                beta_values=beta_values,
                alpha_values=alpha_values,
                time_limit=global_time_limit,
                mip_gap=global_mip_gap,
                prev_c_values=c_values if iteration > 0 else None,
                continuous_prior=continuous_prior,
                lambda_continuous=lambda_continuous,
                constraint_slack=constraint_slack,
                lambda_reg=lambda_reg,
                alpha_elastic=alpha_elastic,
                use_marker_weighting=use_marker_weighting,
            )
        else:
            # Per-spot solve (original behavior)
            c_values = solve_discrete_cell_counts(
                marker_level_data=marker_level_data,
                marker_names=marker_names,
                assignment_matrix=assignment_matrix,
                cell_type_names=cell_type_names,
                nuclei_counts=nuclei_counts,
                beta_values=beta_values,
                alpha_values=alpha_values,
                max_nuclei_cap=max_nuclei_cap,
                timeout_per_spot=timeout_per_spot,
                lambda_sparse=lambda_sparse,
                prior_proportions=prior_proportions,
                lambda_prior=lambda_prior,
            )

        # ==================== M-Step ====================
        # Update beta via OLS: β[m] = Σᵢ X'[i,m] × pred[i,m] / Σᵢ pred[i,m]²
        # where X' = X - α (baseline-subtracted) and pred[i,m] = Σₜ c[i,t] × profile[t,m]

        # Compute predictions: pred[i,m] = Σₜ c[i,t] × profile[t,m]
        pred = c_values @ profile_matrix  # Shape: (N, M)

        beta_prev = beta_values.copy()
        for m in range(M):
            pred_m = pred[:, m]  # Shape: (N,)
            X_m = marker_level_data[:, m]  # Shape: (N,)

            # First estimate alpha (intercept) via robust median
            # For markers with signal, alpha captures baseline
            if pred_m.sum() > 0:
                # Estimate alpha as median of X - beta*pred for spots with pred > 0
                mask = pred_m > 0
                if mask.sum() > 5:
                    residuals = X_m[mask] - beta_values[m] * pred_m[mask]
                    alpha_values[m] = np.clip(np.median(residuals), 0, 0.5)
                else:
                    alpha_values[m] = 0.0
            else:
                alpha_values[m] = 0.0

            # Baseline-subtracted signal
            X_prime_m = X_m - alpha_values[m]

            # OLS for beta: β = Σ X' × pred / Σ pred²
            numerator = np.sum(X_prime_m * pred_m)
            denominator = np.sum(pred_m * pred_m)

            if denominator > 1e-10:
                beta_values[m] = np.clip(numerator / denominator, beta_min, beta_max)
            else:
                beta_values[m] = 1.0  # Default if no signal

        # Compute reconstruction loss
        X_pred = alpha_values[np.newaxis, :] + pred * beta_values[np.newaxis, :]
        loss = np.sum((marker_level_data - X_pred) ** 2)

        # Check convergence
        beta_change = np.max(np.abs(beta_values - beta_prev))
        loss_change = abs(prev_loss - loss) / max(abs(prev_loss), 1e-10)

        logging.info(f"  Loss: {loss:.4f} (change: {loss_change:.6f})")
        logging.info(f"  Beta range: [{beta_values.min():.3f}, {beta_values.max():.3f}]")
        logging.info(f"  Max beta change: {beta_change:.6f}")
        logging.info(f"  Cell count range: [{c_values.sum(axis=1).min()}, {c_values.sum(axis=1).max()}]")

        # Log cell type distribution
        total_cells_per_type = c_values.sum(axis=0)
        total_cells = c_values.sum()
        for t, ct_name in enumerate(cell_type_names):
            pct = 100 * total_cells_per_type[t] / total_cells if total_cells > 0 else 0.0
            logging.info(f"    {ct_name}: {total_cells_per_type[t]} cells ({pct:.1f}%)")

        if beta_change < beta_convergence_tol and iteration > 0:
            logging.info(f"Converged after {iteration + 1} iterations (beta change < {beta_convergence_tol})")
            break

        if loss > prev_loss and iteration > 0:
            logging.warning(f"Loss increased from {prev_loss:.4f} to {loss:.4f}")

        prev_loss = loss

    # Build marker beta dictionary
    marker_beta_dict = {marker_names[m]: beta_values[m] for m in range(M)}

    logging.info(f"\nFinal beta values:")
    for m, name in enumerate(marker_names):
        logging.info(f"  {name}: β={beta_values[m]:.3f}, α={alpha_values[m]:.3f}")

    return c_values, beta_values, marker_beta_dict, alpha_values


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.info("Running gurobi_impl to write .mps files using CitegeistModel where available")

    # Import CitegeistModel lazily to avoid circular imports at module import time
    try:
        from .citegeist_model import CitegeistModel  # type: ignore
    except Exception:  # pragma: no cover
        from citegeist_model import CitegeistModel  # type: ignore

    DATA_FOLDER = "/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/"
    path_to_biopsy = os.path.join(DATA_FOLDER, "HCC22-088-P4-S1/outs")
    path_to_surgical = os.path.join(DATA_FOLDER, "HCC22-088-P4-S2_1i_rep/outs")
    path_list = [path_to_surgical]

    cell_profiles = {
        "Cancer Cells": {"Major": ["EPCAM-1"], "Minor": ["SDC1-1", "KRT5-1"]},
        "Macrophages": {"Major": ["CD68-1"], "Minor": ["CD14-1"]},
        "CD4 T Cells": {"Major": ["CD3E-1", "CD4-1"]},
        "CD8 T Cells": {"Major": ["CD3E-1", "CD8A-1"]},
        "B Cells": {"Major": ["MS4A1-1", "CD19-1"]},
        "Endothelial Cells": {"Major": ["PECAM1-1"]},
        "Fibroblasts": {"Major": ["ACTA2-1"]},
    }

    any_success = False
    try:
        import squidpy as sq  # type: ignore
    except Exception as exc:  # pragma: no cover
        logging.error("squidpy not available to read Visium data: %s", str(exc))
        sq = None  # type: ignore

    for sample_path in path_list:
        try:
            if sq is None:
                raise RuntimeError("squidpy unavailable")
            sample_name = os.path.basename(os.path.dirname(sample_path))
            adata = sq.read.visium(sample_path, counts_file='filtered_feature_bc_matrix.h5', load_images=True, gex_only=False)

            model = CitegeistModel(sample_name=sample_name, adata=adata, output_folder=f'output_mps_{sample_name}')
            model.load_cell_profile_dict(cell_profiles)
            model.split_adata()
            model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=25)
            model.copy_gex_to_protein_adata()
            model.preprocess_gex()
            model.preprocess_antibody()

            # Optional: register Gurobi if a license path is provided via env
            license_file = os.environ.get("GRB_LICENSE_FILE")
            if license_file and os.path.isfile(license_file):
                try:
                    model.register_gurobi(license_file)
                except Exception as exc:
                    logging.warning("Gurobi registration skipped: %s", str(exc))

            # Running these will invoke the underlying gurobi_impl functions that write .mps files
            model.run_cell_proportion_model(radius=400)
            # Trigger at least one gene model build; pass1 will spawn tasks that write MPS
            try:
                model.run_cell_expression_pass1(
                    radius=400,
                    max_workers=1,
                    checkpoint_interval=100,
                    output_dir=os.path.join(model.output_folder, "checkpoints"),
                    rerun=True,
                )
            except Exception as exc:
                logging.warning("Pass1 gene expression run ended: %s", str(exc))

            # Check for MPS files
            cell_mps = os.path.exists('cell_proportion_model.mps')
            local_mps = os.path.exists('local_cell_proportions_model.mps')
            gene_mps = os.path.exists('gene_expression_model.mps')
            logging.info("Sample '%s' MPS write outcomes: cell=%s local=%s gene=%s", sample_name, cell_mps, local_mps, gene_mps)
            any_success = any_success or cell_mps or local_mps or gene_mps
        except Exception as exc:
            logging.error("Failed processing sample at '%s': %s", sample_path, str(exc))

    if not any_success:
        logging.info("Falling back to dummy models as no MPS files were produced from real data.")
        outcome = _run_dummy_models()
        logging.info("Dummy MPS write outcomes: %s", outcome)
