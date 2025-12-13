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
except Exception:  # pragma: no cover - fallback for __main__ execution
    from utils import get_neighbors_with_fixed_radius  # type: ignore
    from checkpoints import CheckpointManager  # type: ignore


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

    Returns:
        Tuple[np.ndarray, np.ndarray]: Y_values (N x T), beta_values (T,)
    """
    import gurobipy as gp
    from gurobipy import GRB

    N, T = profile_based_antibody_data.shape

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

        model.setObjective(total_error + regularization_term + laplacian_term, GRB.MINIMIZE)

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


def assess_redundancy(
    profile_based_antibody_data: np.ndarray,
    Y_values: np.ndarray,
    cell_type_names: List[str],
    usage_threshold: float = 0.02,
    r2_threshold: float = 0.7,  # LOWERED from 0.9 for more aggressive filtering
) -> Dict[str, Any]:
    """
    Detect redundant/overspecified cell types via linear regression of marker patterns.

    A cell type is flagged as redundant if:
    1. It's substantially used (mean proportion ≥ usage_threshold, default 2%)
    2. Its marker pattern across spots is almost completely explained by a linear
       combination of other cell types' patterns (R² ≥ r2_threshold, default 0.7)

    This catches cannibalization: e.g., "fake Neutrophils" with markers from T-cells
    and Myeloid will have R² ≈ 1.0 because its pattern = 0.5×T-cell + 0.5×Myeloid.

    NOTE: r2_threshold lowered from 0.9 to 0.7 to more aggressively filter spurious
    profiles in automatic profile discovery scenarios.

    Args:
        profile_based_antibody_data (np.ndarray): Cell-type aggregated scores (N_spots × T_celltypes)
                                                   Each column is averaged marker signal for that type
        Y_values (np.ndarray): Inferred proportions (N_spots × T_celltypes)
        cell_type_names (List[str]): Cell type names
        usage_threshold (float): Minimum mean proportion to consider (default 0.02 = 2%)
        r2_threshold (float): R² threshold for redundancy (default 0.7 = 70% explained)
    
    Returns:
        Dict with:
        - redundant_fraction: Fraction of used types that are redundant
        - redundant_celltypes: List of dicts with {name, mean_usage, r2_explained_by_others, coeffs}
    
    Example:
        Three real types + one fake built from first two:
        >>> real1, real2, real3 = np.random.rand(3, 100)
        >>> fake = 0.5 * real1 + 0.5 * real2
        >>> scores = np.column_stack([real1, real2, real3, fake])
        >>> Y = np.full((100, 4), 0.25)  # All used equally
        >>> result = assess_redundancy(scores, Y, ["R1", "R2", "R3", "Fake"])
        >>> assert result["redundant_celltypes"][0]["name"] == "Fake"
        >>> assert result["redundant_celltypes"][0]["r2_explained_by_others"] > 0.95
    """
    N, T = profile_based_antibody_data.shape
    redundant = []
    
    for t, celltype in enumerate(cell_type_names):
        # Skip Unknown (has no markers by design)
        if celltype == "Unknown":
            continue
        
        # Check if this type is actually used in the proportions
        mean_usage = np.mean(Y_values[:, t])
        if mean_usage < usage_threshold:
            # Low usage types aren't concerning for redundancy
            continue
        
        # Get this cell type's averaged marker pattern across spots
        m_t = profile_based_antibody_data[:, t]  # (N,)
        
        # Get all OTHER cell types' patterns
        other_idx = [k for k in range(T) if k != t]
        M_others = profile_based_antibody_data[:, other_idx]  # (N, T-1)
        
        # Fit linear regression (no intercept): m_t ≈ M_others × β
        # This checks if this type's pattern is just a linear combo of others
        try:
            beta, residuals, rank, s = np.linalg.lstsq(M_others, m_t, rcond=None)
            m_hat = M_others @ beta
        except np.linalg.LinAlgError:
            logging.warning(f"Could not fit redundancy model for '{celltype}'")
            continue
        
        # Compute R²: fraction of variance in m_t explained by other types
        ss_res = np.sum((m_t - m_hat) ** 2)
        ss_tot = np.sum((m_t - m_t.mean()) ** 2)
        r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-10 else 0.0
        
        # Debug logging
        logging.info(f"REDUNDANCY CHECK: {celltype} usage={mean_usage:.3f} R²={r2:.4f}")
        
        # Flag if this type's pattern is almost fully explained by others
        if r2 >= r2_threshold:
            # Find which other types contribute most (for debugging)
            top_contributors_idx = np.argsort(np.abs(beta))[-3:][::-1]
            top_contributors = [
                (cell_type_names[other_idx[i]], float(beta[i])) 
                for i in top_contributors_idx
            ]
            
            redundant.append({
                "name": celltype,
                "mean_usage": float(mean_usage),
                "r2_explained_by_others": float(r2),
                "top_contributors": top_contributors,  # Most contributing types
            })
    
    # Compute fraction of USED types that are redundant
    used_types = [
        t for t in range(T) 
        if cell_type_names[t] != "Unknown" and np.mean(Y_values[:, t]) >= usage_threshold
    ]
    redundant_fraction = len(redundant) / len(used_types) if used_types else 0.0
    
    return {
        "redundant_fraction": float(redundant_fraction),
        "redundant_celltypes": redundant,
    }


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

    # Check 3: Redundancy validation (if antibody data provided)
    if profile_based_antibody_data is not None:
        result = assess_redundancy(
            profile_based_antibody_data,
            Y_values,
            cell_type_names,
            usage_threshold=0.02,   # Must be used ≥2% to flag
            r2_threshold=0.9,       # R² ≥ 0.9 means 90% explained by others
        )
        
        if result["redundant_fraction"] > redundancy_threshold:
            redundant_list = [
                f"'{ct['name']}' (usage={ct['mean_usage']*100:.1f}%, R²={ct['r2_explained_by_others']:.3f})"
                for ct in result["redundant_celltypes"]
            ]
            error_msg = (
                f"{len(result['redundant_celltypes'])} cell type(s) are redundant "
                f"(marker patterns explained by other types): {', '.join(redundant_list)}. "
                f"These types are likely overspecified or built from combinations of existing types. "
                f"Consider removing them from the cell profile dictionary."
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
    
    # Log redundancy results if available
    if profile_based_antibody_data is not None:
        result = assess_redundancy(profile_based_antibody_data, Y_values, cell_type_names)
        if result["redundant_celltypes"]:
            logging.info(f"  - Redundant types: {len(result['redundant_celltypes'])}/{len([t for t in cell_type_names if t != 'Unknown'])} "
                        f"(fraction: {result['redundant_fraction']:.1%})")
        else:
            logging.info(f"  - No redundant types detected")
    else:
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
            model.setObjective(total_error + reg_term, GRB.MINIMIZE)

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
) -> Optional[np.ndarray]:
    """
    Deconvolute a spot with its neighbors, using both enrichment weights and optional prior.
    """
    model = None
    # Local import to avoid hard dependency at module import time
    import gurobipy as gp  # type: ignore
    from gurobipy import GRB  # type: ignore
    try:
        neighborhood_indices = get_neighbors_with_fixed_radius(spot_idx, adata, radius=int(radius), include_center=True)
        if not neighborhood_indices:
            logging.error(f"No valid neighbors found for spot {spot_idx}.")
            return None

        neighborhood_indices = np.array(
            [int(idx) for idx in neighborhood_indices if isinstance(idx, (int, np.integer))], dtype=int
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
        inverse_frequency_weights = 1.0 / (celltype_frequencies + 1e-10)
        normalized_weights = inverse_frequency_weights / np.max(inverse_frequency_weights)

        # Modified enrichment calculation
        def compute_expression_aware_enrichment(expression_data, cell_type_props, gene_idx):
            """
            Compute expression-aware enrichment scores.

            Args:
                expression_data (np.ndarray): Expression matrix
                cell_type_props (np.ndarray): Cell type proportions
                gene_idx (int): Gene index

            Returns:
                np.ndarray: Enrichment scores for each cell type
            """
            gene_expr = expression_data[:, gene_idx]
            expr_threshold = np.percentile(gene_expr[gene_expr > 0], 50) if np.any(gene_expr > 0) else 0
            high_expr_spots = gene_expr >= expr_threshold

            if not np.any(high_expr_spots):
                return np.ones(cell_type_props.shape[1]) / cell_type_props.shape[1]

            # Normalize cell type proportions by their global frequency
            normalized_props = cell_type_props / (celltype_frequencies + 1e-10)

            high_expr_props = np.mean(normalized_props[high_expr_spots], axis=0)
            background_props = np.mean(normalized_props, axis=0)

            epsilon = 1e-10
            enrichment = high_expr_props / (background_props + epsilon)

            # Apply smoothing to avoid extreme values
            smoothed_enrichment = 0.8 * enrichment + 0.2 * np.ones_like(enrichment)
            return smoothed_enrichment / (np.sum(smoothed_enrichment) + epsilon)

        # Compute expression-aware enrichment for each gene
        gene_specific_enrichment = np.zeros((M, T))

        for k in range(M):
            local_enrich = compute_expression_aware_enrichment(
                neighborhood_expression_data, neighborhood_cell_type_numbers, k
            )
            global_enrich = compute_expression_aware_enrichment(
                deconvolution_expression_data, cell_type_numbers_array, k
            )
            gene_specific_enrichment[k] = (
                local_enrichment_weight * local_enrich + global_enrichment_weight * global_enrich
            )

        # Build Gurobi model
        model = gp.Model(f"discrete_gene_expression_spot_{spot_idx}")
        model.setParam("OutputFlag", 0)
        model.setParam("Threads", 1)
        model.setParam("NodefileStart", 0.5)
        model.setParam("MIPGap", 0.01)
        model.setParam("TimeLimit", 600)
        model.setParam("NodeLimit", 1000000)

        # Variables for count assignment
        X = {}
        center_counts = deconvolution_expression_data[spot_idx, :]

        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    X[j, k] = model.addVar(vtype=GRB.INTEGER, lb=0, ub=total_counts, name=f"X_{j}_{k}")
                # Count conservation constraint
                model.addConstr(gp.quicksum(X[j, k] for j in range(T)) == total_counts, name=f"count_conservation_{k}")

        # Validate prior if asked
        if global_prior is not None:
            if lambda_prior_weight > 0:
                if global_prior is None:
                    raise ValueError("lambda_prior_weight > 0 but no global_prior provided")
            if not isinstance(global_prior, np.ndarray):
                raise ValueError("global_prior must be a numpy array")
            if global_prior.shape != (T, M):
                raise ValueError(f"Prior matrix shape {global_prior.shape} does not match expected shape ({T}, {M})")

        # Modify objective terms to include frequency normalization
        obj_terms = []
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    # Get normalized weights
                    enrichment_weight = gene_specific_enrichment[k, j]
                    cell_type_weight = neighborhood_cell_type_numbers[len(neighborhood_indices) // 2, j]

                    # Apply frequency normalization
                    normalized_weight = cell_type_weight * normalized_weights[j]

                    # Add slight randomness to break ties using seeded RNG for reproducibility
                    rng = np.random.default_rng(42)  # Fixed seed for reproducibility
                    randomness = 0.9 + 0.2 * rng.random()
                    base_term = enrichment_weight * normalized_weight * randomness * X[j, k]
                    obj_terms.append(base_term)

                    # Prior terms remain unchanged
                    if global_prior is not None and lambda_prior_weight > 0:
                        try:
                            prior_value = float(global_prior[j, k])
                            prior_penalty = lambda_prior_weight * (1 - prior_value) * X[j, k]
                            obj_terms.append(-prior_penalty)
                        except Exception as e:
                            logging.warning(f"Error accessing prior at [{j}, {k}]: {str(e)}")
                            continue

        # Maximize the sum of all terms
        model.setObjective(gp.quicksum(obj_terms), GRB.MAXIMIZE)

        model.write('gene_expression_model.mps')

        model.optimize()

        if model.status == GRB.OPTIMAL:
            logging.info(f"Solution found for spot {spot_idx}")
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
) -> Dict[str, Any]:
    """
    Optimize gene expression with enrichment weights and prior guidance.

    Args:
        sample_name (str): Name of the sample
        deconvolution_expression_data (np.ndarray): Gene expression data (N_spots x M_genes)
        cell_type_numbers_array (np.ndarray): Cell type proportions (N_spots x T_celltypes)
        filtered_adata (sc.AnnData): Filtered AnnData object containing gene expression data
        radius (float): Radius for neighbor detection
        alpha (float): Weight for spatial regularization
        lambda_reg_gex (float): Weight for gene expression regularization
        global_enrichment_weight (float): Weight for global expression enrichment (0-1)
        local_enrichment_weight (float): Weight for local expression enrichment (0-1)
        global_prior (np.ndarray, optional): Global prior matrix for guidance
        lambda_prior_weight (float): Weight for prior guidance
        max_workers (int, optional): Maximum number of parallel workers
        checkpoint_interval (int): Number of spots between checkpoints
        output_dir (str): Directory for checkpoints
        rerun (bool): Whether to rerun if results exist

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


# =============================================================================
# MIQP-BASED PROFILE SELECTION FOR AUTO PROFILE DISCOVERY
# =============================================================================


# =============================================================================
# HIERARCHY INFRASTRUCTURE FOR HIERARCHICAL MIQP
# =============================================================================


def compute_hierarchy_levels(
    candidate_profiles: List[set],
) -> Dict[int, List[int]]:
    """
    Organize profiles by cardinality level.

    Level = |profile| - 1, so:
    - 1-marker profiles → level 0
    - 2-marker profiles → level 1
    - 3-marker profiles → level 2
    - etc.

    Args:
        candidate_profiles: List of candidate profiles (each a set of marker indices).

    Returns:
        Dictionary mapping level (int) to list of profile indices at that level.

    Example:
        >>> profiles = [{0}, {1}, {0, 1}, {0, 1, 2}]
        >>> compute_hierarchy_levels(profiles)
        {0: [0, 1], 1: [2], 2: [3]}
    """
    levels: Dict[int, List[int]] = {}

    for p_idx, profile in enumerate(candidate_profiles):
        cardinality = len(profile)
        if cardinality == 0:
            continue  # Skip empty profiles
        level = cardinality - 1
        if level not in levels:
            levels[level] = []
        levels[level].append(p_idx)

    return levels


def find_parent_child_relationships(
    candidate_profiles: List[set],
) -> Tuple[Dict[int, List[int]], Dict[int, List[int]]]:
    """
    Map each profile to its immediate parents (subsets with one fewer marker).

    A profile P is a parent of profile C if:
    1. P ⊂ C (P is a strict subset of C)
    2. |C| - |P| = 1 (C has exactly one more marker than P)

    Args:
        candidate_profiles: List of candidate profiles (each a set of marker indices).

    Returns:
        Tuple of:
        - parents_of: Dict mapping profile index to list of parent profile indices.
        - children_of: Dict mapping profile index to list of child profile indices.

    Example:
        >>> profiles = [{0}, {1}, {0, 1}, {0, 1, 2}]
        >>> parents, children = find_parent_child_relationships(profiles)
        >>> parents[2]  # {0,1} has parents {0} and {1}
        [0, 1]
        >>> children[0]  # {0} has child {0,1}
        [2]
    """
    n_profiles = len(candidate_profiles)
    parents_of: Dict[int, List[int]] = {p: [] for p in range(n_profiles)}
    children_of: Dict[int, List[int]] = {p: [] for p in range(n_profiles)}

    # Group profiles by cardinality for efficient search
    by_cardinality: Dict[int, List[int]] = {}
    for p_idx, profile in enumerate(candidate_profiles):
        card = len(profile)
        if card not in by_cardinality:
            by_cardinality[card] = []
        by_cardinality[card].append(p_idx)

    # For each profile, find parents (profiles with one fewer marker that are subsets)
    for p_idx, profile in enumerate(candidate_profiles):
        card = len(profile)
        if card <= 1:
            continue  # 1-marker profiles have no parents

        # Look for potential parents in the previous cardinality level
        potential_parents = by_cardinality.get(card - 1, [])
        for parent_idx in potential_parents:
            parent_profile = candidate_profiles[parent_idx]
            # Check if parent is a subset
            if parent_profile < profile:  # Strict subset
                parents_of[p_idx].append(parent_idx)
                children_of[parent_idx].append(p_idx)

    return parents_of, children_of


def compute_same_level_overlap_matrix(
    candidate_profiles: List[set],
    levels: Dict[int, List[int]],
) -> np.ndarray:
    """
    Compute overlap coefficient for profiles at the SAME level only.

    overlap[p, q] = |p ∩ q| / min(|p|, |q|) if p and q are at the same level, else 0.

    The overlap coefficient is 1.0 if one profile is a subset of another (but since
    same-level profiles have the same cardinality, this means they are identical),
    and 0.0 if they share no markers.

    Args:
        candidate_profiles: List of candidate profiles (each a set of marker indices).
        levels: Dictionary from compute_hierarchy_levels().

    Returns:
        Symmetric overlap matrix (n_candidates, n_candidates) with values in [0, 1].
        Off-diagonal entries are non-zero only for same-level profile pairs.

    Example:
        >>> profiles = [{0, 1}, {1, 2}, {3, 4}]  # All level 1 (2-marker)
        >>> levels = {1: [0, 1, 2]}
        >>> overlap = compute_same_level_overlap_matrix(profiles, levels)
        >>> overlap[0, 1]  # {0,1} and {1,2} share marker 1
        0.5
        >>> overlap[0, 2]  # {0,1} and {3,4} share nothing
        0.0
    """
    n_candidates = len(candidate_profiles)
    overlap = np.zeros((n_candidates, n_candidates))

    # Only compute overlap for profiles at the same level
    for level, profile_indices in levels.items():
        for i, p_idx in enumerate(profile_indices):
            profile_p = candidate_profiles[p_idx]
            len_p = len(profile_p)
            if len_p == 0:
                continue

            for q_idx in profile_indices[i+1:]:  # Only upper triangle
                profile_q = candidate_profiles[q_idx]
                len_q = len(profile_q)
                if len_q == 0:
                    continue

                intersection = len(profile_p & profile_q)
                min_size = min(len_p, len_q)

                if min_size > 0:
                    overlap_coef = intersection / min_size
                    overlap[p_idx, q_idx] = overlap_coef
                    overlap[q_idx, p_idx] = overlap_coef  # Symmetric

    return overlap


def compute_markers_at_level(
    candidate_profiles: List[set],
    levels: Dict[int, List[int]],
) -> Dict[int, Dict[int, List[int]]]:
    """
    For each level, compute which profiles contain each marker.

    This is used for same-level marker exclusion constraints.

    Args:
        candidate_profiles: List of candidate profiles.
        levels: Dictionary from compute_hierarchy_levels().

    Returns:
        Nested dict: {level: {marker_idx: [profile_indices_at_that_level]}}

    Example:
        >>> profiles = [{0}, {1}, {0, 1}, {0, 2}]
        >>> levels = {0: [0, 1], 1: [2, 3]}
        >>> markers_at_level = compute_markers_at_level(profiles, levels)
        >>> markers_at_level[1][0]  # At level 1, marker 0 is in profiles 2 and 3
        [2, 3]
    """
    markers_at_level: Dict[int, Dict[int, List[int]]] = {}

    for level, profile_indices in levels.items():
        markers_at_level[level] = {}
        for p_idx in profile_indices:
            profile = candidate_profiles[p_idx]
            for marker in profile:
                if marker not in markers_at_level[level]:
                    markers_at_level[level][marker] = []
                markers_at_level[level][marker].append(p_idx)

    return markers_at_level


def compute_profile_spatial_score(
    marker_scores: np.ndarray,
    profile_markers: set,
    aggregation: str = "mean",
) -> float:
    """
    Aggregate marker spatial scores into a profile-level score.

    Args:
        marker_scores: Array of spatial scores for each marker (n_markers,).
        profile_markers: Set of marker indices in this profile.
        aggregation: How to combine scores - "mean", "min", "max", or "geometric".

    Returns:
        Profile-level spatial score in [0, 1].
    """
    if not profile_markers:
        return 0.0

    scores = np.array([marker_scores[m] for m in profile_markers])

    if aggregation == "mean":
        return float(np.mean(scores))
    elif aggregation == "min":
        return float(np.min(scores))
    elif aggregation == "max":
        return float(np.max(scores))
    elif aggregation == "geometric":
        # Geometric mean, with small epsilon to avoid log(0)
        return float(np.exp(np.mean(np.log(scores + 1e-6))))
    else:
        raise ValueError(f"Unknown aggregation method: {aggregation}")


def optimize_profiles_miqp(
    Z: np.ndarray,
    candidate_profiles: List[set],
    marker_spatial_scores: np.ndarray,
    lambda_spatial: float = 0.1,
    lambda_complexity: float = 0.01,
    max_profiles: int = 20,
    min_profiles: int = 2,
    allow_overlap: bool = False,
    time_limit: float = 300.0,
    mip_gap: float = 0.01,
    seed: int = 1234,
    verbose: bool = True,
) -> Tuple[List[set], np.ndarray, Dict[str, Any]]:
    """
    Joint MIQP optimization for profile selection.

    This function jointly optimizes profile selection and reconstruction quality,
    using spatial scores as soft penalties rather than hard filters.

    Args:
        Z: Standardized expression matrix (n_spots, n_markers).
        candidate_profiles: List of candidate profiles (each a set of marker indices).
        marker_spatial_scores: Composite spatial score for each marker (n_markers,).
            Should be in [0, 1], higher = stronger spatial signal.
        lambda_spatial: Weight for spatial penalty (default 0.1 = balanced).
        lambda_complexity: Weight for profile count penalty (default 0.01).
        max_profiles: Maximum number of profiles to select.
        min_profiles: Minimum number of profiles to select.
        allow_overlap: If False, each marker can appear in at most one selected profile.
        time_limit: Solver time limit in seconds.
        mip_gap: Acceptable MIP optimality gap.
        seed: Random seed for reproducibility.
        verbose: Whether to log progress.

    Returns:
        selected_profiles: List of selected profile marker sets.
        proportions: Proportion matrix (n_spots, n_selected).
        metadata: Dict with optimization info (objective, gap, runtime, etc.).
    """
    logger = logging.getLogger(__name__)

    n_spots, n_markers = Z.shape
    n_candidates = len(candidate_profiles)

    if n_candidates == 0:
        if verbose:
            logger.warning("No candidate profiles provided to MIQP optimizer")
        return [], np.zeros((n_spots, 0)), {"status": "no_candidates"}

    # Compute profile-level spatial scores (mean aggregation)
    profile_spatial_scores = np.array([
        compute_profile_spatial_score(marker_spatial_scores, profile, "mean")
        for profile in candidate_profiles
    ])

    # Build profile-marker matrix A (P x M): A[p,m] = 1 if marker m is in profile p
    A = np.zeros((n_candidates, n_markers))
    for p, profile in enumerate(candidate_profiles):
        for m in profile:
            A[p, m] = 1.0

    # Initialize beta (marker weights) - uniform for now
    beta = np.ones(n_markers)

    if verbose:
        logger.info(
            f"MIQP optimization: {n_candidates} candidates, {n_spots} spots, "
            f"{n_markers} markers"
        )
        logger.info(
            f"Parameters: lambda_spatial={lambda_spatial}, "
            f"lambda_complexity={lambda_complexity}, allow_overlap={allow_overlap}"
        )

    try:
        # Create Gurobi model
        model = Model("ProfileMIQP")
        model.setParam("OutputFlag", 0)  # Suppress Gurobi output
        model.setParam("TimeLimit", time_limit)
        model.setParam("MIPGap", mip_gap)
        model.setParam("Seed", seed)

        # === VARIABLES ===

        # z_p: binary selection of profile p
        z = model.addVars(n_candidates, vtype=GRB.BINARY, name="z")

        # Y_sp: proportion of profile p in spot s (continuous 0-1)
        Y = model.addVars(
            n_spots, n_candidates,
            lb=0.0, ub=1.0,
            vtype=GRB.CONTINUOUS,
            name="Y"
        )

        model.update()

        # === CONSTRAINTS ===

        # Profile activation: Y_sp <= z_p (can only use profile if selected)
        for s in range(n_spots):
            for p in range(n_candidates):
                model.addConstr(Y[s, p] <= z[p], name=f"act_{s}_{p}")

        # Coverage constraints per spot (sum of proportions should be reasonable)
        for s in range(n_spots):
            model.addConstr(
                quicksum(Y[s, p] for p in range(n_candidates)) >= 0.5,
                name=f"cov_lb_{s}"
            )
            model.addConstr(
                quicksum(Y[s, p] for p in range(n_candidates)) <= 1.5,
                name=f"cov_ub_{s}"
            )

        # Profile count bounds
        model.addConstr(
            quicksum(z[p] for p in range(n_candidates)) >= min_profiles,
            name="min_prof"
        )
        model.addConstr(
            quicksum(z[p] for p in range(n_candidates)) <= max_profiles,
            name="max_prof"
        )

        # Marker overlap constraint (non-hierarchical mode)
        if not allow_overlap:
            for m in range(n_markers):
                profiles_with_marker = [p for p in range(n_candidates) if A[p, m] > 0]
                if len(profiles_with_marker) > 1:
                    model.addConstr(
                        quicksum(z[p] for p in profiles_with_marker) <= 1,
                        name=f"no_overlap_{m}"
                    )

        # === OBJECTIVE ===
        # We use a linearized approximation since full quadratic over all s,m,p
        # combinations is computationally expensive.
        #
        # Objective = reconstruction_proxy + spatial_penalty + complexity_penalty
        #
        # Reconstruction proxy: For each profile, measure how well it explains
        # the markers it contains. A profile is valuable if its markers are
        # co-expressed across spots.

        # Compute profile "quality" scores based on co-expression
        profile_quality = []
        for p, profile in enumerate(candidate_profiles):
            markers = list(profile)
            if len(markers) == 0:
                profile_quality.append(0.0)
                continue
            # Profile quality = variance explained by mean expression
            marker_data = Z[:, markers]
            profile_expr = marker_data.mean(axis=1)
            predicted = np.outer(profile_expr, np.ones(len(markers)))
            before_error = np.sum(marker_data ** 2)
            after_error = np.sum((marker_data - predicted) ** 2)
            if before_error > 1e-9:
                quality = (before_error - after_error) / before_error
            else:
                quality = 0.0
            profile_quality.append(max(0.0, quality))
        profile_quality = np.array(profile_quality)

        # Objective: maximize quality, minimize spatial penalty and complexity
        # Reformulated as: minimize -quality + spatial_penalty + complexity_penalty
        objective = quicksum(
            -profile_quality[p] * z[p]  # Maximize quality (negated for min)
            + lambda_spatial * (1 - profile_spatial_scores[p]) * z[p]  # Spatial penalty
            + lambda_complexity * z[p]  # Complexity penalty
            for p in range(n_candidates)
        )

        model.setObjective(objective, GRB.MINIMIZE)

        # === SOLVE ===
        model.optimize()

        if model.status not in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL]:
            if verbose:
                logger.warning(f"MIQP solver status: {model.status}")
            return [], np.zeros((n_spots, 0)), {
                "status": "failed",
                "solver_status": model.status
            }

        # === EXTRACT RESULTS ===
        selected_indices = [p for p in range(n_candidates) if z[p].X > 0.5]
        selected_profiles = [candidate_profiles[p] for p in selected_indices]

        # Extract proportions
        proportions = np.zeros((n_spots, len(selected_indices)))
        for i, p in enumerate(selected_indices):
            for s in range(n_spots):
                proportions[s, i] = Y[s, p].X

        # Normalize proportions to sum to 1 per spot
        row_sums = proportions.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-9)  # Avoid division by zero
        proportions = proportions / row_sums

        metadata = {
            "status": "optimal" if model.status == GRB.OPTIMAL else "suboptimal",
            "objective": model.objVal,
            "mip_gap": model.MIPGap if hasattr(model, 'MIPGap') else None,
            "runtime": model.Runtime,
            "n_selected": len(selected_profiles),
            "n_candidates": n_candidates,
            "profile_quality_scores": {
                str(candidate_profiles[p]): profile_quality[p]
                for p in selected_indices
            },
            "profile_spatial_scores": {
                str(candidate_profiles[p]): profile_spatial_scores[p]
                for p in selected_indices
            },
        }

        if verbose:
            logger.info(
                f"MIQP completed: selected {len(selected_profiles)}/{n_candidates} profiles"
            )
            logger.info(f"Objective value: {model.objVal:.4f}")

        return selected_profiles, proportions, metadata

    except gp.GurobiError as e:
        logger.error(f"Gurobi error in MIQP: {e}")
        return [], np.zeros((n_spots, 0)), {
            "status": "gurobi_error",
            "error": str(e)
        }
    except Exception as e:
        logger.error(f"Unexpected error in MIQP: {e}")
        return [], np.zeros((n_spots, 0)), {
            "status": "error",
            "error": str(e)
        }


def compute_spot_profile_fit(
    Z: np.ndarray,
    candidate_profiles: List[set],
    aggregation: str = "mean",
) -> np.ndarray:
    """
    Compute how well each profile fits each spot based on marker expression.

    fit[s, p] = aggregate over markers m in profile p of Z_norm[s, m]

    Higher fit means the profile's markers are well-expressed at this spot.
    Lower fit means the profile shouldn't have high proportion at this spot.

    Args:
        Z: Standardized expression matrix (n_spots, n_markers).
        candidate_profiles: List of candidate profiles (each a set of marker indices).
        aggregation: How to combine marker signals - "mean", "min", or "geometric".
            - "mean": Average of marker expressions (lenient)
            - "min": Minimum of marker expressions (strict - all markers must be present)
            - "geometric": Geometric mean (balanced)

    Returns:
        fit_matrix: (n_spots, n_candidates) with values in [0, 1].
    """
    n_spots = Z.shape[0]
    n_candidates = len(candidate_profiles)
    fit = np.zeros((n_spots, n_candidates))

    # Normalize Z to [0, 1] for fit scoring (per marker)
    Z_min = Z.min(axis=0)
    Z_max = Z.max(axis=0)
    Z_range = Z_max - Z_min + 1e-9
    Z_norm = (Z - Z_min) / Z_range

    for p, profile in enumerate(candidate_profiles):
        markers = list(profile)
        if len(markers) == 0:
            continue

        # Get marker signals at each spot
        marker_signals = Z_norm[:, markers]

        if aggregation == "mean":
            fit[:, p] = marker_signals.mean(axis=1)
        elif aggregation == "min":
            fit[:, p] = marker_signals.min(axis=1)
        elif aggregation == "geometric":
            # Geometric mean with small epsilon to avoid log(0)
            fit[:, p] = np.exp(np.log(marker_signals + 1e-9).mean(axis=1))
        else:
            raise ValueError(f"Unknown aggregation method: {aggregation}")

    return fit


def optimize_profiles_miqp_hierarchical(
    Z: np.ndarray,
    candidate_profiles: List[set],
    marker_spatial_scores: np.ndarray,
    marker_names: Optional[List[str]] = None,
    snr_weights: Optional[Dict[str, float]] = None,
    lambda_spatial: float = 0.1,
    lambda_complexity: float = 0.05,  # INCREASED from 0.01 for stronger penalty
    lambda_overlap: float = 0.7,       # INCREASED from 0.5 for stronger penalty
    lambda_orphan: float = 0.2,
    lambda_sparsity: float = 0.5,      # INCREASED from 0.3 for stronger penalty
    lambda_snr: float = 0.2,           # NEW: SNR-based penalty weight
    min_quality_threshold: float = 0.1,  # NEW: quality floor for profiles
    max_profiles: int = 20,
    min_profiles: int = 2,
    enforce_hierarchy: bool = False,
    sparsity_aggregation: str = "mean",
    time_limit: float = 300.0,
    mip_gap: float = 0.01,
    seed: int = 1234,
    verbose: bool = True,
) -> Tuple[List[set], np.ndarray, Dict[str, Any]]:
    """
    Hierarchical MIQP optimization for profile selection with spot-level sparsity.

    This function extends the basic MIQP with:
    1. Cardinality-based hierarchy: Profiles organized by marker count (level = |p| - 1)
    2. Same-level competition penalty: Quadratic penalty for profiles at same level sharing markers
    3. Orphan penalty: Soft penalty for child profiles without selected parents
    4. Spot-level sparsity: Penalizes Y[s,p] when profile p's markers aren't expressed at spot s
    5. SNR-based penalty: Penalizes profiles with low signal-to-noise ratio markers

    Args:
        Z: Standardized expression matrix (n_spots, n_markers).
        candidate_profiles: List of candidate profiles (each a set of marker indices).
        marker_spatial_scores: Composite spatial score for each marker (n_markers,).
        marker_names: Optional list of marker names (required if snr_weights provided).
        snr_weights: Optional dict mapping marker names to SNR values from SMM.
            Used to penalize profiles containing low-SNR markers.
        lambda_spatial: Weight for spatial penalty (default 0.1).
        lambda_complexity: Weight for profile count penalty (default 0.05).
            Increased from 0.01 for stronger profile count control.
        lambda_overlap: Weight for same-level overlap penalty (default 0.7).
            Increased from 0.5 for stronger competition penalty.
        lambda_orphan: Weight for orphan penalty (default 0.2).
            Penalizes child profiles when no parent is selected.
        lambda_sparsity: Weight for spot-level sparsity penalty (default 0.5).
            Increased from 0.3 for stronger sparsity enforcement.
        lambda_snr: Weight for SNR-based penalty (default 0.2).
            Penalizes profiles with low average SNR.
        min_quality_threshold: Quality floor for profiles (default 0.1).
            Profiles with quality below this threshold are excluded.
        max_profiles: Maximum number of profiles to select.
        min_profiles: Minimum number of profiles to select.
        enforce_hierarchy: If True, use hard constraint (child requires parent).
            If False (default), use soft penalty for orphans.
        sparsity_aggregation: How to compute spot-profile fit - "mean", "min", "geometric".
        time_limit: Solver time limit in seconds.
        mip_gap: Acceptable MIP optimality gap.
        seed: Random seed for reproducibility.
        verbose: Whether to log progress.

    Returns:
        selected_profiles: List of selected profile marker sets.
        proportions: Proportion matrix (n_spots, n_selected).
        metadata: Dict with optimization info including hierarchy details.
    """
    logger = logging.getLogger(__name__)

    n_spots, n_markers = Z.shape
    n_candidates = len(candidate_profiles)

    if n_candidates == 0:
        if verbose:
            logger.warning("No candidate profiles provided to hierarchical MIQP")
        return [], np.zeros((n_spots, 0)), {"status": "no_candidates"}

    # === HIERARCHY COMPUTATION ===
    levels = compute_hierarchy_levels(candidate_profiles)
    parents_of, children_of = find_parent_child_relationships(candidate_profiles)
    overlap_matrix = compute_same_level_overlap_matrix(candidate_profiles, levels)
    markers_at_level = compute_markers_at_level(candidate_profiles, levels)

    if verbose:
        level_counts = {lv: len(indices) for lv, indices in levels.items()}
        logger.info(f"Hierarchy levels: {level_counts}")

    # === PROFILE QUALITY AND SPATIAL SCORES ===
    profile_spatial_scores = np.array([
        compute_profile_spatial_score(marker_spatial_scores, profile, "mean")
        for profile in candidate_profiles
    ])

    # Compute profile quality (variance explained)
    profile_quality = []
    for profile in candidate_profiles:
        markers = list(profile)
        if len(markers) == 0:
            profile_quality.append(0.0)
            continue
        marker_data = Z[:, markers]
        profile_expr = marker_data.mean(axis=1)
        predicted = np.outer(profile_expr, np.ones(len(markers)))
        before_error = np.sum(marker_data ** 2)
        after_error = np.sum((marker_data - predicted) ** 2)
        if before_error > 1e-9:
            quality = (before_error - after_error) / before_error
        else:
            quality = 0.0
        profile_quality.append(max(0.0, quality))
    profile_quality = np.array(profile_quality)

    # === SNR-BASED PROFILE SCORES ===
    # Compute per-profile SNR scores (average of constituent marker SNRs)
    if snr_weights is not None and marker_names is not None:
        profile_snr_scores = []
        for profile in candidate_profiles:
            marker_snrs = [
                snr_weights.get(marker_names[m], 1.0)
                for m in profile
            ]
            profile_snr_scores.append(np.mean(marker_snrs) if marker_snrs else 0.0)
        profile_snr_scores = np.array(profile_snr_scores)
        # Normalize to [0, 1]
        if profile_snr_scores.max() > 0:
            profile_snr_scores = profile_snr_scores / profile_snr_scores.max()
        if verbose:
            logger.info(f"SNR weights applied: mean SNR score = {np.mean(profile_snr_scores):.3f}")
    else:
        # No SNR weights - all profiles get score of 1.0 (no penalty)
        profile_snr_scores = np.ones(n_candidates)

    # === QUALITY FLOOR: Identify low-quality profiles ===
    low_quality_profiles = []
    if min_quality_threshold > 0:
        for p in range(n_candidates):
            if profile_quality[p] < min_quality_threshold:
                low_quality_profiles.append(p)
        if verbose and low_quality_profiles:
            logger.info(
                f"Quality floor: {len(low_quality_profiles)} profiles below "
                f"threshold {min_quality_threshold}"
            )

    # === SPOT-PROFILE FIT FOR SPARSITY ===
    fit_matrix = compute_spot_profile_fit(Z, candidate_profiles, sparsity_aggregation)

    if verbose:
        logger.info(
            f"Hierarchical MIQP: {n_candidates} candidates, {n_spots} spots, "
            f"{n_markers} markers"
        )
        logger.info(
            f"Parameters: lambda_overlap={lambda_overlap}, lambda_orphan={lambda_orphan}, "
            f"lambda_sparsity={lambda_sparsity}, lambda_snr={lambda_snr}, "
            f"enforce_hierarchy={enforce_hierarchy}"
        )

    try:
        # === BUILD MODEL ===
        model = Model("HierarchicalProfileMIQP")
        model.setParam("OutputFlag", 0)
        model.setParam("TimeLimit", time_limit)
        model.setParam("MIPGap", mip_gap)
        model.setParam("Seed", seed)

        # === VARIABLES ===
        # z[p]: binary selection of profile p
        z = model.addVars(n_candidates, vtype=GRB.BINARY, name="z")

        # Y[s,p]: proportion of profile p in spot s
        Y = model.addVars(
            n_spots, n_candidates,
            lb=0.0, ub=1.0,
            vtype=GRB.CONTINUOUS,
            name="Y"
        )

        # Auxiliary variables for linearizing quadratic overlap penalty
        # w[p,q] = z[p] * z[q] for same-level overlapping pairs
        overlap_pairs = []
        for level, profile_indices in levels.items():
            for i, p_idx in enumerate(profile_indices):
                for q_idx in profile_indices[i+1:]:
                    if overlap_matrix[p_idx, q_idx] > 0:
                        overlap_pairs.append((p_idx, q_idx))

        w = {}
        for p_idx, q_idx in overlap_pairs:
            w[p_idx, q_idx] = model.addVar(vtype=GRB.BINARY, name=f"w_{p_idx}_{q_idx}")

        model.update()

        # === CONSTRAINTS ===

        # Profile activation: Y[s,p] <= z[p]
        for s in range(n_spots):
            for p in range(n_candidates):
                model.addConstr(Y[s, p] <= z[p], name=f"act_{s}_{p}")

        # Coverage constraints per spot
        for s in range(n_spots):
            model.addConstr(
                quicksum(Y[s, p] for p in range(n_candidates)) >= 0.5,
                name=f"cov_lb_{s}"
            )
            model.addConstr(
                quicksum(Y[s, p] for p in range(n_candidates)) <= 1.5,
                name=f"cov_ub_{s}"
            )

        # Profile count bounds
        model.addConstr(
            quicksum(z[p] for p in range(n_candidates)) >= min_profiles,
            name="min_prof"
        )
        model.addConstr(
            quicksum(z[p] for p in range(n_candidates)) <= max_profiles,
            name="max_prof"
        )

        # Quality floor constraint: exclude very low quality profiles
        for p in low_quality_profiles:
            model.addConstr(z[p] == 0, name=f"quality_floor_{p}")

        # Same-level marker exclusion (hard constraint for competing profiles)
        for level, marker_dict in markers_at_level.items():
            for marker, profile_list in marker_dict.items():
                if len(profile_list) > 1:
                    # At most one profile at this level can use this marker
                    model.addConstr(
                        quicksum(z[p] for p in profile_list) <= 1,
                        name=f"level_{level}_marker_{marker}_excl"
                    )

        # Linearize quadratic overlap: w[p,q] = z[p] * z[q]
        for p_idx, q_idx in overlap_pairs:
            model.addConstr(w[p_idx, q_idx] <= z[p_idx], name=f"w_ub1_{p_idx}_{q_idx}")
            model.addConstr(w[p_idx, q_idx] <= z[q_idx], name=f"w_ub2_{p_idx}_{q_idx}")
            model.addConstr(
                w[p_idx, q_idx] >= z[p_idx] + z[q_idx] - 1,
                name=f"w_lb_{p_idx}_{q_idx}"
            )

        # Hierarchy constraint (hard or soft)
        if enforce_hierarchy:
            # Hard constraint: child can only be selected if at least one parent is selected
            for p_idx, parent_list in parents_of.items():
                if parent_list:  # Has parents
                    model.addConstr(
                        z[p_idx] <= quicksum(z[parent] for parent in parent_list),
                        name=f"hier_hard_{p_idx}"
                    )

        # === OBJECTIVE ===
        # Base terms: quality, spatial, complexity
        obj_terms = []

        for p in range(n_candidates):
            # Maximize quality (negate for minimization)
            obj_terms.append(-profile_quality[p] * z[p])

            # Spatial penalty
            obj_terms.append(lambda_spatial * (1 - profile_spatial_scores[p]) * z[p])

            # Complexity penalty
            obj_terms.append(lambda_complexity * z[p])

            # SNR penalty: penalize low-SNR profiles
            # profile_snr_scores is normalized to [0, 1], so penalty = 1 - score
            snr_penalty = lambda_snr * (1.0 - profile_snr_scores[p])
            obj_terms.append(snr_penalty * z[p])

        # Same-level overlap penalty (quadratic, linearized)
        for p_idx, q_idx in overlap_pairs:
            overlap_coef = overlap_matrix[p_idx, q_idx]
            obj_terms.append(lambda_overlap * overlap_coef * w[p_idx, q_idx])

        # Orphan penalty (soft hierarchy)
        if not enforce_hierarchy:
            for p_idx, parent_list in parents_of.items():
                if parent_list:  # Has parents
                    # If profile is selected but no parent is selected, add penalty
                    # orphan_penalty = z[p] * (1 - max(z[parent] for parent in parents))
                    # We approximate this by: penalty if z[p]=1 and sum(z[parents])=0
                    # Using auxiliary variable to track if any parent is selected
                    n_parents = len(parent_list)
                    # Simple approximation: penalize proportionally to how few parents selected
                    # penalty = z[p] - (1/n_parents) * sum(z[parent] * z[p])
                    # For simplicity, just add penalty scaled by whether parents exist
                    obj_terms.append(
                        lambda_orphan * z[p_idx] * (1.0 / (1.0 + n_parents))
                    )

        # Spot-level sparsity penalty: penalize Y[s,p] when fit[s,p] is low
        # This is LINEAR in Y[s,p] since fit[s,p] are constants
        for s in range(n_spots):
            for p in range(n_candidates):
                misfit = 1.0 - fit_matrix[s, p]
                if misfit > 0.01:  # Only add meaningful penalties
                    obj_terms.append(lambda_sparsity * misfit * Y[s, p])

        model.setObjective(quicksum(obj_terms), GRB.MINIMIZE)

        # === SOLVE ===
        model.optimize()

        if model.status not in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL]:
            if verbose:
                logger.warning(f"Hierarchical MIQP solver status: {model.status}")
            return [], np.zeros((n_spots, 0)), {
                "status": "failed",
                "solver_status": model.status
            }

        # === EXTRACT RESULTS ===
        selected_indices = [p for p in range(n_candidates) if z[p].X > 0.5]
        selected_profiles = [candidate_profiles[p] for p in selected_indices]

        # Extract proportions
        proportions = np.zeros((n_spots, len(selected_indices)))
        for i, p in enumerate(selected_indices):
            for s in range(n_spots):
                proportions[s, i] = Y[s, p].X

        # Normalize proportions
        row_sums = proportions.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-9)
        proportions = proportions / row_sums

        # Build metadata
        selected_levels = {
            p: [lv for lv, indices in levels.items() if p in indices][0]
            for p in selected_indices
        }

        metadata = {
            "status": "optimal" if model.status == GRB.OPTIMAL else "suboptimal",
            "objective": model.objVal,
            "mip_gap": model.MIPGap if hasattr(model, 'MIPGap') else None,
            "runtime": model.Runtime,
            "n_selected": len(selected_profiles),
            "n_candidates": n_candidates,
            "hierarchy_levels": levels,
            "parents_of": parents_of,
            "children_of": children_of,
            "selected_levels": selected_levels,
            "profile_quality_scores": {
                str(candidate_profiles[p]): profile_quality[p]
                for p in selected_indices
            },
            "profile_spatial_scores": {
                str(candidate_profiles[p]): profile_spatial_scores[p]
                for p in selected_indices
            },
            "overlap_pairs_count": len(overlap_pairs),
        }

        if verbose:
            logger.info(
                f"Hierarchical MIQP completed: selected {len(selected_profiles)}/{n_candidates}"
            )
            logger.info(f"Selected by level: {dict(pd.Series(list(selected_levels.values())).value_counts())}")
            logger.info(f"Objective value: {model.objVal:.4f}")

        return selected_profiles, proportions, metadata

    except gp.GurobiError as e:
        logger.error(f"Gurobi error in hierarchical MIQP: {e}")
        return [], np.zeros((n_spots, 0)), {
            "status": "gurobi_error",
            "error": str(e)
        }
    except Exception as e:
        logger.error(f"Unexpected error in hierarchical MIQP: {e}")
        return [], np.zeros((n_spots, 0)), {
            "status": "error",
            "error": str(e)
        }


def optimize_profiles_miqp_with_fallback(
    Z: np.ndarray,
    candidate_profiles: List[set],
    marker_spatial_scores: np.ndarray,
    fallback_func: callable,
    fallback_kwargs: Dict[str, Any],
    **miqp_kwargs,
) -> Tuple[List[set], np.ndarray, Dict[str, Any]]:
    """
    Run MIQP optimization with fallback to reconstruction-based selection.

    Args:
        Z: Standardized expression matrix.
        candidate_profiles: List of candidate profiles.
        marker_spatial_scores: Spatial scores for each marker.
        fallback_func: Function to call if MIQP fails.
        fallback_kwargs: Kwargs to pass to fallback function.
        **miqp_kwargs: Additional kwargs for optimize_profiles_miqp.

    Returns:
        Same as optimize_profiles_miqp.
    """
    logger = logging.getLogger(__name__)

    try:
        profiles, proportions, metadata = optimize_profiles_miqp(
            Z, candidate_profiles, marker_spatial_scores, **miqp_kwargs
        )

        if metadata.get("status") in ("failed", "gurobi_error", "error"):
            raise RuntimeError(f"MIQP failed with status: {metadata.get('status')}")

        if len(profiles) == 0:
            raise RuntimeError("MIQP returned no profiles")

        return profiles, proportions, metadata

    except Exception as e:
        logger.warning(f"MIQP failed ({e}), falling back to reconstruction-based selection")
        # Call fallback function
        profiles = fallback_func(**fallback_kwargs)
        # Create dummy proportions and metadata
        n_spots = Z.shape[0]
        proportions = np.ones((n_spots, len(profiles))) / max(len(profiles), 1)
        metadata = {
            "status": "fallback",
            "fallback_reason": str(e),
            "n_selected": len(profiles),
        }
        return profiles, proportions, metadata


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


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.info("Running gurobi_impl to write .mps files using CitegeistModel where available")

    # Import CitegeistModel lazily to avoid circular imports at module import time
    try:
        from .citegeist_model import CitegeistModel  # type: ignore
    except Exception:  # pragma: no cover
        from citegeist_model import CitegeistModel  # type: ignore

    DATA_FOLDER = "/bgfs/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/"
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
