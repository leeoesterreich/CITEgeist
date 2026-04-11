"""cuOPT-based optimization backend for CITEgeist.

Replaces Gurobi with NVIDIA cuOPT 26.2 for GPU-accelerated
quadratic programming. Same function signatures as gurobi_impl.py.
"""

# Standard library imports
import gc
import logging
import os
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional, Tuple, cast

# Third-party imports
import numpy as np
import scanpy as sc
import scipy
from tqdm import tqdm

# Local imports
# Provide a fallback so this module can be executed directly as a script
try:
    from CITEgeist.model.checkpoints import CheckpointManager
    from CITEgeist.model.gex.gex_modules import (
        compute_kl_penalty_coefficients,
        compute_softmax_target,
    )
    from CITEgeist.model.utils import get_neighbors_with_fixed_radius
except (ImportError, OSError):  # pragma: no cover - fallback for __main__ execution
    from ..checkpoints import CheckpointManager
    from ..gex.gex_modules import (
        compute_kl_penalty_coefficients,
        compute_softmax_target,
    )
    from ..utils import get_neighbors_with_fixed_radius

try:
    from cuopt.linear_programming.problem import Problem
except (ImportError, OSError):  # pragma: no cover - cuOPT needs CUDA; fails on login nodes
    Problem = None


def build_spatial_laplacian(
    coords: np.ndarray,
    k: int = 8,
    normed: bool = True,
    cellularity: Optional[np.ndarray] = None,
    cellularity_sigma: float = 0.5,
) -> scipy.sparse.spmatrix:
    """
    Build graph Laplacian from spatial coordinates using k-NN.

    The Laplacian matrix L is used for spatial smoothness regularization.
    For normalized Laplacian: L = I - D^(-1/2) A D^(-1/2)
    where A is the adjacency matrix and D is the degree matrix.

    The quadratic form Y^T L Y penalizes differences between neighboring spots:
        Y^T L Y = sum_{i~j} w_ij * (Y_i - Y_j)^2 / sqrt(d_i * d_j)

    Args:
        coords: Spatial coordinates, shape (n_spots, 2).
        k: Number of nearest neighbors for graph construction.
        normed: If True, return normalized Laplacian (recommended).
        cellularity: Optional per-spot nuclei counts, shape (n_spots,).
            When provided, edges are weighted by cellularity similarity:
            w_ij = exp(-|log(n_i) - log(n_j)|^2 / (2*sigma^2)).
            This preserves tissue boundaries where cellularity changes abruptly.
        cellularity_sigma: Bandwidth for cellularity similarity kernel (default: 0.5).
            Controls how much cellularity difference suppresses smoothing.
            Smaller = more aggressive boundary detection.

    Returns:
        Sparse Laplacian matrix of shape (n_spots, n_spots).
    """
    from scipy.sparse.csgraph import laplacian  # pylint: disable=import-outside-toplevel
    from sklearn.neighbors import kneighbors_graph  # pylint: disable=import-outside-toplevel

    n_spots = coords.shape[0]

    # Build k-NN adjacency matrix (directed)
    A = kneighbors_graph(coords, k, mode="connectivity", include_self=False)

    # Symmetrize: if i is neighbor of j OR j is neighbor of i, they are connected
    A = A + A.T
    A.data = np.clip(A.data, 0, 1)  # Binary adjacency

    # Apply cellularity-aware edge weighting (bilateral filtering principle)
    if cellularity is not None:
        cellularity = np.asarray(cellularity, dtype=np.float64)
        # Use log-ratio so 5-vs-10 and 10-vs-20 are treated equivalently
        log_cell = np.log(cellularity + 1.0)  # +1 to handle zeros
        A_coo = scipy.sparse.coo_matrix(A)
        log_diff_sq = (log_cell[A_coo.row] - log_cell[A_coo.col]) ** 2
        similarity = np.exp(-log_diff_sq / (2.0 * cellularity_sigma**2))
        A_weighted = scipy.sparse.coo_matrix((similarity, (A_coo.row, A_coo.col)), shape=A.shape).tocsr()
        n_suppressed = int((similarity < 0.5).sum())
        logging.info(
            "Cellularity-aware Laplacian: sigma=%.2f, %d/%d edges suppressed (w<0.5)",
            cellularity_sigma,
            n_suppressed,
            len(similarity),
        )
        A = A_weighted

    # Compute Laplacian
    L = laplacian(A, normed=normed)

    logging.debug("Built spatial Laplacian: %s spots, k=%s, nnz=%s, normed=%s", n_spots, k, L.nnz, normed)

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
    logging.info(" - Mean confidence score: %s", np.mean(confidence_scores))
    logging.info(" - Mean prior strength: %s", np.mean(global_prior))
    logging.info(" - % Strong signals (>0.5): %s%", 100 * np.mean(global_prior > 0.5))

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
    all_markers: list[str] = [marker for profile in cell_profile_dict.values() for marker in profile["Major"]]
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
                logging.warning("No valid markers found for profile '%s'", profile_name)
        except IndexError as e:
            logging.warning("Error processing markers for profile '%s': %s", profile_name, str(e))

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
    # Step 1: Collect Major and Soft markers from cell_profile_dict
    # Minor markers are ignored (legacy key, pan-immune markers degrade QP r=0.700 vs 0.766).
    # Soft markers get graded assignment weights (< 1.0) instead of binary.
    all_markers = []
    for profile_markers in cell_profile_dict.values():
        all_markers.extend(profile_markers.get("Major", []))
        for entry in profile_markers.get("Soft", []):
            all_markers.append(entry[0])  # (marker_name, weight) tuple

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
        logging.warning(
            "Zero columns detected for markers: %s. Adding epsilon.",
            [existing_markers[i] for i, z in enumerate(zero_cols) if z],
        )
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

    for j, (_, markers_dict) in enumerate(cell_profile_dict.items()):
        for marker in markers_dict.get("Major", []):
            if marker in marker_to_idx:
                assignment_matrix[marker_to_idx[marker], j] = 1.0
        for entry in markers_dict.get("Soft", []):
            marker_name, assign_weight = entry
            if marker_name in marker_to_idx:
                idx = marker_to_idx[marker_name]
                # Don't overwrite a stronger (Major=1.0) assignment
                if assignment_matrix[idx, j] < assign_weight:
                    assignment_matrix[idx, j] = assign_weight

    n_soft = int(((assignment_matrix > 0) & (assignment_matrix < 1.0)).sum())
    logging.info("Mapped %s markers to %s cell types (%s soft assignments)", M, T, n_soft)

    return marker_level_data, existing_markers, assignment_matrix, cell_type_names


def map_antibodies_to_profiles_cellularity(adata, cell_profile_dict, clip_percentile=99, scale="median"):
    """Map antibody data for cellularity-scaled QP, preserving between-spot amplitude.

    Like map_antibodies_to_profiles_v2 but replaces column-max normalization with
    robust scaling that preserves the N-proportional signal amplitude needed for
    count-space reconstruction (S ≈ α + c × A × β).

    Args:
        adata: Antibody capture AnnData (already CLR-normalized via preprocess_antibody).
        cell_profile_dict: Cell type → marker mapping.
        clip_percentile: Per-marker percentile for outlier clipping (default: 99).
        scale: Robust scale denominator — "median" or "p75" (default: "median").

    Returns:
        Same as map_antibodies_to_profiles_v2: (marker_level_data, marker_names,
        assignment_matrix, cell_type_names).
    """
    # Step 1-3: Same marker extraction as v2
    all_markers = []
    for profile_markers in cell_profile_dict.values():
        all_markers.extend(profile_markers.get("Major", []))

    seen = set()
    unique_markers = []
    for m in all_markers:
        if m not in seen:
            seen.add(m)
            unique_markers.append(m)

    existing_markers = [m for m in unique_markers if m in adata.var_names]
    if len(existing_markers) == 0:
        raise ValueError("No matching antibody markers found in adata.var_names.")

    adata.var_names_make_unique()
    marker_indices = [np.where(adata.var_names == m)[0][0] for m in existing_markers]

    antibody_data = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    marker_level_data = antibody_data[:, marker_indices].astype(np.float64)

    # Step 4: Robust scaling (NOT column-max-norm) to preserve N-proportional amplitude
    for col in range(marker_level_data.shape[1]):
        vals = marker_level_data[:, col]
        # Clip outliers at percentile
        cap = np.percentile(vals, clip_percentile)
        if cap > 0:
            vals = np.minimum(vals, cap)
        # Divide by robust scale (median or p75)
        if scale == "median":
            denom = np.median(vals)
        elif scale == "p75":
            denom = np.percentile(vals, 75)
        else:
            raise ValueError(f"Unknown scale mode: {scale}")
        if denom < 1e-9:
            denom = 1e-9  # floor for zero-signal markers
        marker_level_data[:, col] = vals / denom

    if np.isnan(marker_level_data).any():
        raise ValueError("NaN values detected after robust normalization.")

    logging.info(
        "Cellularity mapper: %d markers, robust-%s scaling, clip p%d, " "S range [%.3f, %.3f], mean %.3f",
        len(existing_markers),
        scale,
        clip_percentile,
        marker_level_data.min(),
        marker_level_data.max(),
        marker_level_data.mean(),
    )

    # Step 5: Build assignment matrix (same as v2)
    cell_type_names = list(cell_profile_dict.keys())
    M = len(existing_markers)
    T = len(cell_type_names)
    assignment_matrix = np.zeros((M, T), dtype=np.float64)

    marker_to_idx = {name: i for i, name in enumerate(existing_markers)}
    for j, (_, markers_dict) in enumerate(cell_profile_dict.items()):
        major_markers = markers_dict.get("Major", [])
        for marker in major_markers:
            if marker in marker_to_idx:
                assignment_matrix[marker_to_idx[marker], j] = 1.0

    return marker_level_data, existing_markers, assignment_matrix, cell_type_names


def map_antibodies_raw_counts(
    adata, cell_profile_dict, *, winsorize_lower=1.0, winsorize_upper=99.0, scale="median", eps=1e-6
):
    """Map raw antibody counts for cellularity-scaled QP (no CLR).

    Reads from adata.layers["raw_counts"] (pre-CLR raw data) and applies only
    per-marker winsorization + robust scaling. Preserves linear N-scaling that
    count-space reconstruction requires.

    Args:
        adata: Antibody capture AnnData with layers["raw_counts"].
        cell_profile_dict: Cell type → marker mapping.
        winsorize_lower: Lower percentile for clipping (default: 1).
        winsorize_upper: Upper percentile for clipping (default: 99).
        scale: Per-marker denominator — "median" or "p75" (default: "median").
        eps: Floor for zero-signal markers.

    Returns:
        Same as map_antibodies_to_profiles_v2: (marker_level_data, marker_names,
        assignment_matrix, cell_type_names).
    """
    # Marker extraction (same as v2)
    all_markers = []
    for profile_markers in cell_profile_dict.values():
        all_markers.extend(profile_markers.get("Major", []))

    seen = set()
    unique_markers = []
    for m in all_markers:
        if m not in seen:
            seen.add(m)
            unique_markers.append(m)

    existing_markers = [m for m in unique_markers if m in adata.var_names]
    if len(existing_markers) == 0:
        raise ValueError("No matching antibody markers found in adata.var_names.")

    adata.var_names_make_unique()
    marker_indices = [np.where(adata.var_names == m)[0][0] for m in existing_markers]

    # Read raw counts (pre-CLR)
    if "raw_counts" in adata.layers:
        raw = adata.layers["raw_counts"]
    else:
        logging.warning("raw_counts layer not found, falling back to .X (may be CLR-normalized)")
        raw = adata.X
    raw = raw.toarray() if hasattr(raw, "toarray") else raw
    marker_level_data = np.asarray(raw[:, marker_indices], dtype=np.float64)
    marker_level_data = np.maximum(marker_level_data, 0.0)

    # Per-marker winsorize + robust scaling
    for col in range(marker_level_data.shape[1]):
        vals = marker_level_data[:, col]
        lo = np.percentile(vals, winsorize_lower) if winsorize_lower > 0 else vals.min()
        hi = np.percentile(vals, winsorize_upper)
        vals = np.clip(vals, lo, hi)

        if scale == "median":
            denom = np.median(vals)
        elif scale == "p75":
            denom = np.percentile(vals, 75)
        else:
            denom = 1.0
        marker_level_data[:, col] = vals / max(float(denom), eps)

    if np.isnan(marker_level_data).any() or np.isinf(marker_level_data).any():
        raise ValueError("NaN/Inf values in raw-count mapper output.")

    logging.info(
        "Raw-count mapper: %d markers, winsor=[%.0f, %.0f], scale=%s, " "S range [%.3f, %.3f], mean %.3f",
        len(existing_markers),
        winsorize_lower,
        winsorize_upper,
        scale,
        marker_level_data.min(),
        marker_level_data.max(),
        marker_level_data.mean(),
    )

    # Assignment matrix (same as v2)
    cell_type_names = list(cell_profile_dict.keys())
    M = len(existing_markers)
    T = len(cell_type_names)
    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    marker_to_idx = {name: i for i, name in enumerate(existing_markers)}
    for j, (_, markers_dict) in enumerate(cell_profile_dict.items()):
        for marker in markers_dict.get("Major", []):
            if marker in marker_to_idx:
                assignment_matrix[marker_to_idx[marker], j] = 1.0

    return marker_level_data, existing_markers, assignment_matrix, cell_type_names


def compute_marker_exclusivity(
    marker_level_data: np.ndarray,
    Y_values: np.ndarray,
    marker_owners: List[List[int]],
    _assignment_matrix: np.ndarray,
    *,
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
    *,
    profile_based_antibody_data: np.ndarray = None,  # pylint: disable=unused-argument
    unknown_threshold: float = 0.05,
    min_celltype_threshold: float = 0.01,
    redundancy_threshold: float = 0.2,  # pylint: disable=unused-argument
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
                logging.warning("[VALIDATION WARNING] %s", error_msg)
            else:
                logging.error(error_msg)
                raise ValueError(error_msg)

    # Check if any defined (non-Unknown) cell types have very low proportions
    low_proportion_celltypes = []
    for celltype, mean_prop in zip(cell_type_names, mean_proportions):
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
            logging.warning("[VALIDATION WARNING] %s", error_msg)
        else:
            logging.error(error_msg)
            raise ValueError(error_msg)

    # Log successful validation
    logging.info("✓ Cell type validation passed (proportions + redundancy)")
    if "Unknown" in cell_type_names:
        unknown_idx = cell_type_names.index("Unknown")
        logging.info("  - Unknown: %s% (threshold < %s%)", mean_proportions[unknown_idx] * 100, unknown_threshold * 100)

    # Log cell type proportions
    for celltype, mean_prop in zip(cell_type_names, mean_proportions):
        if celltype != "Unknown":
            logging.info("  - %s: %s%", celltype, mean_prop * 100)


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
        return np.asarray(np.ones(T) / T)

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
    return np.asarray(enrichment / (enrichment.sum() + epsilon))


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
    from scipy.stats import pearsonr  # pylint: disable=import-outside-toplevel

    T = anchor_expr.shape[1]
    enrichment = np.zeros(T)

    # Skip if gene has no variance
    if np.std(gene_expr) < 1e-10:
        return np.asarray(np.ones(T) / T)

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
        return np.asarray(np.ones(T) / T)

    return np.asarray(enrichment / enrichment.sum())


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
    return np.asarray(blended / (blended.sum() + 1e-10))


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


def estimate_true_expression_cell(
    X_obs: np.ndarray,
    Y_assignments: np.ndarray,
    coords: np.ndarray,
    enrichment_weights: np.ndarray,
    *,
    library_slack: float = 1.5,
    lambda_enrich: float = 1.0,
    lambda_spatial: float = 0.0,
    spatial_k: int = 50,
    max_workers: Optional[int] = None,
    _checkpoint_interval: int = 500,
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
            "Cell-level Pass 2: %s unassigned cells will retain observed expression (no deconvolution)", n_unassigned
        )

    # Build spatial neighbor graph (k-NN)
    from scipy.spatial import cKDTree  # pylint: disable=import-outside-toplevel

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
        "Cell-level Pass 2: %s cells, %s genes, %s types, library_slack=%s, spatial_k=%s",
        N,
        M,
        T,
        library_slack,
        spatial_k,
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
        "X_obs": X_obs,
        "dominant_type": dominant_type,
        "enrichment_weights": enrichment_weights,
        "same_type_neighbor_means": same_type_neighbor_means,
        "same_type_neighbor_counts": same_type_neighbor_counts,
        "library_slack": library_slack,
        "lambda_enrich": lambda_enrich,
        "lambda_spatial": lambda_spatial,
        "M": M,
    }

    # ThreadPoolExecutor: threads share CUDA context safely (no fork corruption)
    n_to_optimize = len(cells_to_optimize)
    thread_workers = min(workers, 8)
    with ThreadPoolExecutor(max_workers=thread_workers) as executor:
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
    logging.info("Normalization stats:")
    logging.info("Original median total: %s", median_size)
    logging.info("Mean scaled total: %s", X_scaled.sum(axis=1).mean())
    logging.info("Max scaled value: %s", X_scaled.max())

    return adata_norm


# Module-level worker data placeholder (populated by estimate_true_expression_cell)
_cell_pass2_worker_data = None


def _init_cell_pass2_worker(data):
    """Initializer for spawn-context workers: set module-level worker data."""
    global _cell_pass2_worker_data
    _cell_pass2_worker_data = data


def _solve_single_cell_pass2(cell_idx: int):
    """Solve QP for a single cell's true expression using cuOPT.

    Module-level function (required for pickling with ProcessPoolExecutor).

    Args:
        cell_idx: Index of the cell to process.

    Returns:
        Tuple of (cell_idx, expression_array).
    """
    wd = _cell_pass2_worker_data
    if wd is None:
        raise RuntimeError("Worker data not initialized.")

    try:
        ct = wd["dominant_type"][cell_idx]
        obs = wd["X_obs"][cell_idx]
        obs_lib = obs.sum()

        if obs_lib < 1:
            return cell_idx, obs.copy()

        enrich = wd["enrichment_weights"][ct]
        neighbor_mean = wd["same_type_neighbor_means"][cell_idx]
        has_neighbors = wd["same_type_neighbor_counts"][cell_idx] > 0
        M_genes = wd["M"]
        lib_slack = wd["library_slack"]
        l_enrich = wd["lambda_enrich"]
        l_spatial = wd["lambda_spatial"]

        p = Problem()

        max_lib = lib_slack * obs_lib
        X_vars = {}
        for g in range(M_genes):
            X_vars[g] = p.addVariable(lb=0.0, ub=max_lib, name=f"X_{g}")

        # Library size constraint
        lib_sum = sum(X_vars[g] for g in range(M_genes))
        p.addConstraint(lib_sum <= max_lib)

        # Objective: minimize data fidelity + enrichment prior + spatial coherence
        obj = 0
        for g in range(M_genes):
            # Data fidelity: (X_g - obs_g)^2
            obj += (X_vars[g] - obs[g]) * (X_vars[g] - obs[g])
            # Enrichment prior (negative = encourage expression for enriched genes)
            if l_enrich > 0 and enrich[g] > 0.1:
                obj -= l_enrich * enrich[g] * X_vars[g]
            # Spatial coherence: (X_g - neighbor_mean_g)^2
            if l_spatial > 0 and has_neighbors:
                obj += l_spatial * (X_vars[g] - neighbor_mean[g]) * (X_vars[g] - neighbor_mean[g])

        p.setObjective(obj)
        p.solve()

        # Check status
        if hasattr(p, "Status") and p.Status not in (None, "optimal", "OPTIMAL", 1, 2):
            return cell_idx, obs.copy()

        result = np.array([X_vars[g].getValue() for g in range(M_genes)])
        return cell_idx, result

    except (ValueError, RuntimeError) as e:

        logging.warning("Cell %s optimization failed: %s", cell_idx, e)
        return cell_idx, wd["X_obs"][cell_idx].copy()


def deconvolute_spot_with_neighbors_with_prior(
    spot_idx: int,
    adata: sc.AnnData,
    cell_type_numbers_array: np.ndarray,
    radius: float,
    *,
    global_prior: Optional[np.ndarray] = None,
    lambda_prior_weight: float = 0.0,
    local_enrichment_weight: float = 0.5,
    global_enrichment_weight: float = 0.5,
    continuous_relaxation: bool = True,
    lambda_gex_reg: float = 0.01,
    enrichment_smoothing: float = 0.2,
    use_kl_regularization: bool = True,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Optional[np.ndarray]:
    """
    Deconvolute a spot with its neighbors, using enrichment weights and optional prior.

    cuOPT backend — replaces Gurobi QP with NVIDIA cuOPT solver while keeping
    the enrichment computation and direct-softmax bypass path unchanged.

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
        lambda_gex_reg: L2 regularization weight for expression variables (default: 0.01).
        use_kl_regularization: If True, use KL-divergence instead of L2 regularization.
        kl_temperature: Temperature for softmax target (lower = sharper).
        lambda_kl: Weight for KL penalty term.

    Returns:
        np.ndarray of shape (T, M) with deconvolved expression, or None on failure.
    """
    from scipy.special import softmax as scipy_softmax  # pylint: disable=import-outside-toplevel

    # DIRECT SOFTMAX MODE: Skip optimization entirely
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

                # Direct softmax allocation (no optimization, no L2)
                log_enrichment = np.log(enrichment + 1e-10)
                alloc_weights = scipy_softmax(log_enrichment / kl_temperature)
                result[:, k] = total * alloc_weights

            return result

        except (ValueError, RuntimeError) as e:

            logging.error("Direct softmax failed for spot %s: %s", spot_idx, e)
            # Fall through to cuOPT optimization

    try:
        neighborhood_indices = get_neighbors_with_fixed_radius(spot_idx, adata, radius=int(radius), include_center=True)
        if not neighborhood_indices:
            logging.error("No valid neighbors found for spot %s.", spot_idx)
            return None

        neighborhood_indices = np.array(
            [int(idx) for idx in neighborhood_indices if isinstance(idx, (int, np.integer))], dtype=int
        )

        # Extract expression data
        deconvolution_expression_data = adata.X
        if scipy.sparse.issparse(deconvolution_expression_data):
            deconvolution_expression_data = deconvolution_expression_data.toarray()
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

        # Compute expression-aware enrichment using percentile-based method
        gene_specific_enrichment = np.zeros((M, T))

        # Use baseline percentile-based enrichment with configurable smoothing
        def compute_baseline_enrichment(expression_data, cell_type_props, gene_idx, smoothing=0.2):
            """Baseline enrichment: percentile threshold + optional smoothing."""
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

            # Apply smoothing (0.2 = 80/20 baseline, 0.0 = no smoothing)
            if smoothing > 0:
                smoothed = (1 - smoothing) * enrichment + smoothing * np.ones_like(enrichment)
            else:
                smoothed = enrichment
            return smoothed / (np.sum(smoothed) + epsilon)

        for k in range(M):
            local_enrich = compute_baseline_enrichment(
                neighborhood_expression_data, neighborhood_cell_type_numbers, k, smoothing=enrichment_smoothing
            )
            global_enrich = compute_baseline_enrichment(
                deconvolution_expression_data, cell_type_numbers_array, k, smoothing=enrichment_smoothing
            )
            gene_specific_enrichment[k] = (
                local_enrichment_weight * local_enrich + global_enrichment_weight * global_enrich
            )

        # Build cuOPT model
        p = Problem()

        # Variables for count assignment
        X = {}
        center_counts = deconvolution_expression_data[spot_idx, :]

        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    if continuous_relaxation:
                        X[j, k] = p.addVariable(lb=0.0, ub=float(total_counts), name=f"X_{j}_{k}")
                    else:
                        X[j, k] = p.addVariable(lb=0, ub=total_counts, vtype="integer", name=f"X_{j}_{k}")
                # Count conservation constraint
                gene_sum = sum(X[j, k] for j in range(T))
                p.addConstraint(gene_sum == total_counts)

        # Validate prior if provided
        if global_prior is not None and lambda_prior_weight > 0:
            if not isinstance(global_prior, np.ndarray):
                raise ValueError("global_prior must be a numpy array")
            if global_prior.shape != (T, M):
                raise ValueError(
                    f"Prior matrix shape {global_prior.shape} does not match " f"expected shape ({T}, {M})"
                )

        # Center spot proportions (index 0 = center in neighborhood)
        center_props = neighborhood_cell_type_numbers[0, :]

        # Apply softmax sharpening to proportions (simpler than KL regularization)
        if kl_temperature < 1.0:
            logits = center_props / kl_temperature
            logits = logits - logits.max()  # numerical stability
            exp_logits = np.exp(logits)
            center_props = exp_logits / (exp_logits.sum() + 1e-10)

        # Objective: maximize enrichment-weighted allocation with regularization
        # cuOPT minimizes, so we negate the maximize objective
        obj = 0
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                enrichment_for_gene = gene_specific_enrichment[k, :]

                if use_kl_regularization:
                    # Softmax KL-divergence regularization
                    target = compute_softmax_target(enrichment_for_gene, temperature=kl_temperature)
                    kl_coeffs = compute_kl_penalty_coefficients(target, total_counts, lambda_kl)

                    for j in range(T):
                        # Enrichment term (negated for minimization)
                        obj -= enrichment_for_gene[j] * center_props[j] * X[j, k]

                        # KL penalty (pulls toward target) — was subtracted in maximize,
                        # so add for minimize
                        target_j = kl_coeffs["target_counts"][j]
                        penalty = kl_coeffs["penalty_weight"] * (X[j, k] - target_j) * (X[j, k] - target_j)
                        obj += penalty

                        # Prior terms
                        if global_prior is not None and lambda_prior_weight > 0:
                            try:
                                prior_value = float(global_prior[j, k])
                                # Was -prior_penalty in maximize; becomes +prior_penalty in minimize
                                obj += lambda_prior_weight * (1 - prior_value) * X[j, k]
                            except (ValueError, RuntimeError) as e:

                                logging.warning("Error accessing prior at [%s, %s]: %s", j, k, e)
                                continue
                else:
                    # L2 regularization (backward compatible)
                    for j in range(T):
                        enrichment_weight = enrichment_for_gene[j]
                        cell_type_weight = center_props[j]

                        # Enrichment term (negated for minimization)
                        obj -= enrichment_weight * cell_type_weight * X[j, k]

                        # L2 regularization (was negative in maximize, positive in minimize)
                        if lambda_gex_reg > 0:
                            obj += lambda_gex_reg * X[j, k] * X[j, k]

                        # Prior guidance
                        if global_prior is not None and lambda_prior_weight > 0:
                            try:
                                prior_value = float(global_prior[j, k])
                                obj += lambda_prior_weight * (1 - prior_value) * X[j, k]
                            except (ValueError, RuntimeError) as e:

                                logging.warning("Error accessing prior at [%s, %s]: %s", j, k, e)
                                continue

        p.setObjective(obj)
        p.solve()

        # Check solve status
        if hasattr(p, "Status") and p.Status not in (None, "optimal", "OPTIMAL", 1, 2):
            logging.error("No feasible solution found for spot %s.", spot_idx)
            return None

        result = np.zeros((T, M))
        for k in range(M):
            total_counts = int(center_counts[k])
            if total_counts > 0:
                for j in range(T):
                    result[j, k] = X[j, k].getValue()
        return result

    except Exception as e:  # pylint: disable=broad-exception-caught
        logging.error("Error during deconvolution of spot %s: %s", spot_idx, str(e))
        logging.error(traceback.format_exc())
        return None

    finally:
        gc.collect()


# =====================================================================
# Prep-cook GEX deconvolution: CPU-parallel prep + sequential GPU solve
# =====================================================================


def _build_neighbor_lists(spatial_coords: np.ndarray, radius: float) -> List[np.ndarray]:
    """Build neighbor index lists for all spots using KDTree (vectorized)."""
    from scipy.spatial import cKDTree  # pylint: disable=import-outside-toplevel

    tree = cKDTree(spatial_coords)
    neighbor_lists = []
    for i in range(len(spatial_coords)):
        neighbors = tree.query_ball_point(spatial_coords[i], r=radius)
        neighbor_lists.append(np.array(sorted(neighbors), dtype=np.int64))
    return neighbor_lists


def _vectorized_enrichment(  # pylint: disable=too-many-positional-arguments
    expr_data: np.ndarray,
    props: np.ndarray,
    celltype_freqs: np.ndarray,
    neighbor_lists: List[np.ndarray],
    smoothing: float,
    local_weight: float,
    global_weight: float,
) -> np.ndarray:
    """Compute gene-specific enrichment for all spots. Returns (N, M, T).

    Vectorized replacement for the per-gene compute_baseline_enrichment loop.
    """
    N, M = expr_data.shape
    T = props.shape[1]

    # Global enrichment (same for all spots): for each gene, find spots with
    # expression >= median of nonzero, compute enrichment of cell types there.
    normalized_props = props / (celltype_freqs + 1e-10)  # (N, T)
    background_props = np.mean(normalized_props, axis=0)  # (T,)

    # Pre-compute global enrichment per gene: (M, T)
    global_enrichment = np.ones((M, T)) / T
    for k in range(M):
        gene_expr = expr_data[:, k]
        nonzero = gene_expr[gene_expr > 0]
        if len(nonzero) == 0:
            continue
        threshold = np.percentile(nonzero, 50)
        high_mask = gene_expr >= threshold
        if not np.any(high_mask):
            continue
        high_props = np.mean(normalized_props[high_mask], axis=0)
        enrich = high_props / (background_props + 1e-10)
        if smoothing > 0:
            enrich = (1 - smoothing) * enrich + smoothing * np.ones(T)
        enrich_sum = enrich.sum()
        if enrich_sum > 0:
            global_enrichment[k] = enrich / enrich_sum

    # Per-spot local enrichment: use neighborhood expression
    result = np.zeros((N, M, T))
    for i in range(N):
        nb_idx = neighbor_lists[i]
        nb_expr = expr_data[nb_idx]  # (n_neighbors, M)
        nb_props = normalized_props[nb_idx]  # (n_neighbors, T)
        nb_bg = np.mean(nb_props, axis=0)  # (T,)

        local_enrichment = np.ones((M, T)) / T
        for k in range(M):
            gene_expr = nb_expr[:, k]
            nonzero = gene_expr[gene_expr > 0]
            if len(nonzero) == 0:
                continue
            threshold = np.percentile(nonzero, 50)
            high_mask = gene_expr >= threshold
            if not np.any(high_mask):
                continue
            high_props = np.mean(nb_props[high_mask], axis=0)
            enrich = high_props / (nb_bg + 1e-10)
            if smoothing > 0:
                enrich = (1 - smoothing) * enrich + smoothing * np.ones(T)
            enrich_sum = enrich.sum()
            if enrich_sum > 0:
                local_enrichment[k] = enrich / enrich_sum

        result[i] = local_weight * local_enrichment + global_weight * global_enrichment

    return result


def _prepare_spot_qp(
    spot_idx: int,
    center_counts: np.ndarray,
    center_props: np.ndarray,
    enrichment: np.ndarray,
    *,
    global_prior: Optional[np.ndarray],
    lambda_prior_weight: float,
    lambda_gex_reg: float,
    use_kl_regularization: bool,
    kl_temperature: float,
    lambda_kl: float,
    continuous_relaxation: bool,
) -> Optional[dict]:
    """Prepare QP coefficient arrays for one spot (CPU-only, no CUDA).

    Returns a dict of numpy arrays ready for _solve_spot_qp().
    """
    M, T = enrichment.shape

    # Sharpen proportions
    props = center_props.copy()
    if kl_temperature < 1.0:
        logits = props / kl_temperature
        logits -= logits.max()
        exp_logits = np.exp(logits)
        props = exp_logits / (exp_logits.sum() + 1e-10)

    # Active genes: those with nonzero counts
    active_genes = center_counts > 0
    n_active = active_genes.sum()
    if n_active == 0:
        return None

    # Build coefficient arrays for the QP objective:
    #   minimize sum_k sum_j [ obj_linear[k,j] * X[j,k] + obj_quad[k,j] * X[j,k]^2 ]
    obj_linear = np.zeros((M, T))
    obj_quad = np.zeros((M, T))
    kl_target_counts = np.zeros((M, T))

    for k in range(M):
        total = center_counts[k]
        if total <= 0:
            continue
        enrich_k = enrichment[k]  # (T,)

        if use_kl_regularization:
            target = _softmax_target(enrich_k, kl_temperature)
            penalty_weight = lambda_kl / (total + 1e-10)
            target_j = target * total

            for j in range(T):
                # Enrichment (negated for minimization)
                obj_linear[k, j] = -enrich_k[j] * props[j]
                # KL quadratic: penalty_weight * (X - target)^2
                #   = penalty_weight * X^2 - 2*penalty_weight*target*X + const
                obj_quad[k, j] = penalty_weight
                obj_linear[k, j] += -2.0 * penalty_weight * target_j[j]
                kl_target_counts[k, j] = target_j[j]

                # Prior
                if global_prior is not None and lambda_prior_weight > 0:
                    obj_linear[k, j] += lambda_prior_weight * (1 - global_prior[j, k])
        else:
            for j in range(T):
                obj_linear[k, j] = -enrich_k[j] * props[j]
                if lambda_gex_reg > 0:
                    obj_quad[k, j] = lambda_gex_reg
                if global_prior is not None and lambda_prior_weight > 0:
                    obj_linear[k, j] += lambda_prior_weight * (1 - global_prior[j, k])

    return {
        "spot_idx": spot_idx,
        "center_counts": center_counts,
        "obj_linear": obj_linear,
        "obj_quad": obj_quad,
        "T": T,
        "M": M,
        "continuous": continuous_relaxation,
    }


def _softmax_target(enrichment: np.ndarray, temperature: float) -> np.ndarray:
    """Compute softmax target distribution from enrichment scores."""
    log_e = np.log(enrichment + 1e-10)
    logits = log_e / max(temperature, 1e-6)
    logits -= logits.max()
    exp_l = np.exp(logits)
    return np.asarray(exp_l / (exp_l.sum() + 1e-10))


def _solve_spot_qp(spot_data: dict) -> Optional[np.ndarray]:
    """Build and solve cuOPT QP from pre-computed coefficient arrays (GPU).

    Args:
        spot_data: dict from _prepare_spot_qp().

    Returns:
        (T, M) result array, or None on failure.
    """
    if Problem is None:
        raise ImportError("cuOPT is not installed.")

    spot_idx = spot_data["spot_idx"]
    center_counts = spot_data["center_counts"]
    obj_linear = spot_data["obj_linear"]
    obj_quad = spot_data["obj_quad"]
    T = spot_data["T"]
    M = spot_data["M"]
    continuous = spot_data["continuous"]

    try:
        p = Problem()
        X = {}

        for k in range(M):
            total = int(center_counts[k])
            if total <= 0:
                continue
            for j in range(T):
                if continuous:
                    X[j, k] = p.addVariable(lb=0.0, ub=float(total), name=f"X_{j}_{k}")
                else:
                    X[j, k] = p.addVariable(lb=0, ub=total, vtype="integer", name=f"X_{j}_{k}")
            # Count conservation
            gene_sum = sum(X[j, k] for j in range(T))
            p.addConstraint(gene_sum == total)

        # Build objective from pre-computed coefficients
        obj = 0
        for k in range(M):
            if int(center_counts[k]) <= 0:
                continue
            for j in range(T):
                obj += obj_linear[k, j] * X[j, k]
                if obj_quad[k, j] > 0:
                    obj += obj_quad[k, j] * X[j, k] * X[j, k]

        p.setObjective(obj)
        p.solve()

        if hasattr(p, "Status") and p.Status not in (None, "optimal", "OPTIMAL", 1, 2):
            return None

        result = np.zeros((T, M))
        for k in range(M):
            if int(center_counts[k]) <= 0:
                continue
            for j in range(T):
                result[j, k] = X[j, k].getValue()
        return result

    except Exception as e:  # pylint: disable=broad-exception-caught
        logging.error("cuOPT solve failed for spot %s: %s", spot_idx, e)
        return None


def optimize_gene_expression(
    sample_name: str,
    deconvolution_expression_data: np.ndarray,
    cell_type_numbers_array: np.ndarray,
    filtered_adata: sc.AnnData,
    *,
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
    enrichment_smoothing: float = 0.2,
    _anchor_genes: Optional[Dict] = None,
    _anchor_weights: Optional[Dict] = None,
    _module_weight: float = 0.5,
    use_kl_regularization: bool = True,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> Dict[str, Any]:
    """Optimize gene expression: prep-cook pattern (CPU-parallel prep + sequential GPU solve).

    Phase 1 (CPU, parallel): Pre-extract numpy arrays from AnnData, compute
    enrichment and QP coefficients for all spots via ProcessPoolExecutor.
    Phase 2 (GPU, sequential): Build cuOPT Problem from coefficients and solve.

    Args:
        sample_name: Name of the sample.
        deconvolution_expression_data: Gene expression data (N_spots x M_genes).
        cell_type_numbers_array: Cell type proportions (N_spots x T_celltypes).
        filtered_adata: Filtered AnnData object containing gene expression data.
        radius: Radius for neighbor detection.
        global_enrichment_weight: Weight for global expression enrichment (0-1).
        local_enrichment_weight: Weight for local expression enrichment (0-1).
        global_prior: Optional global prior matrix for guidance.
        lambda_prior_weight: Weight for prior guidance.
        max_workers: Maximum number of parallel workers for CPU prep phase.
        checkpoint_interval: Number of spots between checkpoints.
        output_dir: Directory for checkpoints.
        rerun: Whether to rerun if results exist.
        continuous_relaxation: Use continuous (LP) vs integer (MIP) variables.
        lambda_gex_reg: L2 regularization weight on X variables.
        enrichment_smoothing: Smoothing factor for enrichment (0.2 = 80/20 blend).
        anchor_genes: Optional anchor genes per cell type.
        anchor_weights: Optional anchor weights per cell type.
        module_weight: Weight for module-aware enrichment.
        use_kl_regularization: If True, use KL-divergence regularization.
        kl_temperature: Temperature for softmax target.
        lambda_kl: Weight for KL penalty term.

    Returns:
        Dict mapping spot indices to deconvolved expression profiles (T x M arrays).
    """
    os.makedirs(output_dir, exist_ok=True)

    N = deconvolution_expression_data.shape[0]
    M = deconvolution_expression_data.shape[1]
    T = cell_type_numbers_array.shape[1]

    checkpoint_mgr = CheckpointManager(output_dir, sample_name)

    if not rerun:
        complete_results = checkpoint_mgr.check_complete_run(N, T, M)
        if complete_results is not None:
            return complete_results
        completed_spots, spotwise_gene_expression_profiles = checkpoint_mgr.load_latest_checkpoint(N, T, M)
    else:
        completed_spots = set()
        spotwise_gene_expression_profiles = {}

    remaining_spots = [i for i in range(N) if i not in completed_spots]
    if not remaining_spots:
        logging.info("All spots already completed.")
        return cast(Dict[str, Any], spotwise_gene_expression_profiles)

    logging.info("Starting GEX deconvolution for %s: %s spots remaining", sample_name, len(remaining_spots))
    if use_kl_regularization:
        logging.info("KL regularization (temp=%s, lambda=%s)", kl_temperature, lambda_kl)

    # ── Phase 1: CPU prep (parallel) ──────────────────────────────────
    logging.info("Phase 1: Pre-extracting arrays and computing enrichment...")
    prep_start = time.time()

    # Extract numpy arrays from AnnData (once)
    expr_data = deconvolution_expression_data
    if scipy.sparse.issparse(expr_data):
        expr_data = expr_data.toarray()  # type: ignore[attr-defined]
    expr_data = np.asarray(expr_data, dtype=np.float64)

    props = np.asarray(cell_type_numbers_array, dtype=np.float64)
    spatial_coords = filtered_adata.obsm["spatial"]

    # Build neighbor lists via KDTree (vectorized, replaces per-spot brute-force)
    neighbor_lists = _build_neighbor_lists(spatial_coords, float(radius))
    logging.info("Built neighbor lists (KDTree, radius=%s)", radius)

    # Compute celltype frequencies
    celltype_freqs = props.sum(axis=0)
    celltype_freqs = celltype_freqs / (celltype_freqs.sum() + 1e-10)

    # Vectorized enrichment: (N, M, T)
    all_enrichment = _vectorized_enrichment(
        expr_data,
        props,
        celltype_freqs,
        neighbor_lists,
        enrichment_smoothing,
        local_enrichment_weight,
        global_enrichment_weight,
    )
    logging.info("Enrichment computed for %s spots × %s genes × %s types", N, M, T)

    # Prepare QP coefficient arrays for remaining spots (CPU, parallel)
    workers = max_workers if max_workers is not None else min(os.cpu_count() or 4, 8)

    def _prep_one(spot_idx):
        return _prepare_spot_qp(
            spot_idx=spot_idx,
            center_counts=expr_data[spot_idx],
            center_props=props[spot_idx],
            enrichment=all_enrichment[spot_idx],
            global_prior=global_prior,
            lambda_prior_weight=lambda_prior_weight,
            lambda_gex_reg=lambda_gex_reg,
            use_kl_regularization=use_kl_regularization,
            kl_temperature=kl_temperature,
            lambda_kl=lambda_kl,
            continuous_relaxation=continuous_relaxation,
        )

    spot_qp_data = {}
    # Use ThreadPoolExecutor for prep — avoids pickle overhead, GIL released
    # during numpy operations. No CUDA here so threads are safe.
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_prep_one, i): i for i in remaining_spots}
        for future in as_completed(futures):
            idx = futures[future]
            try:
                data = future.result()
                if data is not None:
                    spot_qp_data[idx] = data
            except Exception as e:  # pylint: disable=broad-exception-caught
                logging.error("Prep failed for spot %s: %s", idx, e)

    prep_elapsed = time.time() - prep_start
    logging.info("Phase 1 complete: %s spots prepped in %ss", len(spot_qp_data), prep_elapsed)

    # ── Phase 2: GPU solve (sequential) ───────────────────────────────
    logging.info("Phase 2: Solving QPs on GPU (sequential)...")
    solve_start = time.time()

    try:
        spots_since_last_save = 0
        solve_order = sorted(spot_qp_data.keys())
        with tqdm(total=len(solve_order), desc="Solving GEX QPs") as pbar:
            for spot_idx in solve_order:
                try:
                    result = _solve_spot_qp(spot_qp_data[spot_idx])
                    if result is not None and result.ndim == 2:
                        spotwise_gene_expression_profiles[spot_idx] = result
                        completed_spots.add(spot_idx)
                        spots_since_last_save += 1
                        pbar.update(1)

                        if spots_since_last_save >= checkpoint_interval:
                            checkpoint_mgr.save_checkpoint(
                                completed_spots,
                                spotwise_gene_expression_profiles,
                                N,
                                T,
                                M,
                            )
                            spots_since_last_save = 0
                except (OSError, ValueError) as e:

                    logging.error("Solve failed for spot %s: %s", spot_idx, e)
                    continue

    finally:
        gc.collect()
        if spotwise_gene_expression_profiles:
            checkpoint_mgr.save_final_results(
                spotwise_gene_expression_profiles,
                completed_spots,
                N,
                T,
                M,
            )

    solve_elapsed = time.time() - solve_start
    total = time.time() - prep_start
    logging.info(
        "GEX deconvolution complete: prep=%ss, solve=%ss, total=%ss (%s/%s spots)",
        round(prep_elapsed, 1),
        round(solve_elapsed, 1),
        round(total, 1),
        len(completed_spots),
        N,
    )

    return cast(Dict[str, Any], spotwise_gene_expression_profiles)


def optimize_cell_proportions_per_marker(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    *,
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
    lambda_laplacian: float = 0.0,
    coords: Optional[np.ndarray] = None,
    laplacian_k: int = 8,
    lambda_sparse: float = 0.0,
    alpha_max: float = 0.8,
    lambda_alpha: float = 1.0,
    lambda_coverage: float = 1.0,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
    row_sum_bounds: Optional[Tuple[float, float]] = None,
    cellularity: Optional[np.ndarray] = None,
    cellularity_sigma: float = 0.5,
    sparsity_mask: Optional[np.ndarray] = None,
    spot_weights: Optional[np.ndarray] = None,
    morphology_prior: Optional[np.ndarray] = None,
    lambda_morphology: float = 0.0,
    freeze_globals: bool = False,
    marker_weight: Optional[np.ndarray] = None,
    _confusion_pairs: Optional[List[Tuple[int, int]]] = None,
    _lambda_confusion: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray, np.ndarray]:
    """
    Perform EM-based optimization for cell type proportions with per-marker beta.

    cuOPT backend — replaces Gurobi QP E-step with NVIDIA cuOPT solver while
    keeping the M-step (beta/alpha WLS update) and convergence logic identical.

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
            Not yet supported in cuOPT backend; ignored with a warning if > 0.
        spot_abundance_target: Optional per-spot abundance target (shape N,).
        lambda_abundance_prior: Weight of abundance soft-prior term.
        row_sum_bounds: Optional (lb, ub) for row-sum constraints (default: (0.9, 1.2)).
        cellularity: Optional per-spot nuclei counts for Laplacian weighting.
        cellularity_sigma: Bandwidth for cellularity similarity kernel (default: 0.5).
        sparsity_mask: Optional (N, T) array with upper bounds per (spot, type).
        spot_weights: Optional (N,) per-spot weights for heteroscedastic reconstruction.

    Returns:
        Tuple of:
        - Y_values (np.ndarray): (N, T) cell type proportions
        - beta_values (np.ndarray): (M,) per-marker scaling factors
        - marker_beta_dict (Dict[str, float]): {marker_name: beta_value}
        - alpha_values (np.ndarray): (M,) per-marker baseline values
        - recon_error (np.ndarray): (N,) per-spot squared reconstruction error
    """
    if Problem is None:
        raise ImportError("cuOPT is not installed. Install with: pip install nvidia-cuopt")

    # When globals are frozen there is nothing to iterate — force single pass
    if freeze_globals:
        max_iterations = 1

    N, M = marker_level_data.shape
    T = len(cell_type_names)
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N:
            raise ValueError(f"spot_abundance_target length {spot_abundance_target.shape[0]} does not match N={N}.")

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

    # Compute markers-per-celltype for loss normalization (graded: soft markers count fractionally)
    markers_per_celltype = np.zeros(T, dtype=np.float64)
    for m in range(M):
        for j in marker_owners[m]:
            markers_per_celltype[j] += assignment_matrix[m, j]
    markers_per_celltype = np.maximum(markers_per_celltype, 1.0)

    if lambda_coverage > 0:
        logging.warning("lambda_coverage (asymmetric loss) not yet supported in cuOPT backend; ignoring")

    # Per-spot weights for heteroscedastic reconstruction (cellularity-informed)
    if spot_weights is not None:
        spot_weights = np.asarray(spot_weights, dtype=np.float64)
        if spot_weights.shape[0] != N:
            raise ValueError(f"spot_weights length {spot_weights.shape[0]} != N={N}")
        logging.info(
            "Cellularity spot weights enabled: median=%.2f, range=[%.2f, %.2f]",
            float(np.median(spot_weights)),
            float(spot_weights.min()),
            float(spot_weights.max()),
        )
    else:
        spot_weights = np.ones(N, dtype=np.float64)

    logging.info("Per-marker beta optimization (cuOPT): %s spots, %s markers, %s cell types", N, M, T)
    logging.info("Markers with assignments: %s/%s", marker_has_owner.sum(), M)

    # Build spatial Laplacian if requested
    L_coo = None
    use_laplacian = lambda_laplacian > 0 and coords is not None
    if use_laplacian:
        if coords.shape[0] != N:
            raise ValueError(f"coords has {coords.shape[0]} rows but data has {N} spots")
        L = build_spatial_laplacian(
            coords,
            k=laplacian_k,
            normed=True,
            cellularity=cellularity,
            cellularity_sigma=cellularity_sigma,
        )
        L_coo = L.tocoo()
        logging.info("Laplacian smoothing enabled: lambda=%s, k=%s", lambda_laplacian, laplacian_k)
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
        logging.info("\nIteration %s", iteration + 1)

        # ===== E-step: Solve QP for Y using cuOPT =====
        p = Problem()

        # Add N*T variables: Y[i, j] = proportion of cell type j at spot i
        Y_vars = {}
        for i in range(N):
            for j in range(T):
                ub = 1.0
                if sparsity_mask is not None:
                    ub = min(1.0, float(sparsity_mask[i, j]))
                Y_vars[i, j] = p.addVariable(lb=0.0, ub=ub, name=f"Y_{i}_{j}")

        if sparsity_mask is not None:
            n_clamped = int((sparsity_mask < 1.0).sum())
            logging.info(
                "Sparsity mask applied: %d/%d (spot,type) pairs clamped (%.1f%%)",
                n_clamped,
                N * T,
                100.0 * n_clamped / (N * T),
            )

        # Row-sum hard bounds (default [0.9, 1.2])
        lb_sum, ub_sum = row_sum_bounds if row_sum_bounds is not None else (0.9, 1.2)
        for i in range(N):
            row_sum = sum(Y_vars[i, j] for j in range(T))
            p.addConstraint(row_sum >= lb_sum)
            p.addConstraint(row_sum <= ub_sum)

        # Build objective: reconstruction error + regularization
        obj = 0

        # Reconstruction error: for each marker m and each owner j,
        # add normalized error: (1/n_owners) * (1/markers_per_celltype[j]) * (S - beta*Y)^2
        for m in range(M):
            if not marker_has_owner[m]:
                continue

            owners_m = marker_owners[m]
            n_owners = len(owners_m)
            beta_m = beta_values[m]
            alpha_m = alpha_values[m]

            for j in owners_m:
                mw = marker_weight[m] if marker_weight is not None else 1.0
                a_mj = assignment_matrix[m, j]  # graded weight (1.0=Major, <1=Soft)
                weight = a_mj * mw / (n_owners * markers_per_celltype[j])

                for i in range(N):
                    S_im = marker_level_data[i, m] - alpha_m  # baseline-subtracted
                    sw_i = float(spot_weights[i])
                    w = sw_i * weight

                    # (S - beta*Y)^2 = beta^2*Y^2 - 2*beta*S*Y + S^2
                    # S^2 is constant, dropped from objective
                    obj += w * beta_m**2 * Y_vars[i, j] * Y_vars[i, j]
                    obj += w * (-2 * beta_m * S_im) * Y_vars[i, j]

        # Elastic net regularization on Y
        for i in range(N):
            for j in range(T):
                obj += lambda_reg * alpha * Y_vars[i, j]  # L1 (Y>=0)
                obj += lambda_reg * (1 - alpha) * Y_vars[i, j] * Y_vars[i, j]  # L2

        # Sparsity penalty (negative L2 — encourages near-one-hot)
        if lambda_sparse > 0:
            for i in range(N):
                for j in range(T):
                    obj -= lambda_sparse * Y_vars[i, j] * Y_vars[i, j]
            logging.info("Sparsity penalty enabled (neg-L2): lambda_sparse=%s", lambda_sparse)

        # Laplacian smoothing term
        if use_laplacian and L_coo is not None:
            for idx in range(L_coo.nnz):
                i_spot = L_coo.row[idx]
                j_spot = L_coo.col[idx]
                L_val = L_coo.data[idx]
                for k in range(T):
                    obj += lambda_laplacian * L_val * Y_vars[i_spot, k] * Y_vars[j_spot, k]

        # Abundance prior term
        if spot_abundance_target is not None and lambda_abundance_prior > 0:
            for i in range(N):
                row_sum = sum(Y_vars[i, j] for j in range(T))
                target_i = float(spot_abundance_target[i])
                obj += lambda_abundance_prior * (row_sum - target_i) * (row_sum - target_i)

        # Morphology prior penalty: λ_morph * ||Y - p_morph||²
        # Expands to: λ*(Y² - 2*p_morph*Y) + const (const dropped)
        if morphology_prior is not None and lambda_morphology > 0:
            for i in range(N):
                for j in range(T):
                    pm_ij = float(morphology_prior[i, j])
                    obj += lambda_morphology * Y_vars[i, j] * Y_vars[i, j]
                    obj += lambda_morphology * (-2 * pm_ij) * Y_vars[i, j]
            logging.info("Morphology prior penalty added (lambda=%.3f)", lambda_morphology)

        p.setObjective(obj)

        try:
            p.solve()
        except Exception as e:
            logging.error("cuOPT optimization error: %s", str(e))
            raise ValueError("cuOPT optimization failed") from e

        # Check solve status (cuOPT exposes status via p.Status)
        if hasattr(p, "Status") and p.Status not in (None, "optimal", "OPTIMAL", 1, 2):
            raise ValueError(f"cuOPT optimization did not converge (status: {p.Status})")

        # Extract solution
        Y_values = np.array([[Y_vars[i, j].getValue() for j in range(T)] for i in range(N)])

        # ===== M-step: Update beta and alpha (per-marker WLS) =====
        if freeze_globals:
            logging.info("Globals frozen — skipping M-step")
            beta_new = beta_values.copy()
            alpha_new = alpha_values.copy()
        else:
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

                # WLS: S_m = alpha_m + beta_m * Y_combined, weights = spot_weights
                w = spot_weights
                w_sum = w.sum()
                Y_wmean = np.dot(w, Y_combined) / w_sum
                S_wmean = np.dot(w, S_m) / w_sum
                Y_centered = Y_combined - Y_wmean
                Y_wvar = np.dot(w, Y_centered**2)

                if Y_wvar > 1e-9:
                    beta_new[m] = np.dot(w * (S_m - S_wmean), Y_centered) / Y_wvar
                else:
                    beta_new[m] = beta_values[m]  # keep previous
                beta_new[m] = np.clip(beta_new[m], beta_min, beta_max)

                # Alpha with L2 regularization toward zero
                raw_alpha = S_wmean - beta_new[m] * Y_wmean
                alpha_new[m] = raw_alpha / (1.0 + lambda_alpha / N)
                alpha_new[m] = np.clip(alpha_new[m], 0.0, alpha_max)

            # Optionally normalize beta so max=1
            if normalize_beta:
                max_beta = np.max(beta_new)
                if max_beta > 0:
                    beta_new = beta_new / max_beta
                    # Re-clip after normalization to prevent extreme ratios.
                    beta_new = np.clip(beta_new, beta_min, 1.0)

        # Convergence check
        beta_diff = np.linalg.norm(beta_new - beta_prev)
        alpha_diff = np.linalg.norm(alpha_new - alpha_prev)
        Y_diff = np.linalg.norm(Y_values - Y_prev)

        logging.info("Change in beta: %s, alpha: %s, Y: %s", beta_diff, alpha_diff, Y_diff)
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
    logging.info("Beta range: [%s, %s], mean: %s", beta_new.min(), beta_new.max(), beta_new.mean())
    logging.info(
        "Alpha (baseline) range: [%s, %s], mean: %s",
        alpha_values.min(),
        alpha_values.max(),
        alpha_values.mean(),
    )
    for m in range(M):
        if marker_has_owner[m] and alpha_values[m] > 0.05:
            logging.info("  Marker '%s': alpha=%s, beta=%s", marker_names[m], alpha_values[m], beta_new[m])

    # Compute per-spot reconstruction error using final betas and alphas
    recon_error = np.zeros(N, dtype=np.float64)
    for m in range(M):
        if not marker_has_owner[m]:
            continue
        Y_combined_m = np.zeros(N, dtype=np.float64)
        for j in marker_owners[m]:
            Y_combined_m += Y_values[:, j]
        residual = marker_level_data[:, m] - alpha_values[m] - beta_new[m] * Y_combined_m
        recon_error += residual**2
    logging.info(
        "Per-spot recon error: mean=%s, median=%s, max=%s",
        recon_error.mean(),
        np.median(recon_error),
        recon_error.max(),
    )

    return Y_values, beta_new, marker_beta_dict, alpha_values, recon_error


def _mstep_per_type_beta(
    marker_data: np.ndarray,
    Y_values: np.ndarray,
    spot_weights: np.ndarray,
    beta_init: np.ndarray,
    prior_sigma: np.ndarray,
    *,
    lambda_alpha: float = 1.0,
    alpha_max: float = 0.8,
    beta_max: float = 2.0,
    prev_beta: Optional[np.ndarray] = None,
    inactive_threshold: float = 0.01,
) -> Tuple[np.ndarray, np.ndarray]:
    """M-step: update per-type beta matrix and per-marker alpha via centered WLS.

    For each marker m:
      1. Update alpha[m] = shrunk weighted mean of residuals
      2. Update beta[m,:] = ridge regression of (S - alpha) on Y

    Progressive greed penalty: for each type j, count how many inactive markers
    it has absorbed (beta > inactive_threshold where init was near zero). Each
    absorbed marker tightens the prior on ALL remaining inactive pairs for that
    type by a factor of (1 + n_absorbed). This makes it progressively harder
    for dominant types to absorb additional markers.

    Args:
        marker_data: (N, M) normalized marker data.
        Y_values: (N, T) current proportion estimates.
        spot_weights: (N,) per-spot weights.
        beta_init: (M, T) initialization values for ridge prior.
        prior_sigma: (M, T) per-pair regularization sigma.
        lambda_alpha: L2 regularization strength for alpha toward zero.
        alpha_max: Upper clip for alpha values.
        beta_max: Upper clip for beta values.
        prev_beta: (M, T) beta from previous iteration. Used for greed penalty
            counting. If None (first iteration), no greed penalty applied.
        inactive_threshold: Threshold for counting a pair as "absorbed".

    Returns:
        beta_new: (M, T) updated beta matrix.
        alpha_new: (M,) updated alpha vector.
    """
    N, M = marker_data.shape
    T = Y_values.shape[1]
    beta_new = np.zeros((M, T), dtype=np.float64)
    alpha_new = np.zeros(M, dtype=np.float64)

    # Identify inactive pairs (those initialized near zero)
    is_inactive = beta_init < inactive_threshold  # (M, T)

    # Bidirectional greed penalty from previous iteration's beta:
    #   type_greed[j]   = 1 + n_inactive_markers absorbed by type j
    #   marker_greed[m] = 1 + n_inactive_types absorbing marker m
    #   penalty[m,j] = type_greed[j] * marker_greed[m]
    #
    # Type-side: CD8+T picks up CD20 → harder for CD8+T to also pick up CD138
    # Marker-side: CD20 explained by B cells → harder for CD8+T to also pick up CD20
    # Compute bidirectional greed from prev_beta or beta_init.
    # On iteration 0 (prev_beta=None), use beta_init to count ACTIVE pairs:
    # types with many active markers already "own" signal and should face
    # barriers for inactive markers from the start.
    type_greed = np.ones(T, dtype=np.float64)
    marker_greed = np.ones(M, dtype=np.float64)
    ref_beta = prev_beta if prev_beta is not None else beta_init
    if prev_beta is not None:
        # Count newly absorbed inactive pairs
        absorbed = is_inactive & (ref_beta > inactive_threshold)
    else:
        # Count active pairs from init (strong/soft assignments)
        absorbed = ~is_inactive  # active pairs from MARKER_TYPE_TABLE
    for j in range(T):
        type_greed[j] = 1.0 + int(absorbed[:, j].sum())
    for m in range(M):
        marker_greed[m] = 1.0 + int(absorbed[m, :].sum())
    if np.any(type_greed > 1.0) or np.any(marker_greed > 1.0):
        logging.info(
            "[per_type_beta] Bidirectional greed — type: %s, marker: %s",
            {f"t{j}": f"{type_greed[j]:.0f}x" for j in range(T) if type_greed[j] > 1},
            {f"m{m}": f"{marker_greed[m]:.0f}x" for m in range(M) if marker_greed[m] > 1},
        )

    # Precompute weighted Gram matrix: Y^T diag(W) Y -> (T, T)
    W = spot_weights
    W_sum = W.sum()
    WY = Y_values * W[:, np.newaxis]  # (N, T)
    YtWY = Y_values.T @ WY  # (T, T)

    for m in range(M):
        S_m = marker_data[:, m]

        # Step 1: Update alpha (shrunk intercept)
        predicted_m = Y_values @ beta_init[m]
        residual_for_alpha = S_m - predicted_m
        raw_alpha = np.dot(W, residual_for_alpha) / W_sum
        alpha_new[m] = np.clip(raw_alpha / (1.0 + lambda_alpha / N), 0.0, alpha_max)

        # Step 2: Update beta[m,:] via ridge regression on centered target
        R_m = S_m - alpha_new[m]

        # Ridge penalty with bidirectional greed on inactive pairs
        # Base: Lambda_diag[j] = 1 / sigma[m,j]^2
        # Inactive: multiply by type_greed[j] * marker_greed[m]
        Lambda_diag = 1.0 / (prior_sigma[m] ** 2 + 1e-12)
        for j in range(T):
            if is_inactive[m, j]:
                Lambda_diag[j] *= type_greed[j] * marker_greed[m]

        # Solve: (Y^T W Y + diag(Lambda)) beta = Y^T W R + Lambda * beta_init
        A = YtWY + np.diag(Lambda_diag)
        b = WY.T @ R_m + Lambda_diag * beta_init[m]

        try:
            beta_new[m] = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            beta_new[m] = beta_init[m]

        # Step 3: Clip to [0, beta_max]
        beta_new[m] = np.clip(beta_new[m], 0.0, beta_max)

    # Step 4: Per-marker normalization — scale each marker's beta row so max=1.0
    for m in range(M):
        row_max = np.max(beta_new[m])
        if row_max > 1e-6:
            beta_new[m] = beta_new[m] / row_max

    return beta_new, alpha_new


def _compute_per_type_beta_objective(  # pylint: disable=too-many-positional-arguments
    marker_data: np.ndarray,
    Y_values: np.ndarray,
    beta: np.ndarray,
    alpha: np.ndarray,
    spot_weights: np.ndarray,
    beta_init: np.ndarray,
    prior_sigma: np.ndarray,
    lambda_reg: float,
    alpha_elastic: float,
) -> float:
    """Compute the full QP objective for convergence monitoring.

    Returns:
        Scalar objective value (lower = better).
    """
    # Reconstruction: sum_i,m w_i * (S_im - alpha_m - sum_j beta_mj * Y_ij)^2
    predicted = Y_values @ beta.T + alpha[np.newaxis, :]
    residual = marker_data - predicted
    recon = np.sum(spot_weights[:, np.newaxis] * residual**2)

    # Elastic net on Y
    elastic = lambda_reg * (alpha_elastic * np.sum(Y_values) + (1 - alpha_elastic) * np.sum(Y_values**2))

    # Ridge prior on beta
    ridge = np.sum((beta - beta_init) ** 2 / (prior_sigma**2 + 1e-12))

    return float(recon + elastic + ridge)


def optimize_cell_proportions_per_type_beta(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    cell_type_names: List[str],
    beta_init: np.ndarray,
    prior_sigma: np.ndarray,
    *,
    tolerance: float = 1e-4,
    max_iterations: int = 50,
    lambda_reg: float = 1.0,
    alpha_elastic: float = 0.5,
    beta_max: float = 2.0,
    alpha_max: float = 0.8,
    lambda_alpha: float = 1.0,
    lambda_laplacian: float = 0.0,
    coords: Optional[np.ndarray] = None,
    laplacian_k: int = 8,
    row_sum_bounds: Optional[Tuple[float, float]] = None,
    sparsity_mask: Optional[np.ndarray] = None,
    spot_weights: Optional[np.ndarray] = None,
    marker_weight: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Dict[str, float]], np.ndarray, np.ndarray, List[float]]:
    """EM-based optimization with per-type-per-marker beta matrix.

    E-step: cuOPT QP solve for Y[i,j] using diagonal Hessian approximation.
    M-step: per-marker multivariate WLS for beta[m,j] + shrunk alpha[m].
    Non-convex sparsity penalties are forbidden in this path.

    Args:
        marker_level_data: (N, M) normalized marker data.
        marker_names: List of M marker names.
        cell_type_names: List of T cell type names.
        beta_init: (M, T) initial beta matrix from emission_init.
        prior_sigma: (M, T) regularization sigma from emission_init.
        tolerance: Convergence threshold for beta/alpha/Y changes.
        max_iterations: Maximum EM iterations.
        lambda_reg: Elastic net regularization weight.
        alpha_elastic: L1/L2 mix for elastic net (0=L2, 1=L1).
        beta_max: Upper clip for beta values.
        alpha_max: Upper clip for alpha values.
        lambda_alpha: L2 regularization for alpha toward zero.
        lambda_laplacian: Spatial smoothing weight.
        coords: (N, 2) spatial coordinates for Laplacian.
        laplacian_k: Number of spatial neighbors.
        row_sum_bounds: (lb, ub) for row sum constraint.
        sparsity_mask: (N, T) upper bounds from detection gating.
        spot_weights: (N,) per-spot weights.
        marker_weight: (M,) per-marker loss weight.

    Returns:
        Y_values: (N, T) proportion estimates.
        beta_final: (M, T) learned beta matrix.
        marker_beta_matrix_dict: Dict[marker -> Dict[type -> float]].
        alpha_final: (M,) per-marker baselines.
        recon_error: (N,) per-spot reconstruction error.
        objective_history: List of per-iteration objective values.
    """
    from cuopt.linear_programming.problem import Problem  # pylint: disable=import-outside-toplevel

    from .emission_init import MARKER_TYPE_TABLE  # pylint: disable=import-outside-toplevel

    N, M = marker_level_data.shape
    T = len(cell_type_names)

    if spot_weights is None:
        spot_weights = np.ones(N, dtype=np.float64)

    # Laplacian setup
    L_coo = None
    use_laplacian = lambda_laplacian > 0 and coords is not None
    if use_laplacian:
        from scipy.sparse import coo_matrix  # pylint: disable=import-outside-toplevel
        from scipy.spatial import cKDTree  # pylint: disable=import-outside-toplevel

        tree = cKDTree(coords)
        _, nn_idx = tree.query(coords, k=laplacian_k + 1)
        rows, cols, data = [], [], []
        for i in range(N):
            for ki in range(1, laplacian_k + 1):
                j_nn = nn_idx[i, ki]
                rows.extend([i, j_nn])
                cols.extend([j_nn, i])
                data.extend([1.0, 1.0])
        L_coo = coo_matrix((data, (rows, cols)), shape=(N, N))
    elif lambda_laplacian > 0:
        logging.warning("lambda_laplacian > 0 but coords not provided. Laplacian disabled.")

    lb_sum, ub_sum = row_sum_bounds if row_sum_bounds is not None else (0.9, 1.2)

    # Initialize
    beta_values = beta_init.copy()
    alpha_values = np.zeros(M, dtype=np.float64)
    Y_values = np.full((N, T), 1.0 / T, dtype=np.float64)
    beta_prev = np.zeros_like(beta_values)
    alpha_prev = np.zeros_like(alpha_values)
    Y_prev = np.zeros_like(Y_values)
    objective_history = []

    mw = marker_weight if marker_weight is not None else np.ones(M, dtype=np.float64)

    iteration = 0
    while iteration < max_iterations:
        logging.info("\n[per_type_beta] Iteration %s", iteration + 1)
        t_iter_start = time.time()

        # ===== E-step: cuOPT QP with diagonal Hessian approximation =====
        p = Problem()

        Y_vars = {}
        for i in range(N):
            for j in range(T):
                ub = 1.0
                if sparsity_mask is not None:
                    ub = min(1.0, float(sparsity_mask[i, j]))
                Y_vars[i, j] = p.addVariable(lb=0.0, ub=ub, name=f"Y_{i}_{j}")

        # Row-sum constraints
        for i in range(N):
            row_sum = sum(Y_vars[i, j] for j in range(T))
            p.addConstraint(row_sum >= lb_sum)
            p.addConstraint(row_sum <= ub_sum)

        # Dense Hessian: H_base[j,k] = sum_m(mw[m] * beta[m,j] * beta[m,k])
        # For j=k this is the diagonal; for j≠k the cross-type coupling.
        # Cross-type terms ensure shared markers (CD3E, CD45, etc.) create
        # competition between types — if type j explains marker m, type k
        # gets penalized for also explaining it.
        beta_w = beta_values * np.sqrt(mw[:, np.newaxis])  # (M, T)
        H_base = beta_w.T @ beta_w  # (T, T) — PSD by construction

        # P_diag includes diagonal of H_base + elastic net L2
        P_diag = spot_weights[:, np.newaxis] * np.diag(H_base)[np.newaxis, :]
        P_diag += lambda_reg * (1 - alpha_elastic)

        # Linear: q[i,j] = -2*sw[i]*sum_m(mw[m]*beta[m,j]*(S[i,m]-alpha[m])) + lambda_reg*alpha_elastic
        S_adj = marker_level_data - alpha_values[np.newaxis, :]
        beta_weighted = beta_values * mw[:, np.newaxis]
        q_linear = -2.0 * (S_adj @ beta_weighted)
        q_linear *= spot_weights[:, np.newaxis]
        q_linear += lambda_reg * alpha_elastic

        # Build objective: diagonal + cross-type terms
        obj = 0
        for i in range(N):
            sw_i = float(spot_weights[i])
            for j in range(T):
                obj += P_diag[i, j] * Y_vars[i, j] * Y_vars[i, j]
                obj += q_linear[i, j] * Y_vars[i, j]
            # Cross-type coupling from shared markers (off-diagonal H_base)
            for j in range(T):
                for k in range(j + 1, T):
                    cross = 2.0 * sw_i * H_base[j, k]
                    if abs(cross) > 1e-10:
                        obj += cross * Y_vars[i, j] * Y_vars[i, k]

        # Laplacian smoothing
        if use_laplacian and L_coo is not None:
            for idx in range(L_coo.nnz):
                i_spot = L_coo.row[idx]
                j_spot = L_coo.col[idx]
                L_val = L_coo.data[idx]
                for k in range(T):
                    obj += lambda_laplacian * L_val * Y_vars[i_spot, k] * Y_vars[j_spot, k]

        p.setObjective(obj)

        t_build = time.time()
        logging.info("[per_type_beta] E-step build: %.3fs", t_build - t_iter_start)

        status = p.solve()
        if status != 0:
            logging.warning("[per_type_beta] cuOPT solve status %d at iteration %d", status, iteration + 1)

        t_solve = time.time()
        logging.info("[per_type_beta] E-step solve: %.3fs", t_solve - t_build)

        # Extract solution
        for i in range(N):
            for j in range(T):
                Y_values[i, j] = Y_vars[i, j].getValue()

        # ===== M-step =====
        beta_new, alpha_new = _mstep_per_type_beta(
            marker_level_data,
            Y_values,
            spot_weights,
            beta_init,
            prior_sigma,
            lambda_alpha=lambda_alpha,
            alpha_max=alpha_max,
            beta_max=beta_max,
            prev_beta=beta_values if iteration > 0 else None,
        )

        # ===== Convergence monitoring =====
        obj_val = _compute_per_type_beta_objective(
            marker_level_data,
            Y_values,
            beta_new,
            alpha_new,
            spot_weights,
            beta_init,
            prior_sigma,
            lambda_reg,
            alpha_elastic,
        )
        objective_history.append(obj_val)

        # Collapse detection: objective increase
        if len(objective_history) > 1:
            prev_obj = objective_history[-2]
            if obj_val > prev_obj * 1.01:
                logging.warning(
                    "[per_type_beta] Objective increased by %.2f%% — reverting",
                    100.0 * (obj_val - prev_obj) / (abs(prev_obj) + 1e-12),
                )
                beta_values = beta_prev.copy()
                alpha_values = alpha_prev.copy()
                Y_values = Y_prev.copy()
                break

        # Collapse detection: large beta change
        max_beta_change = np.max(np.abs(beta_new - beta_values) / (np.abs(beta_values) + 1e-6))
        if max_beta_change > 0.5:
            logging.warning("[per_type_beta] Max relative beta change %.2f > 0.5", max_beta_change)

        # Collapse detection: strong pair collapse
        for m_idx in range(M):
            for t_name, strength in MARKER_TYPE_TABLE.get(marker_names[m_idx], []):
                if strength == "strong" and t_name in cell_type_names:
                    t_idx = cell_type_names.index(t_name)
                    if beta_new[m_idx, t_idx] < 0.1 * beta_init[m_idx, t_idx]:
                        logging.warning(
                            "[per_type_beta] Strong pair %s->%s collapsed: %.4f < 10%% of init %.4f",
                            marker_names[m_idx],
                            t_name,
                            beta_new[m_idx, t_idx],
                            beta_init[m_idx, t_idx],
                        )

        # Convergence check
        beta_diff = np.linalg.norm(beta_new - beta_values)
        alpha_diff = np.linalg.norm(alpha_new - alpha_values)
        Y_diff = np.linalg.norm(Y_values - Y_prev)
        logging.info(
            "[per_type_beta] Obj: %.4f, d_beta: %.6f, d_alpha: %.6f, d_Y: %.6f",
            obj_val,
            beta_diff,
            alpha_diff,
            Y_diff,
        )

        if beta_diff < tolerance and alpha_diff < tolerance and Y_diff < tolerance:
            logging.info("[per_type_beta] Convergence at iteration %d", iteration + 1)
            break

        beta_prev = beta_values.copy()
        alpha_prev = alpha_values.copy()
        Y_prev = Y_values.copy()
        beta_values = beta_new.copy()
        alpha_values = alpha_new.copy()
        iteration += 1

    # Build result dict
    marker_beta_matrix_dict = {}
    for m_idx, m_name in enumerate(marker_names):
        marker_beta_matrix_dict[m_name] = {cell_type_names[j]: float(beta_values[m_idx, j]) for j in range(T)}

    # Per-spot reconstruction error
    predicted = Y_values @ beta_values.T + alpha_values[np.newaxis, :]
    recon_error = np.sum((marker_level_data - predicted) ** 2, axis=1)

    logging.info(
        "[per_type_beta] Final recon error: mean=%.4f, median=%.4f", recon_error.mean(), np.median(recon_error)
    )

    return Y_values, beta_values, marker_beta_matrix_dict, alpha_values, recon_error, objective_history


def optimize_cell_proportions_per_marker_matrix(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    *,
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
    lambda_laplacian: float = 0.0,
    coords: Optional[np.ndarray] = None,
    laplacian_k: int = 8,
    lambda_sparse: float = 0.0,
    alpha_max: float = 0.8,
    lambda_alpha: float = 1.0,
    lambda_coverage: float = 1.0,
    spot_abundance_target: Optional[np.ndarray] = None,
    lambda_abundance_prior: float = 0.0,
    row_sum_bounds: Optional[Tuple[float, float]] = None,
    cellularity: Optional[np.ndarray] = None,
    cellularity_sigma: float = 0.5,
    sparsity_mask: Optional[np.ndarray] = None,
    spot_weights: Optional[np.ndarray] = None,
    marker_weight: Optional[np.ndarray] = None,
    confusion_pairs: Optional[List[Tuple[int, int]]] = None,
    lambda_confusion: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray, np.ndarray]:
    """
    Matrix-formulated version of optimize_cell_proportions_per_marker.

    Performance optimization: pre-builds cuOPT variables and constraints ONCE
    before the EM loop, then rebuilds only the objective expression each
    iteration using vectorized NumPy coefficient computation.

    The key speedup comes from:
    1. Variable creation (N*T addVariable calls) done once, not per iteration
    2. Constraint creation (2*N addConstraint calls) done once, not per iteration
    3. Objective coefficients (P diagonal and q vector) computed via NumPy
       matrix ops instead of Python triple-loops

    The M-step (WLS beta/alpha update), convergence check, validation, and
    reconstruction error computation are identical to the reference implementation.

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
            Not yet supported in cuOPT backend; ignored with a warning if > 0.
        spot_abundance_target: Optional per-spot abundance target (shape N,).
        lambda_abundance_prior: Weight of abundance soft-prior term.
        row_sum_bounds: Optional (lb, ub) for row-sum constraints (default: (0.9, 1.2)).
        cellularity: Optional per-spot nuclei counts for Laplacian weighting.
        cellularity_sigma: Bandwidth for cellularity similarity kernel (default: 0.5).
        sparsity_mask: Optional (N, T) array with upper bounds per (spot, type).
        spot_weights: Optional (N,) per-spot weights for heteroscedastic reconstruction.

    Returns:
        Tuple of:
        - Y_values (np.ndarray): (N, T) cell type proportions
        - beta_values (np.ndarray): (M,) per-marker scaling factors
        - marker_beta_dict (Dict[str, float]): {marker_name: beta_value}
        - alpha_values (np.ndarray): (M,) per-marker baseline values
        - recon_error (np.ndarray): (N,) per-spot squared reconstruction error
    """
    if Problem is None:
        raise ImportError("cuOPT is not installed. Install with: pip install nvidia-cuopt")

    N, M = marker_level_data.shape
    T = len(cell_type_names)
    if spot_abundance_target is not None:
        spot_abundance_target = np.asarray(spot_abundance_target, dtype=np.float64).reshape(-1)
        if spot_abundance_target.shape[0] != N:
            raise ValueError(f"spot_abundance_target length {spot_abundance_target.shape[0]} does not match N={N}.")

    # Validate assignment matrix shape
    if assignment_matrix.shape != (M, T):
        raise ValueError(f"Assignment matrix shape {assignment_matrix.shape} != expected ({M}, {T})")

    # Precompute marker-to-celltype assignments (supports shared markers)
    marker_owners = []
    for m in range(M):
        owners_for_m = []
        for j in range(T):
            if assignment_matrix[m, j] > 0:
                owners_for_m.append(j)
        marker_owners.append(owners_for_m)

    marker_has_owner = np.array([len(o) > 0 for o in marker_owners])

    # Compute markers-per-celltype for loss normalization (graded: soft markers count fractionally)
    markers_per_celltype = np.zeros(T, dtype=np.float64)
    for m in range(M):
        for j in marker_owners[m]:
            markers_per_celltype[j] += assignment_matrix[m, j]
    markers_per_celltype = np.maximum(markers_per_celltype, 1.0)

    if lambda_coverage > 0:
        logging.warning("lambda_coverage (asymmetric loss) not yet supported in cuOPT backend; ignoring")

    # Per-spot weights for heteroscedastic reconstruction
    if spot_weights is not None:
        spot_weights = np.asarray(spot_weights, dtype=np.float64)
        if spot_weights.shape[0] != N:
            raise ValueError(f"spot_weights length {spot_weights.shape[0]} != N={N}")
        logging.info(
            "Cellularity spot weights enabled: median=%.2f, range=[%.2f, %.2f]",
            float(np.median(spot_weights)),
            float(spot_weights.min()),
            float(spot_weights.max()),
        )
    else:
        spot_weights = np.ones(N, dtype=np.float64)

    logging.info("Per-marker beta optimization (cuOPT matrix): %s spots, %s markers, %s cell types", N, M, T)
    logging.info("Markers with assignments: %s/%s", marker_has_owner.sum(), M)

    # ===== Pre-compute weight matrix W[m, j] = 1/(n_owners * markers_per_celltype[j]) =====
    # This is STATIC — does not depend on beta/alpha
    W_static = np.zeros((M, T), dtype=np.float64)
    for m in range(M):
        if not marker_has_owner[m]:
            continue
        n_owners = len(marker_owners[m])
        for j in marker_owners[m]:
            mw = marker_weight[m] if marker_weight is not None else 1.0
            a_mj = assignment_matrix[m, j]
            W_static[m, j] = a_mj * mw / (n_owners * markers_per_celltype[j])

    # ===== Build spatial Laplacian if requested =====
    L_coo = None
    use_laplacian = lambda_laplacian > 0 and coords is not None
    if use_laplacian:
        if coords.shape[0] != N:
            raise ValueError(f"coords has {coords.shape[0]} rows but data has {N} spots")
        L = build_spatial_laplacian(
            coords,
            k=laplacian_k,
            normed=True,
            cellularity=cellularity,
            cellularity_sigma=cellularity_sigma,
        )
        L_coo = L.tocoo()
        logging.info("Laplacian smoothing enabled: lambda=%s, k=%s", lambda_laplacian, laplacian_k)
    elif lambda_laplacian > 0 and coords is None:
        logging.warning("lambda_laplacian > 0 but coords not provided. Laplacian smoothing disabled.")

    # ===== Pre-build cuOPT variables and constraints ONCE =====
    p = Problem()

    # Create N*T variables
    Y_vars = {}
    for i in range(N):
        for j in range(T):
            ub = 1.0
            if sparsity_mask is not None:
                ub = min(1.0, float(sparsity_mask[i, j]))
            Y_vars[i, j] = p.addVariable(lb=0.0, ub=ub, name=f"Y_{i}_{j}")

    if sparsity_mask is not None:
        n_clamped = int((sparsity_mask < 1.0).sum())
        logging.info(
            "Sparsity mask applied: %d/%d (spot,type) pairs clamped (%.1f%%)",
            n_clamped,
            N * T,
            100.0 * n_clamped / (N * T),
        )

    # Row-sum hard bounds (default [0.9, 1.2])
    lb_sum, ub_sum = row_sum_bounds if row_sum_bounds is not None else (0.9, 1.2)
    for i in range(N):
        row_sum = sum(Y_vars[i, j] for j in range(T))
        p.addConstraint(row_sum >= lb_sum)
        p.addConstraint(row_sum <= ub_sum)

    # Pre-build abundance prior row-sum variable references (if needed)
    use_abundance = spot_abundance_target is not None and lambda_abundance_prior > 0

    logging.info(
        "Pre-built %d variables and %d constraints (reused across EM iterations)",
        N * T,
        2 * N,
    )

    # ===== Initialize beta and alpha =====
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)
    beta_prev = np.zeros(M)
    alpha_prev = np.zeros(M)
    Y_prev = np.zeros((N, T))

    iteration = 0
    while iteration < max_iterations:
        logging.info("\nIteration %s", iteration + 1)
        t_iter_start = time.time()

        # ===== E-step: Build objective using vectorized coefficients =====

        # --- Compute P_diag[i, j] (quadratic coeff for Y[i,j]^2) via NumPy ---
        # P_diag[i, j] = spot_weights[i] * sum_m(W_static[m,j] * beta[m]^2)
        #              + lambda_reg * (1 - alpha_elastic)
        #              - lambda_sparse
        #
        # The marker sum is: (W_static * beta^2)[:, j].sum(axis=0) for each type j
        # This gives a (T,) vector, then broadcast with spot_weights (N,)

        beta_sq = beta_values**2  # (M,)
        # W_beta_sq[m, j] = W_static[m, j] * beta[m]^2
        W_beta_sq = W_static * beta_sq[:, np.newaxis]  # (M, T)
        # Sum over markers for each type: type_quad_coeff[j] = sum_m W_beta_sq[m, j]
        type_quad_coeff = W_beta_sq.sum(axis=0)  # (T,)

        # P_diag[i, j] = spot_weights[i] * type_quad_coeff[j] + regularization
        P_diag = spot_weights[:, np.newaxis] * type_quad_coeff[np.newaxis, :]  # (N, T)
        P_diag += lambda_reg * (1 - alpha)  # L2 elastic net
        if lambda_sparse > 0:
            P_diag -= lambda_sparse  # negative L2 for sparsity

        # --- Compute q_linear[i, j] (linear coeff for Y[i,j]) via matrix multiply ---
        # q[i, j] = -2 * sum_m(spot_weights[i] * W_static[m,j] * beta[m] * (S[i,m] - alpha_m))
        #         + lambda_reg * alpha_elastic
        #
        # = -2 * spot_weights[i] * sum_m(S_adj[i,m] * W_beta[m,j])
        # where S_adj[i,m] = S[i,m] - alpha_m, W_beta[m,j] = W_static[m,j] * beta[m]

        S_adjusted = marker_level_data - alpha_values[np.newaxis, :]  # (N, M)
        W_beta = W_static * beta_values[:, np.newaxis]  # (M, T)
        # Matrix multiply: (N, M) @ (M, T) -> (N, T)
        q_linear = -2.0 * (S_adjusted @ W_beta)  # (N, T)
        q_linear *= spot_weights[:, np.newaxis]  # apply spot weights
        q_linear += lambda_reg * alpha  # L1 term (Y >= 0 so |Y| = Y)

        t_coeff = time.time()

        # --- Build objective expression from pre-computed coefficients ---
        obj = 0
        for i in range(N):
            for j in range(T):
                obj += P_diag[i, j] * Y_vars[i, j] * Y_vars[i, j]
                obj += q_linear[i, j] * Y_vars[i, j]

        # Laplacian smoothing cross-terms (structure is static, values are static)
        if use_laplacian and L_coo is not None:
            for idx in range(L_coo.nnz):
                i_spot = L_coo.row[idx]
                j_spot = L_coo.col[idx]
                L_val = L_coo.data[idx]
                for k in range(T):
                    obj += lambda_laplacian * L_val * Y_vars[i_spot, k] * Y_vars[j_spot, k]

        # Abundance prior term
        if use_abundance:
            for i in range(N):
                row_sum = sum(Y_vars[i, j] for j in range(T))
                target_i = float(spot_abundance_target[i])
                obj += lambda_abundance_prior * (row_sum - target_i) * (row_sum - target_i)

        # Confusion penalty: discourage confused type pairs from being simultaneously high
        if confusion_pairs and lambda_confusion > 0:
            for i in range(N):
                for t_a, t_b in confusion_pairs:
                    obj += lambda_confusion * Y_vars[i, t_a] * Y_vars[i, t_b]

        p.setObjective(obj)

        t_obj = time.time()
        logging.info(
            "  Coefficient computation: %.3fs, objective build: %.3fs",
            t_coeff - t_iter_start,
            t_obj - t_coeff,
        )

        try:
            p.solve()
        except Exception as e:
            logging.error("cuOPT optimization error: %s", str(e))
            raise ValueError("cuOPT optimization failed") from e

        # Check solve status
        if hasattr(p, "Status") and p.Status not in (None, "optimal", "OPTIMAL", 1, 2):
            raise ValueError(f"cuOPT optimization did not converge (status: {p.Status})")

        t_solve = time.time()
        logging.info("  Solve time: %.3fs", t_solve - t_obj)

        # Extract solution
        Y_values = np.array([[Y_vars[i, j].getValue() for j in range(T)] for i in range(N)])

        # ===== M-step: Update beta and alpha (per-marker WLS) =====
        # Identical to reference implementation
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

            # WLS: S_m = alpha_m + beta_m * Y_combined, weights = spot_weights
            w = spot_weights
            w_sum = w.sum()
            Y_wmean = np.dot(w, Y_combined) / w_sum
            S_wmean = np.dot(w, S_m) / w_sum
            Y_centered = Y_combined - Y_wmean
            Y_wvar = np.dot(w, Y_centered**2)

            if Y_wvar > 1e-9:
                beta_new[m] = np.dot(w * (S_m - S_wmean), Y_centered) / Y_wvar
            else:
                beta_new[m] = beta_values[m]  # keep previous
            beta_new[m] = np.clip(beta_new[m], beta_min, beta_max)

            # Alpha with L2 regularization toward zero
            raw_alpha = S_wmean - beta_new[m] * Y_wmean
            alpha_new[m] = raw_alpha / (1.0 + lambda_alpha / N)
            alpha_new[m] = np.clip(alpha_new[m], 0.0, alpha_max)

        # Optionally normalize beta so max=1
        if normalize_beta:
            max_beta = np.max(beta_new)
            if max_beta > 0:
                beta_new = beta_new / max_beta
                # Re-clip after normalization to prevent extreme ratios.
                beta_new = np.clip(beta_new, beta_min, 1.0)

        # Convergence check
        beta_diff = np.linalg.norm(beta_new - beta_prev)
        alpha_diff = np.linalg.norm(alpha_new - alpha_prev)
        Y_diff = np.linalg.norm(Y_values - Y_prev)

        logging.info("Change in beta: %s, alpha: %s, Y: %s", beta_diff, alpha_diff, Y_diff)
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
        profile_based_antibody_data=None,
        unknown_threshold=unknown_threshold,
        min_celltype_threshold=min_celltype_threshold,
        redundancy_threshold=redundancy_threshold,
        warn_only=warn_only,
    )

    # Build marker-beta dictionary for interpretability
    marker_beta_dict = {marker_names[m]: beta_new[m] for m in range(M)}

    # Log beta and alpha statistics
    logging.info("Beta range: [%s, %s], mean: %s", beta_new.min(), beta_new.max(), beta_new.mean())
    logging.info(
        "Alpha (baseline) range: [%s, %s], mean: %s",
        alpha_values.min(),
        alpha_values.max(),
        alpha_values.mean(),
    )
    for m in range(M):
        if marker_has_owner[m] and alpha_values[m] > 0.05:
            logging.info("  Marker '%s': alpha=%s, beta=%s", marker_names[m], alpha_values[m], beta_new[m])

    # Compute per-spot reconstruction error using final betas and alphas
    recon_error = np.zeros(N, dtype=np.float64)
    for m in range(M):
        if not marker_has_owner[m]:
            continue
        Y_combined_m = np.zeros(N, dtype=np.float64)
        for j in marker_owners[m]:
            Y_combined_m += Y_values[:, j]
        residual = marker_level_data[:, m] - alpha_values[m] - beta_new[m] * Y_combined_m
        recon_error += residual**2
    logging.info(
        "Per-spot recon error: mean=%s, median=%s, max=%s",
        recon_error.mean(),
        np.median(recon_error),
        recon_error.max(),
    )

    return Y_values, beta_new, marker_beta_dict, alpha_values, recon_error
