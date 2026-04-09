"""
Joint optimization for profile discovery and proportion estimation.

This module implements a Gurobi-based optimization that simultaneously learns:
- W[k,m]: Soft profile weights (how much marker m defines cell type k)
- Y[s,k]: Cell type proportions at each spot
- beta[m]: Per-marker scaling factors

This replaces the sequential approach of:
1. discover_profiles() -> profile selection
2. optimize_cell_proportions() -> proportion estimation

With a single joint optimization using Block Coordinate Descent (alternating minimization).
"""

import logging
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import gurobipy as gp
import numpy as np
import scipy.sparse as sp
from gurobipy import GRB
from sklearn.decomposition import NMF

from .gurobi_impl import build_spatial_laplacian


@dataclass
class JointOptimizationResult:
    """Container for joint optimization results."""

    W: np.ndarray  # (K, n_markers) profile weights
    Y: np.ndarray  # (n_spots, K) proportions
    beta: np.ndarray  # (n_markers,) scaling factors
    profiles: Dict[str, Dict[str, List[str]]]  # CITEgeist-compatible format
    K: int  # Selected number of cell types
    bic: float  # BIC score for selected K
    marker_names: List[str]  # Marker names in order
    index_to_name: Dict[int, str] = field(default_factory=dict)  # Mapping k -> profile name
    metadata: Dict[str, Any] = field(default_factory=dict)  # Diagnostics


def _initialize_nmf(
    X: np.ndarray,
    K: int,
    seed: int = 1234,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Initialize W, Y, beta using Non-negative Matrix Factorization.

    Args:
        X: Antibody data (n_spots, n_markers), should be non-negative.
        K: Number of cell types/components.
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (W, Y, beta) initial values.
    """
    n_spots, n_markers = X.shape

    # Ensure non-negative (CLR can produce negatives)
    X_nn = np.maximum(X, 0)

    # Add small constant to avoid all-zero columns/rows
    X_nn = X_nn + 1e-6

    try:
        nmf = NMF(
            n_components=K,
            init="nndsvd",
            random_state=seed,
            max_iter=200,
        )
        Y_init = nmf.fit_transform(X_nn)  # (n_spots, K)
        W_init = nmf.components_  # (K, n_markers)
    except Exception as e:
        logging.warning(f"NMF initialization failed: {e}. Using random initialization.")
        rng = np.random.default_rng(seed)
        Y_init = rng.uniform(0.1, 0.9, size=(n_spots, K))
        W_init = rng.uniform(0.1, 0.9, size=(K, n_markers))

    # Normalize Y rows to sum to ~1
    Y_row_sums = Y_init.sum(axis=1, keepdims=True)
    Y_row_sums = np.maximum(Y_row_sums, 1e-6)
    Y_init = Y_init / Y_row_sums

    # Normalize W to [0, 1] range per row
    W_row_max = W_init.max(axis=1, keepdims=True)
    W_row_max = np.maximum(W_row_max, 1e-6)
    W_init = W_init / W_row_max

    # Initialize beta = 1
    beta_init = np.ones(n_markers)

    return W_init, Y_init, beta_init


def _initialize_random(
    n_spots: int,
    n_markers: int,
    K: int,
    seed: int = 1234,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Random initialization for W, Y, beta."""
    rng = np.random.default_rng(seed)

    # Random Y, normalized to sum to 1
    Y_init = rng.uniform(0.1, 0.9, size=(n_spots, K))
    Y_init = Y_init / Y_init.sum(axis=1, keepdims=True)

    # Random W in [0, 1]
    W_init = rng.uniform(0.0, 1.0, size=(K, n_markers))

    # Beta = 1
    beta_init = np.ones(n_markers)

    return W_init, Y_init, beta_init


def _optimize_Y(
    X: np.ndarray,
    W: np.ndarray,
    beta: np.ndarray,
    L_coo: sp.coo_matrix,
    lambda_spatial: float,
    Y_prev: Optional[np.ndarray] = None,
    max_y_change: float = 0.4,
    time_limit: float = 120.0,
) -> np.ndarray:
    """
    Optimize Y given fixed W and beta (convex QP).

    Objective:
        min ||X - Y @ diag(beta) @ W^T||^2_F + lambda_spatial * sum_k Y[:,k]^T L Y[:,k]

    Constraints:
        Y[s,k] in [0, 1]
        0.9 <= sum_k Y[s,k] <= 1.2

    Args:
        X: Antibody data (n_spots, n_markers).
        W: Profile weights (K, n_markers), fixed.
        beta: Marker scaling (n_markers,), fixed.
        L_coo: Spatial Laplacian in COO format.
        lambda_spatial: Weight for Laplacian smoothing.
        Y_prev: Previous Y values for bounded updates.
        max_y_change: Maximum change allowed from Y_prev.
        time_limit: Gurobi time limit in seconds.

    Returns:
        Optimized Y array (n_spots, K).
    """
    n_spots, n_markers = X.shape
    K = W.shape[0]

    model = gp.Model("Y_subproblem")
    model.setParam("OutputFlag", 1)  # Enable verbose output for debugging
    model.setParam("LogToConsole", 1)
    model.setParam("TimeLimit", time_limit)
    logging.info(f"Y_subproblem: {n_spots} spots, {K} cell types, time_limit={time_limit}s")

    # Decision variables: Y[s, k]
    Y_vars = {}
    for s in range(n_spots):
        for k in range(K):
            lb = 0.0
            ub = 1.0
            if Y_prev is not None:
                lb = max(0.0, Y_prev[s, k] - max_y_change)
                ub = min(1.0, Y_prev[s, k] + max_y_change)
            Y_vars[s, k] = model.addVar(lb=lb, ub=ub, vtype=GRB.CONTINUOUS, name=f"Y_{s}_{k}")

    model.update()

    # Precompute coefficients: C[s, k] = sum_m beta[m] * W[k, m]
    # X_pred[s, m] = sum_k Y[s, k] * beta[m] * W[k, m]
    # Error = sum_{s,m} (X[s,m] - X_pred[s,m])^2

    # Expand: (X[s,m] - sum_k Y[s,k] * beta[m] * W[k,m])^2
    # = X[s,m]^2 - 2*X[s,m]*sum_k(...) + (sum_k(...))^2

    # Build objective
    obj = gp.QuadExpr()

    # Reconstruction error term
    for s in range(n_spots):
        for m in range(n_markers):
            x_val = X[s, m]
            # Linear coefficients for Y[s, k]
            coeff = np.array([beta[m] * W[k, m] for k in range(K)])

            # Add linear term: -2 * X[s,m] * sum_k coeff[k] * Y[s,k]
            for k in range(K):
                obj.addTerms(-2 * x_val * coeff[k], Y_vars[s, k])

            # Add quadratic term: (sum_k coeff[k] * Y[s,k])^2
            for k1 in range(K):
                for k2 in range(K):
                    obj.addTerms(coeff[k1] * coeff[k2], Y_vars[s, k1], Y_vars[s, k2])

            # Constant term X[s,m]^2 doesn't affect optimization

    # Laplacian smoothing term: lambda_spatial * sum_k Y[:,k]^T L Y[:,k]
    if lambda_spatial > 0 and L_coo.nnz > 0:
        for idx in range(L_coo.nnz):
            i = L_coo.row[idx]
            j = L_coo.col[idx]
            L_val = L_coo.data[idx]
            for k in range(K):
                obj.addTerms(lambda_spatial * L_val, Y_vars[i, k], Y_vars[j, k])

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraints: sum-to-one with slack
    for s in range(n_spots):
        model.addConstr(
            gp.quicksum(Y_vars[s, k] for k in range(K)) >= 0.9,
            name=f"sum_lb_{s}",
        )
        model.addConstr(
            gp.quicksum(Y_vars[s, k] for k in range(K)) <= 1.2,
            name=f"sum_ub_{s}",
        )

    model.optimize()

    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
        Y_result = np.array([[Y_vars[s, k].X for k in range(K)] for s in range(n_spots)])
        return Y_result
    else:
        logging.warning(f"Y optimization status: {model.status}. Returning previous Y.")
        return Y_prev if Y_prev is not None else np.ones((n_spots, K)) / K


def _optimize_W(
    X: np.ndarray,
    Y: np.ndarray,
    beta: np.ndarray,
    lambda_sparsity: float,
    lambda_distinct: float,
    max_markers_per_type: int,
    time_limit: float = 120.0,
) -> np.ndarray:
    """
    Optimize W given fixed Y and beta (convex QP).

    Objective:
        min ||X - Y @ diag(beta) @ W^T||^2_F
            + lambda_sparsity * ||W||_1
            + lambda_distinct * sum_{k1<k2} W[k1] . W[k2]

    Constraints:
        W[k,m] in [0, 1]
        sum_m W[k,m] >= 0.5 (each cell type has at least some markers)
        sum_m W[k,m] <= max_markers_per_type

    Args:
        X: Antibody data (n_spots, n_markers).
        Y: Proportions (n_spots, K), fixed.
        beta: Marker scaling (n_markers,), fixed.
        lambda_sparsity: Weight for L1 sparsity on W.
        lambda_distinct: Weight for profile distinctness.
        max_markers_per_type: Maximum markers per cell type.
        time_limit: Gurobi time limit in seconds.

    Returns:
        Optimized W array (K, n_markers).
    """
    n_spots, n_markers = X.shape
    K = Y.shape[1]

    model = gp.Model("W_subproblem")
    model.setParam("OutputFlag", 1)  # Enable verbose output for debugging
    model.setParam("LogToConsole", 1)
    model.setParam("TimeLimit", time_limit)
    logging.info(f"W_subproblem: {K} cell types, {n_markers} markers, time_limit={time_limit}s")

    # Decision variables: W[k, m]
    W_vars = {}
    for k in range(K):
        for m in range(n_markers):
            W_vars[k, m] = model.addVar(lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name=f"W_{k}_{m}")

    model.update()

    # Build objective
    obj = gp.QuadExpr()

    # Reconstruction error: ||X - Y @ diag(beta) @ W^T||^2_F
    # For each (s, m): (X[s,m] - sum_k Y[s,k] * beta[m] * W[k,m])^2
    for m in range(n_markers):
        for s in range(n_spots):
            x_val = X[s, m]
            # Coefficient for W[k, m] is Y[s, k] * beta[m]
            coeff = np.array([Y[s, k] * beta[m] for k in range(K)])

            # Linear term: -2 * X[s,m] * sum_k coeff[k] * W[k,m]
            for k in range(K):
                obj.addTerms(-2 * x_val * coeff[k], W_vars[k, m])

            # Quadratic term: (sum_k coeff[k] * W[k,m])^2
            for k1 in range(K):
                for k2 in range(K):
                    obj.addTerms(coeff[k1] * coeff[k2], W_vars[k1, m], W_vars[k2, m])

    # Sparsity term: lambda_sparsity * sum_{k,m} W[k,m]
    # Since W >= 0, L1 norm is just the sum
    if lambda_sparsity > 0:
        for k in range(K):
            for m in range(n_markers):
                obj.addTerms(lambda_sparsity, W_vars[k, m])

    # Distinctness term: Penalize markers used by multiple cell types
    # Original bilinear form: sum_{k1<k2, m} W[k1,m] * W[k2,m] is NON-CONVEX
    # Convex alternative: penalize (sum_k W[k,m])^2 per marker
    # This encourages each marker to be used by fewer cell types
    if lambda_distinct > 0:
        for m in range(n_markers):
            # (sum_k W[k,m])^2 = sum_{k1,k2} W[k1,m] * W[k2,m]
            # This is convex because the Hessian is positive semidefinite
            for k1 in range(K):
                for k2 in range(K):
                    obj.addTerms(lambda_distinct, W_vars[k1, m], W_vars[k2, m])

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraints: marker count bounds per cell type
    for k in range(K):
        model.addConstr(
            gp.quicksum(W_vars[k, m] for m in range(n_markers)) >= 0.5,
            name=f"W_lb_{k}",
        )
        model.addConstr(
            gp.quicksum(W_vars[k, m] for m in range(n_markers)) <= max_markers_per_type,
            name=f"W_ub_{k}",
        )

    model.optimize()

    if model.status == GRB.OPTIMAL or model.status == GRB.TIME_LIMIT:
        W_result = np.array([[W_vars[k, m].X for m in range(n_markers)] for k in range(K)])
        return W_result
    else:
        logging.warning(f"W optimization status: {model.status}. Returning uniform W.")
        return np.ones((K, n_markers)) / n_markers


def _update_beta(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
) -> np.ndarray:
    """
    Update beta given fixed Y and W (closed-form least squares).

    For each marker m:
        beta[m] = (X[:,m]^T @ (Y @ W[:,m])) / ((Y @ W[:,m])^T @ (Y @ W[:,m]))

    Args:
        X: Antibody data (n_spots, n_markers).
        Y: Proportions (n_spots, K).
        W: Profile weights (K, n_markers).
        beta_min: Minimum beta value.
        beta_max: Maximum beta value.

    Returns:
        Updated beta array (n_markers,).
    """
    n_markers = X.shape[1]
    beta = np.zeros(n_markers)

    for m in range(n_markers):
        # YW = Y @ W[:, m] is the predicted signal for marker m (before scaling)
        YW = Y @ W[:, m]  # (n_spots,)

        numerator = np.dot(X[:, m], YW)
        denominator = np.dot(YW, YW) + 1e-9

        beta[m] = np.clip(numerator / denominator, beta_min, beta_max)

    return beta


def _compute_reconstruction_error(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    beta: np.ndarray,
) -> float:
    """Compute reconstruction error ||X - Y @ diag(beta) @ W^T||^2_F."""
    X_pred = Y @ (W * beta[np.newaxis, :])  # Broadcasting: (n_spots, K) @ (K, n_markers)
    error = np.sum((X - X_pred) ** 2)
    return error


def _compute_bic(
    X: np.ndarray,
    Y: np.ndarray,
    W: np.ndarray,
    beta: np.ndarray,
) -> float:
    """
    Compute Bayesian Information Criterion (BIC) for model selection.

    BIC = n * log(SSE/n) + k * log(n)

    where:
        n = number of observations (n_spots * n_markers)
        SSE = reconstruction error
        k = number of parameters (K * n_markers + n_spots * K + n_markers)

    Lower BIC is better.
    """
    n_spots, n_markers = X.shape
    K = Y.shape[1]

    n = n_spots * n_markers
    sse = _compute_reconstruction_error(X, Y, W, beta)

    # Number of parameters
    # W: K * n_markers (soft weights)
    # Y: n_spots * K (proportions)
    # beta: n_markers (scaling)
    n_params = K * n_markers + n_spots * K + n_markers

    # BIC formula
    bic = n * np.log(sse / n + 1e-10) + n_params * np.log(n)

    return bic


def _alternating_minimization(
    X: np.ndarray,
    W_init: np.ndarray,
    Y_init: np.ndarray,
    beta_init: np.ndarray,
    L_coo: sp.coo_matrix,
    lambda_spatial: float,
    lambda_sparsity: float,
    lambda_distinct: float,
    max_markers_per_type: int,
    max_iterations: int,
    tolerance: float,
    verbose: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, int]:
    """
    Main alternating minimization loop.

    Args:
        X: Antibody data (n_spots, n_markers).
        W_init, Y_init, beta_init: Initial values.
        L_coo: Spatial Laplacian.
        lambda_spatial, lambda_sparsity, lambda_distinct: Regularization weights.
        max_markers_per_type: Hard constraint on markers per cell type.
        max_iterations: Maximum iterations.
        tolerance: Convergence tolerance on objective change.
        verbose: Whether to log progress.

    Returns:
        Tuple of (W, Y, beta, final_objective, n_iterations).
    """
    W = W_init.copy()
    Y = Y_init.copy()
    beta = beta_init.copy()

    prev_obj = float("inf")
    total_start = time.time()

    for iteration in range(max_iterations):
        iter_start = time.time()

        # Step 1: Optimize Y (fix W, beta)
        print(f"  Iter {iteration + 1}/{max_iterations}: Optimizing Y...", flush=True)
        y_start = time.time()
        Y = _optimize_Y(X, W, beta, L_coo, lambda_spatial, Y_prev=Y)
        y_time = time.time() - y_start

        # Step 2: Optimize W (fix Y, beta)
        print(f"  Iter {iteration + 1}/{max_iterations}: Optimizing W...", flush=True)
        w_start = time.time()
        W = _optimize_W(X, Y, beta, lambda_sparsity, lambda_distinct, max_markers_per_type)
        w_time = time.time() - w_start

        # Step 3: Update beta (closed form)
        beta = _update_beta(X, Y, W)

        # Compute objective
        obj = _compute_reconstruction_error(X, Y, W, beta)
        iter_time = time.time() - iter_start

        if verbose:
            msg = f"  Iteration {iteration + 1}: error={obj:.4f}, Y={y_time:.1f}s, W={w_time:.1f}s, total={iter_time:.1f}s"
            logging.info(msg)
            print(msg, flush=True)

        # Check convergence
        if abs(prev_obj - obj) < tolerance:
            total_time = time.time() - total_start
            logging.info(f"  Converged after {iteration + 1} iterations ({total_time:.1f}s total)")
            return W, Y, beta, obj, iteration + 1

        prev_obj = obj

    total_time = time.time() - total_start
    logging.info(f"  Reached max iterations ({max_iterations}) in {total_time:.1f}s")
    return W, Y, beta, prev_obj, max_iterations


def _format_profiles(
    W: np.ndarray,
    marker_names: List[str],
    threshold: float = 0.5,
    max_markers: int = 3,
) -> Tuple[Dict[str, Dict[str, List[str]]], Dict[int, str]]:
    """
    Convert soft W weights to CITEgeist-compatible profile dictionary.

    For each cell type (row of W), identify markers with weight > threshold.
    Name the cell type using the top markers.

    Args:
        W: Profile weights (K, n_markers).
        marker_names: List of marker names.
        threshold: Minimum weight to include marker in profile.
        max_markers: Maximum markers to include per profile (limits profile size).

    Returns:
        Tuple of:
        - Profile dictionary in CITEgeist format.
        - Mapping from cell type index k to profile name (for active profiles).
    """
    K, n_markers = W.shape
    profiles = {}
    index_to_name: Dict[int, str] = {}

    for k in range(K):
        # Get markers above threshold
        marker_indices = np.where(W[k, :] > threshold)[0]

        if len(marker_indices) == 0:
            # No markers above threshold - this cell type is inactive
            # Still assign a name for proportion tracking
            index_to_name[k] = f"Inactive_{k}"
            continue

        # Sort by weight (descending)
        sorted_indices = marker_indices[np.argsort(-W[k, marker_indices])]

        # Limit to max_markers (prevents overly large profiles)
        sorted_indices = sorted_indices[:max_markers]

        # Get marker names
        profile_markers = [marker_names[i] for i in sorted_indices]

        # Create name from top 2 markers (or fewer if less available)
        name_markers = profile_markers[:2]
        profile_name = "_".join(sorted(name_markers))

        # Handle duplicate names by appending index
        original_name = profile_name
        suffix = 1
        while profile_name in profiles:
            profile_name = f"{original_name}_{suffix}"
            suffix += 1

        profiles[profile_name] = {"Major": profile_markers}
        index_to_name[k] = profile_name

    # Add Unknown cell type (required by CITEgeist)
    profiles["Unknown"] = {"Major": []}

    return profiles, index_to_name


def optimize_profiles_jointly(
    X: np.ndarray,
    marker_names: List[str],
    coords: np.ndarray,
    min_K: int = 2,
    max_K: int = 12,
    max_markers_per_type: int = 3,
    lambda_spatial: float = 0.1,
    lambda_sparsity: float = 0.1,
    lambda_distinct: float = 0.1,  # Reduced from 0.5 (convex penalty is more aggressive)
    laplacian_k: int = 8,
    max_iterations: int = 10,  # Reduced from 50 for faster debugging
    tolerance: float = 1e-4,
    n_restarts: int = 3,
    seed: int = 1234,
    profile_threshold: float = 0.5,  # Increased from 0.3 for cleaner profiles
    verbose: bool = True,
) -> JointOptimizationResult:
    """
    Joint optimization for profile discovery and proportion estimation.

    This function simultaneously learns:
    - W[k,m]: Soft profile weights (how much marker m defines cell type k)
    - Y[s,k]: Cell type proportions at each spot
    - beta[m]: Per-marker scaling factors

    Uses Block Coordinate Descent (alternating minimization) with BIC-based K selection.

    Args:
        X: Antibody capture data (n_spots, n_markers). Should be normalized (e.g., CLR).
        marker_names: List of marker names corresponding to columns of X.
        coords: Spatial coordinates (n_spots, 2) for Laplacian smoothing.
        min_K: Minimum number of cell types to try.
        max_K: Maximum number of cell types to try.
        max_markers_per_type: Maximum markers per cell type (hard constraint).
        lambda_spatial: Weight for Laplacian spatial smoothing on Y.
        lambda_sparsity: Weight for L1 sparsity on W.
        lambda_distinct: Weight for profile distinctness (penalizes overlap).
        laplacian_k: Number of neighbors for Laplacian graph construction.
        max_iterations: Maximum alternating minimization iterations per K.
        tolerance: Convergence tolerance on objective change.
        n_restarts: Number of random restarts per K for robustness.
        seed: Random seed for reproducibility.
        profile_threshold: Minimum W weight to include marker in profile.
        verbose: Whether to log detailed progress.

    Returns:
        JointOptimizationResult with optimized W, Y, beta, profiles, and metadata.
    """
    n_spots, n_markers = X.shape

    if len(marker_names) != n_markers:
        raise ValueError(f"marker_names length ({len(marker_names)}) != n_markers ({n_markers})")

    if coords.shape[0] != n_spots:
        raise ValueError(f"coords rows ({coords.shape[0]}) != n_spots ({n_spots})")

    logging.info(f"Joint optimization: {n_spots} spots, {n_markers} markers")
    logging.info(f"K range: [{min_K}, {max_K}], restarts: {n_restarts}")

    # Handle negative values (e.g., from CLR transformation)
    # NMF-style factorization requires non-negative data
    X_min = X.min()
    if X_min < 0:
        logging.warning(
            f"Input data contains negative values (min={X_min:.4f}). "
            "Shifting to non-negative range for NMF-style factorization."
        )
        X = X - X_min  # Shift so min becomes 0
        logging.info(f"After shift: min={X.min():.4f}, max={X.max():.4f}")

    # Build spatial Laplacian (same for all K)
    L = build_spatial_laplacian(coords, k=laplacian_k, normed=True)
    L_coo = sp.coo_matrix(L)

    # Store results for each K
    results_by_K: Dict[int, Tuple[np.ndarray, np.ndarray, np.ndarray, float, float]] = {}

    for K in range(min_K, max_K + 1):
        logging.info(f"Trying K={K}...")

        best_obj_for_K = float("inf")
        best_result_for_K = None

        for restart in range(n_restarts):
            restart_seed = seed + K * 100 + restart

            # Initialize
            W_init, Y_init, beta_init = _initialize_nmf(X, K, seed=restart_seed)

            # Run alternating minimization
            W, Y, beta, obj, n_iter = _alternating_minimization(
                X,
                W_init,
                Y_init,
                beta_init,
                L_coo,
                lambda_spatial=lambda_spatial,
                lambda_sparsity=lambda_sparsity,
                lambda_distinct=lambda_distinct,
                max_markers_per_type=max_markers_per_type,
                max_iterations=max_iterations,
                tolerance=tolerance,
                verbose=verbose,
            )

            if obj < best_obj_for_K:
                best_obj_for_K = obj
                best_result_for_K = (W, Y, beta)

        # Compute BIC for best result at this K
        W, Y, beta = best_result_for_K
        bic = _compute_bic(X, Y, W, beta)
        results_by_K[K] = (W, Y, beta, best_obj_for_K, bic)

        logging.info(f"  K={K}: best_obj={best_obj_for_K:.4f}, BIC={bic:.4f}")

    # Select best K by BIC
    best_K = min(results_by_K.keys(), key=lambda k: results_by_K[k][4])
    W_best, Y_best, beta_best, obj_best, bic_best = results_by_K[best_K]

    logging.info(f"Selected K={best_K} with BIC={bic_best:.4f}")

    # Format profiles and get index-to-name mapping
    profiles, index_to_name = _format_profiles(
        W_best, marker_names, threshold=profile_threshold, max_markers=max_markers_per_type
    )

    # Collect metadata
    metadata = {
        "seed": seed,
        "min_K": min_K,
        "max_K": max_K,
        "n_restarts": n_restarts,
        "lambda_spatial": lambda_spatial,
        "lambda_sparsity": lambda_sparsity,
        "lambda_distinct": lambda_distinct,
        "max_markers_per_type": max_markers_per_type,
        "laplacian_k": laplacian_k,
        "profile_threshold": profile_threshold,
        "final_reconstruction_error": obj_best,
        "bic_trace": {k: results_by_K[k][4] for k in results_by_K},
        "n_profiles_discovered": len(profiles) - 1,  # Exclude Unknown
    }

    return JointOptimizationResult(
        W=W_best,
        Y=Y_best,
        beta=beta_best,
        profiles=profiles,
        K=best_K,
        bic=bic_best,
        marker_names=marker_names,
        index_to_name=index_to_name,
        metadata=metadata,
    )
