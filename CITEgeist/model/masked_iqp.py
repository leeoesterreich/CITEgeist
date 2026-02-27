# CITEgeist/model/masked_iqp.py
"""Masked Integer Quadratic Programming solver for cell count estimation.

This module implements Stage 2 of the two-stage detection + estimation model.
Given a detection mask from Stage 1, it solves for integer cell counts while:
- Enforcing true zeros for non-detected cell types
- Learning baseline alpha[m] (background signal per marker)
- Learning beta[m] (signal-per-cell per marker)
- Respecting nuclei count sum constraint

Uses block coordinate descent:
1. Fix beta → solve IQP for counts and alpha
2. Fix counts → solve OLS for alpha and beta
"""
import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


def estimate_beta_direct(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    detected: np.ndarray,
    min_single_spots: int = 3,
) -> np.ndarray:
    """
    Estimate beta directly from single-type detected spots.

    For spots where ONLY one cell type is detected, all nuclei belong to that type.
    This allows direct estimation: beta_m = (X_m - alpha_m) / nuclei_count

    This avoids the chicken-and-egg problem where OLS learns beta from bad IQP
    allocations, which then collapse beta to zero for common cell types.

    Args:
        X: (n_spots, n_markers) observed antibody signal.
        nuclei_counts: (n_spots,) integer nuclei count per spot.
        profile: (n_types, n_markers) binary matrix.
        detected: (n_spots, n_types) boolean detection mask.
        min_single_spots: Minimum single-type spots needed for direct estimate.

    Returns:
        beta_init: (n_markers,) initial beta estimates.
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    beta = np.ones(n_markers)  # default initialization
    alpha = np.zeros(n_markers)

    # Count how many types are detected per spot
    n_detected_per_spot = detected.sum(axis=1)

    for k in range(n_types):
        # Spots where ONLY this type is detected
        single_detected = detected[:, k] & (n_detected_per_spot == 1)
        not_detected = ~detected[:, k]

        # Markers used by this cell type
        marker_indices = np.where(profile[k] > 0)[0]

        if single_detected.sum() < min_single_spots:
            # Not enough single-type spots - use dynamic range heuristic
            for m in marker_indices:
                if not_detected.sum() > 10:
                    # Estimate alpha from non-detected spots
                    alpha[m] = np.percentile(X[not_detected, m], 50)
                    # Estimate beta from dynamic range
                    detected_signal = X[detected[:, k], m]
                    if len(detected_signal) > 0:
                        signal_above = np.maximum(detected_signal - alpha[m], 0)
                        # Assume average ~5 cells per spot for scaling
                        beta[m] = max(1.0, np.percentile(signal_above, 75) / 5)
            continue

        # Direct estimation from single-type spots
        for m in marker_indices:
            # Alpha from non-detected spots
            if not_detected.sum() > 10:
                alpha[m] = np.percentile(X[not_detected, m], 50)
            else:
                alpha[m] = np.percentile(X[:, m], 10)

            # Beta = (signal - alpha) / nuclei_count in single-type spots
            signal = X[single_detected, m]
            n_cells = nuclei_counts[single_detected]

            per_cell_signal = (signal - alpha[m]) / np.maximum(n_cells, 1)
            valid = per_cell_signal > 0

            if valid.sum() > 0:
                beta[m] = np.median(per_cell_signal[valid])
            else:
                # Fallback to dynamic range
                beta[m] = max(1.0, (signal.max() - alpha[m]) / max(n_cells.max(), 1))

    logger.info(f"Direct beta initialization: {beta}")
    return beta


def solve_masked_iqp(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    detected: np.ndarray,
    weights: np.ndarray,
    beta_init: np.ndarray = None,
    beta_min: float = 1e-3,
    timeout: float = 300.0,
    max_iter: int = 10,
    convergence_tol: float = 0.01,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Solve IQP for cell counts with detection mask using block coordinate descent.

    Block 1 (IQP): Fix beta, solve for counts and alpha
    Block 2 (OLS): Fix counts, solve for alpha and beta via least squares

    Minimizes weighted sum of squared residuals:
        sum_i sum_m w[m] * (X[i,m] - alpha[m] - sum_k c[i,k]*profile[k,m]*beta[m])^2

    Subject to:
        c[i,k] = 0 if detected[i,k] = False
        c[i,k] in {0,1,...,N_i} if detected[i,k] = True
        sum_k c[i,k] = N_i for spots with detected types
        alpha[m] >= 0
        beta[m] >= beta_min

    Args:
        X: (n_spots, n_markers) observed antibody signal.
        nuclei_counts: (n_spots,) integer nuclei count per spot.
        profile: (n_types, n_markers) binary matrix where profile[k,m]=1
            if marker m defines cell type k.
        detected: (n_spots, n_types) boolean detection mask from Stage 1.
        weights: (n_markers,) inverse variance weights for each marker.
        beta_init: (n_markers,) initial beta values. If None, initialize to 1.0.
        beta_min: Minimum value for beta (signal-per-cell), default 1e-3.
        timeout: Solver timeout in seconds per IQP solve.
        max_iter: Maximum block coordinate descent iterations.
        convergence_tol: Relative change in beta to declare convergence.

    Returns:
        counts: (n_spots, n_types) integer cell counts.
        alpha: (n_markers,) learned baseline per marker.
        beta: (n_markers,) learned signal-per-cell per marker.
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    # Initialize beta - use direct estimation if not provided
    if beta_init is None:
        beta = estimate_beta_direct(X, nuclei_counts, profile, detected)
    else:
        beta = beta_init.copy()

    # Save initial beta as FLOOR - OLS can refine upward but NEVER collapse below
    # This is the key robustness mechanism: no hyperparameter tuning needed
    beta_floor = beta.copy()

    # Block coordinate descent
    counts = None
    alpha = None

    for iteration in range(max_iter):
        logger.debug(f"BCD iteration {iteration + 1}/{max_iter}")

        # Block 1: Fix beta, solve IQP for counts and alpha
        counts, alpha = _solve_iqp_fixed_beta(
            X, nuclei_counts, profile, detected, weights, beta, timeout
        )

        # Block 2: Fix counts, update alpha and beta
        # Use simple OLS (no ridge needed since we have floor constraint)
        alpha_new, beta_ols = _solve_ols_fixed_counts(
            X, counts, profile, weights, beta_min,
            beta_prior=None,  # No ridge regularization
            ridge_lambda=0.0,  # Disabled
        )

        # FLOOR CONSTRAINT: Beta can only go UP from direct estimate
        # This prevents collapse without requiring hyperparameter tuning
        beta_new = np.maximum(beta_ols, beta_floor)

        # Check convergence
        beta_change = np.abs(beta_new - beta).max() / (np.abs(beta).max() + 1e-8)
        logger.debug(f"  beta change: {beta_change:.4f}")

        alpha = alpha_new
        beta = beta_new

        if beta_change < convergence_tol:
            logger.debug(f"Converged at iteration {iteration + 1}")
            break

    return counts, alpha, beta


def _solve_iqp_fixed_beta(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    detected: np.ndarray,
    weights: np.ndarray,
    beta: np.ndarray,
    timeout: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve IQP for counts and alpha with fixed beta.

    Uses per-spot optimization for speed (each spot is independent given alpha).
    Alpha is estimated from residuals after count optimization.
    """
    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    # Precompute effective signal contribution: profile[k,m] * beta[m]
    signal_coeff = profile * beta  # (n_types, n_markers)

    # Initialize counts
    counts = np.zeros((n_spots, n_types), dtype=int)

    # Estimate initial alpha from spots with minimal signal
    # (where few types are detected)
    n_detected_per_spot = detected.sum(axis=1)
    background_mask = n_detected_per_spot <= 1
    if background_mask.sum() > 10:
        alpha = np.percentile(X[background_mask], 10, axis=0)
    else:
        alpha = np.percentile(X, 5, axis=0)
    alpha = np.maximum(alpha, 0)

    # Solve each spot independently
    for i in range(n_spots):
        detected_types = np.where(detected[i])[0]
        n_i = int(nuclei_counts[i])

        if len(detected_types) == 0 or n_i == 0:
            continue

        if len(detected_types) == 1:
            # Only one type detected - assign all nuclei to it
            counts[i, detected_types[0]] = n_i
            continue

        # Multiple types detected - solve small IQP
        counts[i] = _solve_single_spot_iqp(
            X[i], n_i, signal_coeff, detected_types, weights, alpha, timeout=5.0
        )

    # Update alpha estimate from residuals
    # For each marker, alpha = median(X - predicted) where predicted uses counts
    predicted = (counts @ profile) * beta
    residuals = X - predicted
    alpha = np.maximum(np.median(residuals, axis=0), 0)

    return counts, alpha


def _solve_single_spot_iqp(
    x: np.ndarray,
    n_cells: int,
    signal_coeff: np.ndarray,
    detected_types: np.ndarray,
    weights: np.ndarray,
    alpha: np.ndarray,
    timeout: float = 5.0,
) -> np.ndarray:
    """
    Solve IQP for a single spot.

    Much faster than global IQP since only ~7 integer variables per spot.
    """
    import gurobipy as gp
    from gurobipy import GRB

    n_types = signal_coeff.shape[0]
    n_markers = len(x)
    n_detected = len(detected_types)

    # For very small problems, enumerate solutions
    if n_detected <= 2 and n_cells <= 20:
        return _enumerate_best_counts(x, n_cells, signal_coeff, detected_types, weights, alpha)

    model = gp.Model("spot_iqp")
    model.setParam("OutputFlag", 0)
    model.setParam("TimeLimit", timeout)

    # Variables: counts for detected types only
    c = model.addVars(n_detected, vtype=GRB.INTEGER, lb=0, ub=n_cells, name="c")

    model.update()

    # Objective: weighted sum of squared residuals
    obj = gp.QuadExpr()
    for m in range(n_markers):
        # predicted = alpha[m] + sum over detected types of c[k] * signal_coeff[type,m]
        pred = alpha[m]
        for idx, k in enumerate(detected_types):
            pred = pred + c[idx] * signal_coeff[k, m]

        residual = x[m] - pred
        obj += weights[m] * residual * residual

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraint: counts sum to n_cells
    model.addConstr(gp.quicksum(c[idx] for idx in range(n_detected)) == n_cells, "sum")

    model.optimize()

    # Extract solution
    counts = np.zeros(n_types, dtype=int)
    if model.status in [GRB.OPTIMAL, GRB.TIME_LIMIT]:
        for idx, k in enumerate(detected_types):
            counts[k] = int(round(c[idx].X))

        # Ensure sum constraint (rounding may break it)
        diff = n_cells - counts.sum()
        if diff != 0:
            # Adjust largest count
            largest_idx = detected_types[np.argmax(counts[detected_types])]
            counts[largest_idx] += diff
    else:
        # Fallback: distribute evenly
        per_type = n_cells // n_detected
        remainder = n_cells % n_detected
        for idx, k in enumerate(detected_types):
            counts[k] = per_type + (1 if idx < remainder else 0)

    return counts


def _enumerate_best_counts(
    x: np.ndarray,
    n_cells: int,
    signal_coeff: np.ndarray,
    detected_types: np.ndarray,
    weights: np.ndarray,
    alpha: np.ndarray,
) -> np.ndarray:
    """
    For small problems, enumerate all valid count combinations.
    """
    n_types = signal_coeff.shape[0]
    n_detected = len(detected_types)

    best_counts = np.zeros(n_types, dtype=int)
    best_cost = float('inf')

    if n_detected == 1:
        best_counts[detected_types[0]] = n_cells
        return best_counts

    # Two types: enumerate all splits
    for c0 in range(n_cells + 1):
        c1 = n_cells - c0

        # Compute cost
        counts_vec = np.zeros(n_types)
        counts_vec[detected_types[0]] = c0
        counts_vec[detected_types[1]] = c1

        predicted = alpha + counts_vec @ signal_coeff
        residuals = x - predicted
        cost = (weights * residuals ** 2).sum()

        if cost < best_cost:
            best_cost = cost
            best_counts = counts_vec.astype(int)

    return best_counts


def _solve_iqp_fixed_beta_DEPRECATED(
    X: np.ndarray,
    nuclei_counts: np.ndarray,
    profile: np.ndarray,
    detected: np.ndarray,
    weights: np.ndarray,
    beta: np.ndarray,
    timeout: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    DEPRECATED: Global IQP solver - too slow for large datasets.
    Kept for reference.
    """
    import gurobipy as gp
    from gurobipy import GRB

    n_spots, n_markers = X.shape
    n_types = profile.shape[0]

    model = gp.Model("iqp_fixed_beta")
    model.setParam("OutputFlag", 0)
    model.setParam("TimeLimit", timeout)

    # Variables: cell counts (integer where detected)
    c = {}
    for i in range(n_spots):
        for k in range(n_types):
            if detected[i, k]:
                c[i, k] = model.addVar(
                    vtype=GRB.INTEGER,
                    lb=0,
                    ub=int(nuclei_counts[i]),
                    name=f"c_{i}_{k}"
                )

    # Variables: alpha (baseline)
    alpha_vars = model.addVars(n_markers, lb=0, name="alpha")

    model.update()

    # Precompute effective signal contribution: profile[k,m] * beta[m]
    signal_coeff = profile * beta  # (n_types, n_markers)

    # Objective: weighted sum of squared residuals
    obj = gp.QuadExpr()
    for i in range(n_spots):
        for m in range(n_markers):
            # predicted = alpha[m] + sum_k c[i,k] * signal_coeff[k,m]
            pred = alpha_vars[m]
            for k in range(n_types):
                if detected[i, k] and signal_coeff[k, m] > 0:
                    pred = pred + c[i, k] * signal_coeff[k, m]

            residual = X[i, m] - pred
            obj += weights[m] * residual * residual

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraints: cell counts sum to nuclei count
    for i in range(n_spots):
        n_i = int(nuclei_counts[i])
        detected_types = [k for k in range(n_types) if detected[i, k]]

        if n_i > 0 and detected_types:
            model.addConstr(
                gp.quicksum(c[i, k] for k in detected_types) == n_i,
                name=f"sum_{i}"
            )

    # Solve
    model.optimize()

    if model.status != GRB.OPTIMAL:
        logger.warning(f"IQP solver status: {model.status}")

    # Extract solution
    counts = np.zeros((n_spots, n_types), dtype=int)
    for i in range(n_spots):
        for k in range(n_types):
            if detected[i, k]:
                counts[i, k] = int(round(c[i, k].X))

    alpha = np.array([alpha_vars[m].X for m in range(n_markers)])

    return counts, alpha


def _solve_alpha_fixed_counts_beta(
    X: np.ndarray,
    counts: np.ndarray,
    profile: np.ndarray,
    beta: np.ndarray,
) -> np.ndarray:
    """
    Solve for alpha (baseline) with fixed counts and beta.

    Model: X[i,m] = alpha[m] + (counts @ profile)[i,m] * beta[m]
    Solution: alpha[m] = median(X[:,m] - predicted[:,m])
    """
    # Compute predicted signal from counts
    predicted = (counts @ profile) * beta  # (n_spots, n_markers)

    # Alpha is the median residual (robust to outliers)
    residuals = X - predicted
    alpha = np.maximum(np.median(residuals, axis=0), 0)

    return alpha


def _solve_ols_fixed_counts(
    X: np.ndarray,
    counts: np.ndarray,
    profile: np.ndarray,
    weights: np.ndarray,
    beta_min: float,
    beta_prior: np.ndarray = None,
    ridge_lambda: float = 1.0,  # Strong regularization to prevent collapse
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve OLS for alpha and beta with fixed counts.

    With counts fixed, the model is linear:
        X[i,m] = alpha[m] + sum_k counts[i,k] * profile[k,m] * beta[m] + noise

    This is equivalent to per-marker regression:
        X[:,m] = alpha[m] + beta[m] * (counts @ profile[:,m])

    Uses ridge regularization toward beta_prior to prevent collapse.
    """
    n_spots, n_markers = X.shape

    alpha = np.zeros(n_markers)
    beta = np.zeros(n_markers)

    # Default prior is uniform
    if beta_prior is None:
        beta_prior = np.ones(n_markers)

    for m in range(n_markers):
        # Effective cell count for this marker: sum_k counts[i,k] * profile[k,m]
        effective_counts = counts @ profile[:, m]  # (n_spots,)

        # Skip markers not used by any detected cell type
        if effective_counts.sum() == 0:
            beta[m] = beta_prior[m]
            alpha[m] = np.median(X[:, m])
            continue

        # Weighted least squares: X[:,m] = alpha_m + beta_m * effective_counts
        # Design matrix: [1, effective_counts]
        A = np.column_stack([np.ones(n_spots), effective_counts])
        b = X[:, m]

        # Weighted normal equations with ridge regularization on beta
        # Regularize toward beta_prior to prevent collapse
        W = np.diag(np.full(n_spots, weights[m]))
        AtWA = A.T @ W @ A

        # Add ridge penalty on beta (second parameter), scaled by data magnitude
        ridge_penalty = ridge_lambda * np.diag([0, 1]) * AtWA[1, 1]
        AtWA_reg = AtWA + ridge_penalty

        # Adjust RHS to regularize toward prior
        AtWb = A.T @ W @ b
        AtWb[1] += ridge_lambda * AtWA[1, 1] * beta_prior[m]

        try:
            params = np.linalg.solve(AtWA_reg, AtWb)
            alpha[m] = max(0, params[0])  # alpha >= 0
            beta[m] = max(beta_min, params[1])  # beta >= beta_min
        except np.linalg.LinAlgError:
            # Fallback to prior estimates
            alpha[m] = np.median(X[:, m])
            beta[m] = beta_prior[m]

    return alpha, beta
