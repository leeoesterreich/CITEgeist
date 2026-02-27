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

    # Initialize beta
    if beta_init is None:
        beta = np.ones(n_markers)
    else:
        beta = beta_init.copy()

    # Block coordinate descent
    counts = None
    alpha = None

    for iteration in range(max_iter):
        logger.debug(f"BCD iteration {iteration + 1}/{max_iter}")

        # Block 1: Fix beta, solve IQP for counts and alpha
        counts, alpha = _solve_iqp_fixed_beta(
            X, nuclei_counts, profile, detected, weights, beta, timeout
        )

        # Block 2: Fix counts, solve OLS for alpha and beta
        alpha_new, beta_new = _solve_ols_fixed_counts(
            X, counts, profile, weights, beta_min
        )

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

    With beta fixed, the problem is:
        min sum_i sum_m w[m] * (X[i,m] - alpha[m] - sum_k c[i,k]*profile[k,m]*beta[m])^2

    This is quadratic in (c, alpha) since beta is constant.
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


def _solve_ols_fixed_counts(
    X: np.ndarray,
    counts: np.ndarray,
    profile: np.ndarray,
    weights: np.ndarray,
    beta_min: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve OLS for alpha and beta with fixed counts.

    With counts fixed, the model is linear:
        X[i,m] = alpha[m] + sum_k counts[i,k] * profile[k,m] * beta[m] + noise

    This is equivalent to per-marker regression:
        X[:,m] = alpha[m] + beta[m] * (counts @ profile[:,m])
    """
    n_spots, n_markers = X.shape

    alpha = np.zeros(n_markers)
    beta = np.zeros(n_markers)

    for m in range(n_markers):
        # Effective cell count for this marker: sum_k counts[i,k] * profile[k,m]
        effective_counts = counts @ profile[:, m]  # (n_spots,)

        # Weighted least squares: X[:,m] = alpha_m + beta_m * effective_counts
        # Design matrix: [1, effective_counts]
        A = np.column_stack([np.ones(n_spots), effective_counts])
        b = X[:, m]

        # Weighted normal equations
        W = np.diag(np.full(n_spots, weights[m]))
        AtWA = A.T @ W @ A
        AtWb = A.T @ W @ b

        try:
            params = np.linalg.solve(AtWA, AtWb)
            alpha[m] = max(0, params[0])  # alpha >= 0
            beta[m] = max(beta_min, params[1])  # beta >= beta_min
        except np.linalg.LinAlgError:
            # Fallback to simple estimates
            alpha[m] = X[:, m].min()
            beta[m] = beta_min

    return alpha, beta
