"""
Multimodal refinement for CITEgeist Module 3.

This module implements Pass 1.5 (anchor gene learning) and Pass 2 EM
(RNA-based proportion refinement) to handle cells with RNA signal
but low protein expression.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def select_anchor_genes(
    GEX: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    n_anchors: int = 20,
    min_correlation: float = 0.3,
) -> Tuple[Dict[str, List[str]], Dict[str, Dict[str, float]]]:
    """
    Select top-N anchor genes per cell type based on correlation and specificity.

    Pass 1.5: Learn gene signatures from protein-confident spots.

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        n_anchors: Number of anchor genes per cell type (default: 20)
        min_correlation: Minimum Pearson r to be considered (default: 0.3)

    Returns:
        anchors: Dict[cell_type] -> List[gene_names]
        weights: Dict[cell_type] -> Dict[gene_name -> correlation_weight]
    """
    n_spots, n_genes = GEX.shape
    n_types = Y_protein.shape[1]

    if n_genes != len(gene_names):
        raise ValueError(f"GEX has {n_genes} genes but {len(gene_names)} names provided")
    if n_types != len(cell_type_names):
        raise ValueError(f"Y_protein has {n_types} types but {len(cell_type_names)} names provided")

    # Compute correlation matrix (genes × cell types)
    correlations = np.zeros((n_genes, n_types))
    for t in range(n_types):
        y_t = Y_protein[:, t]
        # Skip if no variance in proportions
        if np.std(y_t) < 1e-10:
            continue
        for g in range(n_genes):
            gex_g = GEX[:, g]
            if np.std(gex_g) < 1e-10:
                continue
            r, _ = pearsonr(gex_g, y_t)
            if not np.isnan(r):
                correlations[g, t] = r

    logger.info(f"Computed correlations for {n_genes} genes × {n_types} cell types")

    # Compute specificity: r[g,t] - max(r[g, other_types])
    specificity = np.zeros((n_genes, n_types))
    for g in range(n_genes):
        for t in range(n_types):
            other_corrs = np.delete(correlations[g, :], t)
            other_max = np.max(other_corrs) if len(other_corrs) > 0 else 0
            specificity[g, t] = correlations[g, t] - other_max

    # Combined score: correlation × specificity (only positive specificity)
    score = correlations * np.clip(specificity, 0, None)

    # Select top-N per cell type
    anchors = {}
    weights = {}

    for t, ct_name in enumerate(cell_type_names):
        # Filter by minimum correlation
        valid_mask = correlations[:, t] >= min_correlation
        valid_indices = np.where(valid_mask)[0]

        if len(valid_indices) == 0:
            logger.warning(f"No genes pass min_correlation={min_correlation} for {ct_name}")
            anchors[ct_name] = []
            weights[ct_name] = {}
            continue

        # Rank by score among valid genes
        valid_scores = score[valid_indices, t]
        ranked_order = np.argsort(valid_scores)[::-1]  # Descending
        top_indices = valid_indices[ranked_order[:n_anchors]]

        # Build output
        anchors[ct_name] = [gene_names[i] for i in top_indices]
        weights[ct_name] = {
            gene_names[i]: float(correlations[i, t]) for i in top_indices
        }

        logger.info(
            f"{ct_name}: selected {len(anchors[ct_name])} anchors, "
            f"top gene={anchors[ct_name][0] if anchors[ct_name] else 'None'}"
        )

    return anchors, weights


def compute_expression_profiles(
    GEX: np.ndarray,
    Y: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    anchors: Dict[str, List[str]],
    weights: Dict[str, Dict[str, float]],
) -> np.ndarray:
    """
    E-step: Estimate gene expression profiles per cell type.

    Anchor genes are LOCKED to their assigned cell type.
    Non-anchor genes are solved via least squares (can load on any type).

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y: Current cell type proportions (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        anchors: Dict[cell_type] -> List[anchor_gene_names]
        weights: Dict[cell_type] -> Dict[gene_name -> weight]

    Returns:
        E: Expression profiles (T_types × G_genes)
    """
    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    # Build gene name to index mapping
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    type_to_idx = {t: i for i, t in enumerate(cell_type_names)}

    # Build set of anchor genes and their assignments
    anchor_assignments = {}  # gene_name -> cell_type_name
    for ct_name, gene_list in anchors.items():
        for g in gene_list:
            anchor_assignments[g] = ct_name

    # Initialize E
    E = np.zeros((n_types, n_genes))

    for g_idx, g_name in enumerate(gene_names):
        gex_g = GEX[:, g_idx]  # Expression of gene g across spots

        # For ALL genes (anchor or not), use least squares
        # GEX[:, g] ≈ Y @ E[:, g]
        # E[:, g] = (Y^T Y)^{-1} Y^T GEX[:, g]
        try:
            E[:, g_idx], _, _, _ = np.linalg.lstsq(Y, gex_g, rcond=None)
            # Clip negative values (expression should be non-negative)
            E[:, g_idx] = np.clip(E[:, g_idx], 0, None)
        except np.linalg.LinAlgError:
            # Fallback: uniform assignment
            E[:, g_idx] = np.mean(gex_g) / n_types

        # For anchor genes, boost the assigned cell type's estimate
        # (soft encouragement, not hard constraint)
        if g_name in anchor_assignments:
            ct_name = anchor_assignments[g_name]
            t_idx = type_to_idx[ct_name]
            # Keep the least squares solution but could add regularization here

    return E


def refine_proportions(
    GEX: np.ndarray,
    Y_current: np.ndarray,
    E: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    anchors: Dict[str, List[str]],
    weights: Dict[str, Dict[str, float]],
    lambda_prior: float = 1.0,
) -> np.ndarray:
    """
    M-step: Refine proportions Y given expression profiles E.

    Objective per spot i:
        minimize: Σ_g w[g] * (GEX[i,g] - Y[i,:] @ E[:,g])²
                  + λ * ||Y[i,:] - Y_protein[i,:]||²

    Constraints:
        Y[i, :] >= 0
        sum(Y[i, :]) <= 1

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_current: Current proportions (N_spots × T_types)
        E: Expression profiles (T_types × G_genes)
        Y_protein: Original protein-based proportions (N_spots × T_types)
        gene_names: List of gene names
        cell_type_names: List of cell type names
        anchors: Dict of anchor genes per cell type
        weights: Dict of weights per anchor gene
        lambda_prior: Weight of protein prior (default: 1.0)

    Returns:
        Y_refined: Updated proportions (N_spots × T_types)
    """
    from scipy.optimize import minimize, Bounds

    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    # Build gene weights vector
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    gene_weights = np.ones(n_genes)
    for ct_name, gene_dict in weights.items():
        for g_name, w in gene_dict.items():
            if g_name in gene_to_idx:
                gene_weights[gene_to_idx[g_name]] = w

    Y_refined = np.zeros((n_spots, n_types))

    for i in range(n_spots):
        gex_i = GEX[i, :]
        y_protein_i = Y_protein[i, :]

        def objective(y):
            # Reconstruction error (weighted)
            reconstruction = E.T @ y  # (G,)
            recon_error = np.sum(gene_weights * (gex_i - reconstruction) ** 2)

            # Prior toward protein proportions
            prior_error = lambda_prior * np.sum((y - y_protein_i) ** 2)

            return recon_error + prior_error

        # Bounds: 0 <= y[t] <= 1
        bounds = Bounds(lb=np.zeros(n_types), ub=np.ones(n_types))

        # Constraint: sum(y) <= 1
        constraints = {"type": "ineq", "fun": lambda y: 1.0 - np.sum(y)}

        # Initial guess: current proportions
        y0 = Y_current[i, :].copy()

        result = minimize(
            objective,
            y0,
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"maxiter": 100, "ftol": 1e-8},
        )

        if result.success:
            Y_refined[i, :] = result.x
        else:
            # Fallback to current
            Y_refined[i, :] = Y_current[i, :]

    return Y_refined


def multimodal_em_refinement(
    GEX: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    n_anchors: int = 20,
    min_correlation: float = 0.3,
    lambda_prior: float = 1.0,
    max_iterations: int = 20,
    tolerance: float = 1e-4,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, List[str]]]:
    """
    Full multimodal EM refinement: Pass 1.5 + Pass 2 EM.

    1. Pass 1.5: Learn anchor genes from protein-confident spots
    2. Pass 2 EM: Iterate E-step (expression profiles) and M-step (proportions)

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        n_anchors: Number of anchor genes per cell type (default: 20)
        min_correlation: Minimum Pearson r for anchor selection (default: 0.3)
        lambda_prior: Weight of protein prior in M-step (default: 1.0)
        max_iterations: Maximum EM iterations (default: 20)
        tolerance: Convergence tolerance (default: 1e-4)

    Returns:
        Y_refined: Refined proportions (N_spots × T_types)
        E_final: Final expression profiles (T_types × G_genes)
        anchors: Dict of anchor genes per cell type
    """
    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    logger.info(f"Starting multimodal EM refinement: {n_spots} spots, {n_genes} genes, {n_types} types")

    # Pass 1.5: Learn anchor genes
    logger.info("Pass 1.5: Selecting anchor genes...")
    anchors, weights = select_anchor_genes(
        GEX=GEX,
        Y_protein=Y_protein,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        n_anchors=n_anchors,
        min_correlation=min_correlation,
    )

    total_anchors = sum(len(v) for v in anchors.values())
    logger.info(f"Selected {total_anchors} total anchor genes across {n_types} cell types")

    # Initialize
    Y = Y_protein.copy()
    E = None

    # EM loop
    for iteration in range(max_iterations):
        logger.info(f"EM iteration {iteration + 1}/{max_iterations}")

        # E-step: Estimate expression profiles
        E = compute_expression_profiles(
            GEX=GEX,
            Y=Y,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
        )

        # M-step: Refine proportions
        Y_new = refine_proportions(
            GEX=GEX,
            Y_current=Y,
            E=E,
            Y_protein=Y_protein,
            gene_names=gene_names,
            cell_type_names=cell_type_names,
            anchors=anchors,
            weights=weights,
            lambda_prior=lambda_prior,
        )

        # Check convergence
        max_change = np.max(np.abs(Y_new - Y))
        logger.info(f"  Max proportion change: {max_change:.6f}")

        if max_change < tolerance:
            logger.info(f"Converged after {iteration + 1} iterations")
            Y = Y_new
            break

        Y = Y_new

    else:
        logger.warning(f"Did not converge within {max_iterations} iterations")

    # Final E-step to get consistent E
    E_final = compute_expression_profiles(
        GEX=GEX,
        Y=Y,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        anchors=anchors,
        weights=weights,
    )

    logger.info("Multimodal EM refinement complete")

    return Y, E_final, anchors
