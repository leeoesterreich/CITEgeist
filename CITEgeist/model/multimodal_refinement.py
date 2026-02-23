"""
Multimodal refinement for CITEgeist Module 3.

This module implements Pass 1.5 (anchor gene learning) and Pass 2 EM
(RNA-based proportion refinement) to handle cells with RNA signal
but low protein expression.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.stats import spearmanr

logger = logging.getLogger(__name__)


def select_anchor_genes(
    GEX: np.ndarray,
    Y_protein: np.ndarray,
    gene_names: List[str],
    cell_type_names: List[str],
    n_anchors: int = 20,
    min_correlation: float = 0.3,
    min_anchors_per_type: int = 5,
    min_expressing_spots: int = 20,
    sparse_aware: bool = True,
    threshold_step: float = 0.05,
    min_threshold: float = 0.05,
) -> Tuple[Dict[str, List[str]], Dict[str, Dict[str, float]]]:
    """
    Select top-N anchor genes per cell type based on Spearman correlation and specificity.

    Pass 1.5: Learn gene signatures from protein-confident spots.

    Uses Spearman (rank-based) correlation for robustness to outliers.
    Adaptive thresholding ensures minimum anchor coverage per cell type.

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        n_anchors: Maximum anchor genes per cell type (default: 20, cap)
        min_correlation: Initial minimum Spearman rho to be considered (default: 0.3)
        min_anchors_per_type: Minimum anchors required per cell type (default: 5, floor)
        min_expressing_spots: Minimum spots with non-zero expression (default: 20)
        sparse_aware: If True, compute correlations only on expressing spots (default: True)
        threshold_step: Amount to lower threshold when floor not met (default: 0.05)
        min_threshold: Absolute minimum threshold to try (default: 0.05)

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

    # Check sparsity
    sparsity = (GEX == 0).sum() / GEX.size
    logger.info(f"GEX sparsity: {sparsity:.1%}, sparse_aware={sparse_aware}, using Spearman correlation")

    # Compute correlation matrix (genes × cell types) using Spearman
    correlations = np.zeros((n_genes, n_types))
    n_expressing = np.zeros(n_genes, dtype=int)  # Track expressing spots per gene

    for g in range(n_genes):
        gex_g = GEX[:, g]

        if sparse_aware:
            # Only use spots where gene is expressed (non-zero)
            expressing_mask = gex_g > 0
            n_expr = expressing_mask.sum()
            n_expressing[g] = n_expr

            if n_expr < min_expressing_spots:
                continue  # Skip genes with too few expressing spots

            gex_g_subset = gex_g[expressing_mask]
            if np.std(gex_g_subset) < 1e-10:
                continue

            for t in range(n_types):
                y_t_subset = Y_protein[expressing_mask, t]
                if np.std(y_t_subset) < 1e-10:
                    continue
                rho, _ = spearmanr(gex_g_subset, y_t_subset)
                if not np.isnan(rho):
                    correlations[g, t] = rho
        else:
            # Original behavior: use all spots
            n_expressing[g] = n_spots
            if np.std(gex_g) < 1e-10:
                continue

            for t in range(n_types):
                y_t = Y_protein[:, t]
                if np.std(y_t) < 1e-10:
                    continue
                rho, _ = spearmanr(gex_g, y_t)
                if not np.isnan(rho):
                    correlations[g, t] = rho

    n_valid_genes = (n_expressing >= min_expressing_spots).sum()
    logger.info(f"Computed Spearman correlations for {n_valid_genes}/{n_genes} genes × {n_types} cell types")

    # Log correlation statistics per cell type
    for t, ct_name in enumerate(cell_type_names):
        ct_corrs = correlations[:, t]
        positive_corrs = ct_corrs[ct_corrs > 0]
        if len(positive_corrs) > 0:
            logger.info(
                f"  {ct_name}: max_rho={ct_corrs.max():.3f}, "
                f"genes with rho>0.1: {(ct_corrs > 0.1).sum()}, "
                f"genes with rho>0.2: {(ct_corrs > 0.2).sum()}, "
                f"genes with rho>0.3: {(ct_corrs > 0.3).sum()}"
            )

    # Compute specificity: r[g,t] - max(r[g, other_types])
    specificity = np.zeros((n_genes, n_types))
    for g in range(n_genes):
        for t in range(n_types):
            other_corrs = np.delete(correlations[g, :], t)
            other_max = np.max(other_corrs) if len(other_corrs) > 0 else 0
            specificity[g, t] = correlations[g, t] - other_max

    # Combined score: correlation × specificity (only positive specificity)
    score = correlations * np.clip(specificity, 0, None)

    # Select anchors with adaptive thresholding
    anchors = {}
    weights = {}

    for t, ct_name in enumerate(cell_type_names):
        # Start with initial threshold, lower if needed to meet floor
        current_threshold = min_correlation
        selected_indices = None

        while current_threshold >= min_threshold:
            # Filter by current threshold
            valid_mask = correlations[:, t] >= current_threshold
            valid_indices = np.where(valid_mask)[0]

            if len(valid_indices) >= min_anchors_per_type:
                # We have enough, rank by score and take top n_anchors
                valid_scores = score[valid_indices, t]
                ranked_order = np.argsort(valid_scores)[::-1]  # Descending
                selected_indices = valid_indices[ranked_order[:n_anchors]]

                if current_threshold < min_correlation:
                    logger.info(
                        f"{ct_name}: lowered threshold to {current_threshold:.2f} "
                        f"to get {len(selected_indices)} anchors (min={min_anchors_per_type})"
                    )
                break

            # Not enough anchors, lower threshold
            current_threshold -= threshold_step

        # If we still don't have enough even at min_threshold, take what we can
        if selected_indices is None:
            valid_mask = correlations[:, t] > 0  # Any positive correlation
            valid_indices = np.where(valid_mask)[0]

            if len(valid_indices) > 0:
                valid_scores = score[valid_indices, t]
                ranked_order = np.argsort(valid_scores)[::-1]
                selected_indices = valid_indices[ranked_order[:n_anchors]]
                logger.warning(
                    f"{ct_name}: only {len(selected_indices)} genes with positive correlation "
                    f"(floor={min_anchors_per_type} not met)"
                )
            else:
                logger.warning(f"{ct_name}: no genes with positive correlation!")
                anchors[ct_name] = []
                weights[ct_name] = {}
                continue

        # Build output
        anchors[ct_name] = [gene_names[i] for i in selected_indices]
        weights[ct_name] = {
            gene_names[i]: float(correlations[i, t]) for i in selected_indices
        }

        logger.info(
            f"{ct_name}: selected {len(anchors[ct_name])} anchors (cap={n_anchors}, floor={min_anchors_per_type}), "
            f"top gene={anchors[ct_name][0]} (rho={weights[ct_name][anchors[ct_name][0]]:.3f})"
        )

    total_anchors = sum(len(v) for v in anchors.values())
    logger.info(f"Total anchors selected: {total_anchors} across {n_types} cell types")

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

        if g_name in anchor_assignments:
            # LOCKED: anchor gene assigned to one cell type only
            ct_name = anchor_assignments[g_name]
            t_idx = type_to_idx[ct_name]

            # Weighted mean of expression where cell type is present
            y_t = Y[:, t_idx]
            if np.sum(y_t) > 1e-10:
                E[t_idx, g_idx] = np.sum(gex_g * y_t) / np.sum(y_t)
            else:
                E[t_idx, g_idx] = 0.0

            # Other cell types get 0 for this anchor gene (already initialized to 0)

        else:
            # FREE: non-anchor gene can load on any cell type via least squares
            try:
                E[:, g_idx], _, _, _ = np.linalg.lstsq(Y, gex_g, rcond=None)
                E[:, g_idx] = np.clip(E[:, g_idx], 0, None)
            except np.linalg.LinAlgError:
                E[:, g_idx] = np.mean(gex_g) / n_types

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
    min_anchors_per_type: int = 5,
    lambda_prior: float = 1.0,
    max_iterations: int = 20,
    tolerance: float = 1e-4,
    sparse_aware: bool = True,
    min_expressing_spots: int = 20,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, List[str]]]:
    """
    Full multimodal EM refinement: Pass 1.5 + Pass 2 EM.

    1. Pass 1.5: Learn anchor genes from protein-confident spots (Spearman correlation)
    2. Pass 2 EM: Iterate E-step (expression profiles) and M-step (proportions)

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        gene_names: List of gene names (length G)
        cell_type_names: List of cell type names (length T)
        n_anchors: Maximum anchor genes per cell type (default: 20, cap)
        min_correlation: Initial minimum Spearman rho for anchor selection (default: 0.3)
        min_anchors_per_type: Minimum anchors required per cell type (default: 5, floor)
        lambda_prior: Weight of protein prior in M-step (default: 1.0)
        max_iterations: Maximum EM iterations (default: 20)
        tolerance: Convergence tolerance (default: 1e-4)
        sparse_aware: If True, compute correlations only on expressing spots (default: True)
        min_expressing_spots: Minimum spots with expression for anchor selection (default: 20)

    Returns:
        Y_refined: Refined proportions (N_spots × T_types)
        E_final: Final expression profiles (T_types × G_genes)
        anchors: Dict of anchor genes per cell type
    """
    n_spots, n_genes = GEX.shape
    n_types = len(cell_type_names)

    logger.info(f"Starting multimodal EM refinement: {n_spots} spots, {n_genes} genes, {n_types} types")

    # Pass 1.5: Learn anchor genes using Spearman correlation
    logger.info("Pass 1.5: Selecting anchor genes (Spearman, adaptive threshold)...")
    anchors, weights = select_anchor_genes(
        GEX=GEX,
        Y_protein=Y_protein,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        n_anchors=n_anchors,
        min_correlation=min_correlation,
        min_anchors_per_type=min_anchors_per_type,
        min_expressing_spots=min_expressing_spots,
        sparse_aware=sparse_aware,
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
