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

            # Other cell types get 0 for this anchor gene
            # (already initialized to 0)

        else:
            # FREE: non-anchor gene can load on any cell type
            # Solve least squares: GEX[:, g] ≈ Y @ E[:, g]
            # E[:, g] = (Y^T Y)^{-1} Y^T GEX[:, g]
            try:
                E[:, g_idx], _, _, _ = np.linalg.lstsq(Y, gex_g, rcond=None)
                # Clip negative values (expression should be non-negative)
                E[:, g_idx] = np.clip(E[:, g_idx], 0, None)
            except np.linalg.LinAlgError:
                # Fallback: uniform assignment
                E[:, g_idx] = np.mean(gex_g) / n_types

    return E
