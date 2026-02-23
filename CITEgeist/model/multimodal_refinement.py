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
