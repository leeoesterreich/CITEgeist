"""
GEX module-aware enrichment and KL regularization for gene expression deconvolution.

This module implements:
1. Anchor gene discovery - find genes correlated with cell type proportions
2. Module-aware enrichment - boost gene enrichment based on anchor correlation
3. Softmax KL allocation - replace L2 regularization with KL-divergence
"""

import logging
from typing import Dict, List, Tuple

import numpy as np
from scipy.stats import pearsonr

logger = logging.getLogger(__name__)


def discover_anchor_genes(
    gene_expression: np.ndarray,
    cell_proportions: np.ndarray,
    min_anchors: int = 5,
    max_anchors: int = 10,
    initial_min_correlation: float = 0.3,
    min_expressing_spots: float = 0.1,
) -> Tuple[Dict[int, List[int]], Dict[int, float]]:
    """
    Discover anchor genes for each cell type based on correlation with proportions.

    Anchor genes are genes whose spatial expression pattern correlates strongly
    with cell type proportions. These serve as reference markers for module-aware
    enrichment, helping allocate weakly-expressed genes to the appropriate cell types.

    The algorithm uses adaptive thresholding: starting at a high correlation threshold
    (0.30), it progressively lowers the threshold until at least min_anchors are found
    for each cell type. This ensures even cell types with weak marker signals get
    reasonable anchor sets.

    Args:
        gene_expression: (N_spots x M_genes) expression matrix (CPM or similar)
        cell_proportions: (N_spots x T_celltypes) proportion matrix (rows sum to 1)
        min_anchors: Minimum anchors per cell type (threshold adapts to reach this)
        max_anchors: Maximum anchors per cell type
        initial_min_correlation: Starting correlation threshold (default 0.3)
        min_expressing_spots: Minimum fraction of spots where gene must be expressed
            (filters out sparse genes that may produce spurious correlations)

    Returns:
        anchors: Dict mapping cell type index to list of gene indices
            {cell_type_idx: [gene_indices]} sorted by correlation descending
        thresholds_used: Dict mapping cell type index to the correlation threshold
            actually used to select anchors for that type
            {cell_type_idx: correlation_threshold}

    Example:
        >>> proportions = model.cell_proportions  # from Module 3 Pass 1
        >>> gene_expr = adata.X  # normalized expression
        >>> anchors, thresholds = discover_anchor_genes(gene_expr, proportions)
        >>> print(f"Type 0 anchors: {anchors[0][:5]}")  # top 5 anchor genes
    """
    N, M = gene_expression.shape
    T = cell_proportions.shape[1]

    # Validate inputs
    if cell_proportions.shape[0] != N:
        raise ValueError(
            f"Mismatch: gene_expression has {N} spots, "
            f"cell_proportions has {cell_proportions.shape[0]}"
        )

    # Threshold sequence: 0.30 -> 0.25 -> 0.20 -> 0.15 -> 0.10
    # Progressively lower to find anchors for weak cell types
    threshold_sequence = [0.30, 0.25, 0.20, 0.15, 0.10]

    anchors = {}
    thresholds_used = {}

    for t in range(T):
        prop_vector = cell_proportions[:, t]

        # Skip if proportions have no variance (cell type absent)
        if np.std(prop_vector) < 1e-10:
            anchors[t] = []
            thresholds_used[t] = threshold_sequence[-1]
            logger.warning(f"Cell type {t}: zero variance in proportions, no anchors")
            continue

        # Compute correlations for all genes that pass filters
        all_correlations = []
        for g in range(M):
            gene_vector = gene_expression[:, g]

            # Filter: gene must be expressed in enough spots
            expressing_fraction = (gene_vector > 0).mean()
            if expressing_fraction < min_expressing_spots:
                continue

            # Filter: need variance for correlation
            if np.std(gene_vector) < 1e-10:
                continue

            # Compute Pearson correlation
            r, p = pearsonr(gene_vector, prop_vector)

            # Only keep significant positive correlations
            if p < 0.05 and not np.isnan(r) and r > 0:
                all_correlations.append((g, r))

        # Sort by correlation descending
        all_correlations.sort(key=lambda x: -x[1])

        # Step through thresholds until we have min_anchors
        selected_threshold = threshold_sequence[-1]  # default to floor
        for threshold in threshold_sequence:
            candidates = [g for g, r in all_correlations if r >= threshold]
            if len(candidates) >= min_anchors:
                selected_threshold = threshold
                break

        # Select anchors at chosen threshold (capped at max)
        final_candidates = [g for g, r in all_correlations if r >= selected_threshold]
        anchors[t] = final_candidates[:max_anchors]
        thresholds_used[t] = selected_threshold

        # Log info about anchor discovery
        n_anchors = len(anchors[t])
        if n_anchors < min_anchors:
            logger.warning(
                f"Cell type {t}: only {n_anchors} anchors found "
                f"(min={min_anchors}, threshold={selected_threshold:.2f})"
            )
        else:
            logger.debug(
                f"Cell type {t}: {n_anchors} anchors at threshold={selected_threshold:.2f}"
            )

    return anchors, thresholds_used
