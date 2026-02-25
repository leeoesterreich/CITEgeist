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


def compute_module_aware_enrichment(
    spot_idx: int,
    neighborhood_expression: np.ndarray,
    base_enrichment: np.ndarray,
    anchor_genes: Dict[int, List[int]],
    module_weight: float = 0.5,
    min_neighbors_for_corr: int = 10,
) -> np.ndarray:
    """
    Compute module-aware enrichment by correlating genes with anchors in neighborhood.

    This function adjusts base enrichment scores based on how each gene correlates
    with anchor genes in the local spatial neighborhood. Genes that co-express with
    a cell type's anchor genes get boosted enrichment for that cell type, helping
    allocate weakly-expressed genes to the appropriate cell types.

    Args:
        spot_idx: Index of center spot (for logging)
        neighborhood_expression: (N_neighbors x M_genes) expression in neighborhood
        base_enrichment: (M_genes x T_celltypes) base enrichment scores
        anchor_genes: {cell_type_idx: [gene_indices]} anchor genes per type
        module_weight: Weight for module signal vs base enrichment (0-1)
        min_neighbors_for_corr: Minimum neighbors for reliable correlation

    Returns:
        adjusted_enrichment: (M_genes x T_celltypes) module-adjusted enrichment
    """
    M, T = base_enrichment.shape
    adjusted = base_enrichment.copy()

    # Skip if neighborhood too small for reliable correlation
    if neighborhood_expression.shape[0] < min_neighbors_for_corr:
        return adjusted

    for g in range(M):
        gene_expr = neighborhood_expression[:, g]

        # Skip if gene has no variance in neighborhood
        if np.std(gene_expr) < 1e-6:
            continue

        # Compute correlation with each cell type's anchor genes
        module_scores = np.zeros(T)
        for t in range(T):
            if t not in anchor_genes or not anchor_genes[t]:
                continue

            # Average correlation with this cell type's anchors
            anchor_corrs = []
            for anchor_idx in anchor_genes[t]:
                if anchor_idx >= neighborhood_expression.shape[1]:
                    continue
                anchor_expr = neighborhood_expression[:, anchor_idx]
                if np.std(anchor_expr) > 1e-6:
                    r, _ = pearsonr(gene_expr, anchor_expr)
                    if not np.isnan(r):
                        anchor_corrs.append(max(0, r))  # only positive correlations

            if anchor_corrs:
                module_scores[t] = np.mean(anchor_corrs)

        # Normalize module scores to sum to 1 (if any signal)
        if module_scores.sum() > 0:
            module_scores = module_scores / module_scores.sum()
        else:
            module_scores = np.ones(T) / T  # fallback to uniform

        # Combine base enrichment with module signal
        adjusted[g, :] = (
            (1 - module_weight) * base_enrichment[g, :] +
            module_weight * module_scores
        )

        # Re-normalize to sum to 1
        row_sum = adjusted[g, :].sum()
        if row_sum > 0:
            adjusted[g, :] = adjusted[g, :] / row_sum

    return adjusted


def compute_softmax_target(
    enrichment: np.ndarray,
    temperature: float = 0.3,
) -> np.ndarray:
    """
    Compute softmax target distribution from enrichment scores.

    Converts enrichment scores into a probability distribution using temperature-
    scaled softmax. Lower temperatures produce sharper distributions that concentrate
    on high-enrichment cell types; higher temperatures produce softer distributions.

    This is used to create target distributions for KL-divergence regularization,
    replacing uniform L2 regularization with enrichment-aware allocation targets.

    Args:
        enrichment: (T,) enrichment scores per cell type (can be any non-negative values)
        temperature: Softmax temperature parameter (default 0.3)
            - Lower values (0.1-0.3): Sharp, concentrate on top cell type
            - Higher values (0.5-1.0): Soft, more uniform distribution

    Returns:
        target: (T,) probability distribution summing to 1.0

    Example:
        >>> enrichment = np.array([0.4, 0.3, 0.1, 0.1, 0.1])
        >>> target = compute_softmax_target(enrichment, temperature=0.3)
        >>> print(target)  # [0.63, 0.28, 0.03, 0.03, 0.03] approximately
    """
    # Apply temperature scaling
    logits = enrichment / temperature

    # Numerical stability: subtract max before exp
    logits = logits - logits.max()

    # Softmax
    exp_logits = np.exp(logits)
    target = exp_logits / exp_logits.sum()

    # Clip to avoid zeros (for KL stability in downstream Gurobi objective)
    target = np.clip(target, 1e-6, 1.0)

    # Re-normalize after clipping
    target = target / target.sum()

    return target


def compute_kl_penalty_coefficients(
    target_distribution: np.ndarray,
    total_counts: int,
    lambda_kl: float = 0.1,
) -> Dict[str, np.ndarray]:
    """
    Compute coefficients for KL-divergence penalty in Gurobi objective.

    The KL penalty is approximated as a quadratic penalty that pulls allocations
    toward the target distribution: lambda * (X[j] - target[j])^2

    This differs from standard L2 regularization which uses uniform targets.
    By using enrichment-based targets from compute_softmax_target(), we guide
    gene allocations toward biologically plausible cell types.

    Args:
        target_distribution: (T,) target probabilities (must sum to 1)
        total_counts: Total counts to allocate for this gene
        lambda_kl: Penalty weight (default 0.1)
            - Higher values: Stronger pull toward target, less variance
            - Lower values: Weaker pull, more data-driven allocation

    Returns:
        Dict containing:
            - 'target_counts': (T,) target allocation in count space
            - 'penalty_weight': Normalized penalty weight for Gurobi

    Example:
        >>> target = np.array([0.5, 0.3, 0.2])
        >>> coeffs = compute_kl_penalty_coefficients(target, total_counts=100, lambda_kl=0.1)
        >>> print(coeffs['target_counts'])  # [50, 30, 20]
    """
    # Convert probability distribution to count space
    target_counts = target_distribution * total_counts

    # Normalize penalty by total_counts to make lambda_kl scale-invariant
    # This ensures consistent regularization behavior across genes with different
    # expression levels (highly expressed genes don't dominate the penalty)
    penalty_weight = lambda_kl / (total_counts + 1)

    return {
        'target_counts': target_counts,
        'penalty_weight': penalty_weight,
    }
