"""
Unified spatial colocalization that works for both spot-level and cell-level data.

The key insight: markers that define the same cell type have SIMILAR SPATIAL PATTERNS
regardless of whether we're measuring spots (mixed cells) or individual cells.

This module provides resolution-agnostic colocalization metrics:
1. Spatial pattern correlation: Do markers have similar spatial distributions?
2. Neighborhood-aware similarity: Do markers co-occur in spatial neighborhoods?
3. Permutation-based significance: Is the observed pattern greater than chance?

Works identically for:
- Spot-level data (Visium): Each spot has multiple cells, markers from same cell type colocalize
- Cell-level data (Visium HD, Xenium): Each cell is one type, markers in same cells cluster spatially
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from joblib import Parallel, delayed

logger = logging.getLogger(__name__)


@dataclass
class UnifiedColocalizationResult:
    """Result of unified spatial colocalization analysis."""

    marker_a: str
    marker_b: str

    # Primary metric: spatial pattern correlation
    spatial_pattern_correlation: float  # How similar are their spatial distributions
    spatial_pattern_pvalue: float  # Permutation p-value

    # Secondary metrics
    neighborhood_cosimilarity: float  # Neighborhood-aggregated pattern similarity
    expression_correlation: float  # Direct expression correlation (Spearman)

    # Co-expression metrics (for cell type separation)
    coexpression_jaccard: float  # Jaccard of high-expressing locations
    ratio_consistency: float  # Consistency of A/B ratio across locations
    neighborhood_coexpression: float  # Neighborhood-based co-expression (handles dropout)
    coexpression_score: float  # Combined co-expression score

    # Combined score (weighted)
    unified_score: float


def _build_neighborhood_matrix(
    coords: NDArray[np.floating],
    k: int = 15,
) -> NDArray[np.floating]:
    """
    Build row-normalized neighborhood weight matrix.

    Args:
        coords: Spatial coordinates (n_spots, 2).
        k: Number of neighbors per location.

    Returns:
        Row-normalized weight matrix W where W[i,j] = 1/k if j is neighbor of i.
    """
    n = coords.shape[0]
    tree = cKDTree(coords)

    # Query k+1 to include self, then exclude
    _, indices = tree.query(coords, k=min(k + 1, n))

    # Build weight matrix
    W = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        neighbors = indices[i, 1:] if indices.ndim > 1 else []
        if len(neighbors) > 0:
            W[i, neighbors] = 1.0 / len(neighbors)

    return W


def _compute_coexpression_jaccard(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    threshold_percentile: float = 75.0,
) -> float:
    """
    Compute Jaccard index of high-expressing locations.

    For single-cell: measures what fraction of expressing cells express BOTH markers.
    For spot-level: measures what fraction of high spots are high for BOTH markers.

    Args:
        values_a: Expression for marker A.
        values_b: Expression for marker B.
        threshold_percentile: Percentile threshold for "high" expression.

    Returns:
        Jaccard index [0, 1].
    """
    # Get thresholds from non-zero values
    nonzero_a = values_a[values_a > 0]
    nonzero_b = values_b[values_b > 0]

    if len(nonzero_a) < 10 or len(nonzero_b) < 10:
        return 0.0

    thresh_a = np.percentile(nonzero_a, threshold_percentile)
    thresh_b = np.percentile(nonzero_b, threshold_percentile)

    high_a = values_a > thresh_a
    high_b = values_b > thresh_b

    intersection = (high_a & high_b).sum()
    union = (high_a | high_b).sum()

    return intersection / union if union > 0 else 0.0


def _compute_ratio_consistency(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    min_coexpress: int = 10,
) -> float:
    """
    Compute ratio consistency in co-expressing locations.

    Same cell type markers have CONSTANT ratio across locations (both track
    the same underlying cell type abundance).

    Different cell type markers have VARIABLE ratio (different sources).

    Args:
        values_a: Expression for marker A.
        values_b: Expression for marker B.
        min_coexpress: Minimum co-expressing locations required.

    Returns:
        Ratio consistency score [0, 1]. High = consistent ratio = same cell type.
    """
    # Find locations where both are expressed
    coexpress = (values_a > 0) & (values_b > 0)

    if coexpress.sum() < min_coexpress:
        # Too few co-expressing locations - likely different cell types
        return 0.0

    # Compute log ratio (handles different expression scales)
    log_ratio = np.log1p(values_a[coexpress]) - np.log1p(values_b[coexpress])

    # Standard deviation of log ratio
    # Low std = consistent ratio = same cell type
    # High std = variable ratio = different cell types
    std_ratio = np.std(log_ratio)

    # Convert to [0, 1] score where 1 = perfectly consistent
    # Use sigmoid-like transformation
    consistency = 1.0 / (1.0 + std_ratio)

    return float(consistency)


def _compute_neighborhood_coexpression(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    W: NDArray[np.floating],
) -> float:
    """
    Compute neighborhood-based co-expression for sparse single-cell data.

    Instead of asking "do A and B appear in the SAME cell?" (fails with dropout),
    ask "do cells expressing A have neighbors expressing B?"

    This works because:
    - Same cell type: T cells cluster → CD3E+ cells have CD4+ neighbors
    - Different cell types: T cells separate from B cells → CD3E+ cells don't have CD19+ neighbors

    Args:
        values_a: Expression for marker A.
        values_b: Expression for marker B.
        W: Row-normalized neighborhood weight matrix.

    Returns:
        Neighborhood co-expression score [0, 1].
    """
    # Cells expressing each marker (any expression)
    cells_a = values_a > 0
    cells_b = values_b > 0

    n_a = cells_a.sum()
    n_b = cells_b.sum()

    if n_a < 5 or n_b < 5:
        return 0.0

    # For each cell, compute fraction of neighbors expressing the other marker
    neighbor_b_frac = W @ cells_b.astype(float)  # Fraction of neighbors expressing B
    neighbor_a_frac = W @ cells_a.astype(float)  # Fraction of neighbors expressing A

    # Score 1: Among A-expressing cells, what's the avg fraction of B-expressing neighbors?
    a_to_b = np.mean(neighbor_b_frac[cells_a]) if n_a > 0 else 0.0

    # Score 2: Among B-expressing cells, what's the avg fraction of A-expressing neighbors?
    b_to_a = np.mean(neighbor_a_frac[cells_b]) if n_b > 0 else 0.0

    # Symmetric score (geometric mean)
    if a_to_b > 0 and b_to_a > 0:
        return float(np.sqrt(a_to_b * b_to_a))
    else:
        return 0.0


def _compute_unified_coexpression(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    W: Optional[NDArray[np.floating]] = None,
    threshold_percentile: float = 75.0,
) -> Tuple[float, float, float, float]:
    """
    Compute unified co-expression score that works for both resolutions.

    Combines:
    1. Jaccard of high-expressing locations
    2. Ratio consistency (for spot-level separation)
    3. Neighborhood co-expression (for sparse single-cell data)

    Args:
        values_a: Expression for marker A.
        values_b: Expression for marker B.
        W: Optional neighborhood weight matrix (required for neighborhood score).
        threshold_percentile: Percentile for Jaccard threshold.

    Returns:
        Tuple of (jaccard, ratio_consistency, neighborhood_coexpr, combined_score).
    """
    jaccard = _compute_coexpression_jaccard(values_a, values_b, threshold_percentile)
    ratio_consistency = _compute_ratio_consistency(values_a, values_b)

    # Neighborhood co-expression (handles sparse single-cell dropout)
    if W is not None:
        neighborhood_coexpr = _compute_neighborhood_coexpression(values_a, values_b, W)
    else:
        neighborhood_coexpr = 0.0

    # Combined score strategy:
    # - If we have good within-cell co-expression (jaccard > 0.1), use jaccard + ratio
    # - If sparse (jaccard < 0.1), rely on neighborhood co-expression
    if jaccard > 0.1 and ratio_consistency > 0:
        # Dense data: use traditional approach
        combined = np.sqrt(jaccard * ratio_consistency)
    elif neighborhood_coexpr > 0:
        # Sparse data: use neighborhood approach
        combined = neighborhood_coexpr
    else:
        combined = 0.0

    return jaccard, ratio_consistency, neighborhood_coexpr, combined


def _compute_spatial_fingerprint(
    values: NDArray[np.floating],
    W: NDArray[np.floating],
    self_weight: float = 0.5,
) -> NDArray[np.floating]:
    """
    Compute spatial fingerprint: weighted combination of local and neighborhood expression.

    fingerprint_i = self_weight * value_i + (1 - self_weight) * neighbor_avg_i

    This captures both:
    - Direct expression at each location
    - Expression context from spatial neighborhood

    Args:
        values: Expression vector (n_locations,).
        W: Row-normalized neighborhood weight matrix.
        self_weight: Weight for self vs neighborhood (default: 0.5).

    Returns:
        Spatial fingerprint vector (n_locations,).
    """
    # Neighborhood average
    neighbor_avg = W @ values

    # Weighted combination
    fingerprint = self_weight * values + (1 - self_weight) * neighbor_avg

    return fingerprint


def _compute_spatial_pattern_correlation(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    W: NDArray[np.floating],
    self_weight: float = 0.5,
) -> float:
    """
    Compute correlation between spatial fingerprints of two markers.

    This measures: "Do markers A and B have similar spatial distributions?"

    Works for both spot and cell data because:
    - Spot data: Markers from same cell type have similar spatial patterns
      (high where that cell type is abundant)
    - Cell data: Markers from same cell type have similar spatial patterns
      (high where those cells are located, including neighborhood context)

    Args:
        values_a: Expression for marker A (n_locations,).
        values_b: Expression for marker B (n_locations,).
        W: Row-normalized neighborhood weight matrix.
        self_weight: Weight for self vs neighborhood in fingerprint.

    Returns:
        Spearman correlation of spatial fingerprints.
    """
    # Compute spatial fingerprints
    fp_a = _compute_spatial_fingerprint(values_a, W, self_weight)
    fp_b = _compute_spatial_fingerprint(values_b, W, self_weight)

    # Handle constant vectors
    if np.std(fp_a) < 1e-10 or np.std(fp_b) < 1e-10:
        return 0.0

    # Spearman correlation (robust to outliers)
    corr, _ = spearmanr(fp_a, fp_b)

    return float(corr) if not np.isnan(corr) else 0.0


def _compute_neighborhood_cosimilarity(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    W: NDArray[np.floating],
) -> float:
    """
    Compute neighborhood co-similarity.

    For each location, compute the neighborhood average of both markers,
    then correlate these neighborhood profiles.

    This captures: "Are neighborhoods that are high in A also high in B?"

    This is related to bivariate Moran's I but uses both markers'
    neighborhood context rather than one lagged against the other.

    Args:
        values_a: Expression for marker A.
        values_b: Expression for marker B.
        W: Row-normalized neighborhood weight matrix.

    Returns:
        Correlation of neighborhood profiles.
    """
    # Neighborhood averages
    neigh_a = W @ values_a
    neigh_b = W @ values_b

    # Handle constant
    if np.std(neigh_a) < 1e-10 or np.std(neigh_b) < 1e-10:
        return 0.0

    corr, _ = spearmanr(neigh_a, neigh_b)

    return float(corr) if not np.isnan(corr) else 0.0


def _permutation_test_spatial_pattern(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    W: NDArray[np.floating],
    observed_corr: float,
    n_perm: int = 499,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """
    Permutation test for spatial pattern correlation significance.

    Null hypothesis: No spatial relationship between markers.
    Test: Shuffle spatial positions of one marker and recompute correlation.

    Args:
        values_a: Expression for marker A.
        values_b: Expression for marker B.
        W: Neighborhood weight matrix.
        observed_corr: Observed spatial pattern correlation.
        n_perm: Number of permutations.
        rng: Random generator.

    Returns:
        P-value (fraction of null >= observed).
    """
    if rng is None:
        rng = np.random.default_rng(42)

    null_corrs = np.zeros(n_perm)

    for i in range(n_perm):
        # Shuffle B's values across locations
        perm_b = rng.permutation(values_b)
        null_corrs[i] = _compute_spatial_pattern_correlation(values_a, perm_b, W)

    # Two-sided test (both positive and negative associations)
    p_value = (np.sum(np.abs(null_corrs) >= np.abs(observed_corr)) + 1) / (n_perm + 1)

    return float(p_value)


def _process_pair_unified(
    i: int,
    j: int,
    marker_names: List[str],
    X: NDArray[np.floating],
    W: NDArray[np.floating],
    n_perm: int,
    seed: int,
    scoring_method: str = "geometric",
) -> UnifiedColocalizationResult:
    """Process a single marker pair with unified colocalization."""

    pair_seed = seed + i * 1000 + j
    rng = np.random.default_rng(pair_seed)

    name_a = marker_names[i]
    name_b = marker_names[j]

    values_a = X[:, i]
    values_b = X[:, j]

    # Primary: Spatial pattern correlation
    spatial_corr = _compute_spatial_pattern_correlation(values_a, values_b, W)
    spatial_pval = _permutation_test_spatial_pattern(
        values_a, values_b, W, spatial_corr, n_perm=n_perm, rng=rng
    )

    # Secondary: Neighborhood co-similarity
    neigh_cosim = _compute_neighborhood_cosimilarity(values_a, values_b, W)

    # Direct expression correlation
    if np.std(values_a) > 1e-10 and np.std(values_b) > 1e-10:
        expr_corr, _ = spearmanr(values_a, values_b)
        expr_corr = float(expr_corr) if not np.isnan(expr_corr) else 0.0
    else:
        expr_corr = 0.0

    # Co-expression metrics (CRITICAL for cell type separation at all resolutions)
    coexpr_jaccard, ratio_consist, neigh_coexpr, coexpr_score = _compute_unified_coexpression(
        values_a, values_b, W=W
    )

    # Normalize correlations to [0, 1] range
    spatial_norm = (spatial_corr + 1) / 2
    expr_norm = (expr_corr + 1) / 2

    # Scoring methods
    if scoring_method == "neighborhood":
        # NEW: Use neighborhood co-expression (handles sparse single-cell dropout)
        # Do cells expressing A have neighbors expressing B?
        unified_score = neigh_coexpr
    elif scoring_method == "ratio":
        # Use Jaccard + Ratio Consistency (works for both resolutions)
        # Must have high co-expression Jaccard AND consistent ratio
        unified_score = coexpr_score
    elif scoring_method == "spatial_ratio":
        # Combine spatial pattern with ratio consistency
        # Must be spatially similar AND have consistent ratio
        unified_score = np.sqrt(spatial_norm * coexpr_score) if coexpr_score > 0 else 0.0
    elif scoring_method == "spatial_neighborhood":
        # Combine spatial pattern with neighborhood co-expression
        # Best for sparse single-cell data
        unified_score = np.sqrt(spatial_norm * neigh_coexpr) if neigh_coexpr > 0 else 0.0
    elif scoring_method == "geometric":
        # Geometric mean of spatial and expression correlation
        unified_score = np.sqrt(spatial_norm * expr_norm)
    elif scoring_method == "min":
        # Bottleneck: limited by weaker signal
        unified_score = min(spatial_norm, expr_norm)
    elif scoring_method == "harmonic":
        # Harmonic mean: strongly penalizes low values
        if spatial_norm > 0 and expr_norm > 0:
            unified_score = 2 * spatial_norm * expr_norm / (spatial_norm + expr_norm)
        else:
            unified_score = 0.0
    else:
        # Legacy weighted average (for comparison)
        unified_score = 0.5 * spatial_norm + 0.3 * (neigh_cosim + 1) / 2 + 0.2 * expr_norm

    return UnifiedColocalizationResult(
        marker_a=name_a,
        marker_b=name_b,
        spatial_pattern_correlation=spatial_corr,
        spatial_pattern_pvalue=spatial_pval,
        neighborhood_cosimilarity=neigh_cosim,
        expression_correlation=expr_corr,
        coexpression_jaccard=coexpr_jaccard,
        ratio_consistency=ratio_consist,
        neighborhood_coexpression=neigh_coexpr,
        coexpression_score=coexpr_score,
        unified_score=unified_score,
    )


def analyze_colocalization_unified(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    *,
    neighbor_k: int = 15,
    n_permutations: int = 499,
    n_jobs: int = -1,
    seed: int = 42,
    verbose: bool = True,
    scoring_method: str = "neighborhood",
) -> List[UnifiedColocalizationResult]:
    """
    Analyze pairwise spatial colocalization using unified resolution-agnostic approach.

    This method works identically for spot-level and cell-level data by measuring
    spatial pattern similarity rather than within-location co-occurrence.

    NOTE: Default scoring_method="neighborhood" was selected based on cross-resolution
    benchmarking showing best average Jaccard (0.384) across single-cell and spot-level
    data. See tests/test_unified_both_resolutions.py for comparison.

    Key metrics:
    1. Spatial pattern correlation: Do markers have similar spatial distributions?
       - Computed on spatial fingerprints (local + neighborhood context)
       - Tested via permutation (shuffle spatial positions)

    2. Neighborhood co-similarity: Are neighborhoods high in A also high in B?

    3. Expression correlation: Standard Spearman correlation.

    Args:
        X: Expression matrix (n_locations, n_markers).
        coords: Spatial coordinates (n_locations, 2).
        marker_names: Names for each marker.
        neighbor_k: Number of neighbors for spatial context (default: 15).
        n_permutations: Permutations per pair for significance (default: 499).
        n_jobs: Parallel jobs (-1 = all cores).
        seed: Random seed.
        verbose: Log progress.
        scoring_method: How to combine spatial and expression metrics:
            - "neighborhood" (RECOMMENDED): Uses neighborhood co-expression.
              Best cross-resolution performance (Avg Jaccard=0.384).
              Dramatically better on spot-level (+0.302 vs weighted).
            - "spatial_neighborhood": sqrt(spatial * neighborhood co-expression)
            - "ratio": Jaccard × ratio consistency (for single-cell)
            - "geometric": sqrt(spatial * expr)
            - "weighted": legacy (0.5*spatial + 0.3*neigh + 0.2*expr)

    Returns:
        List of UnifiedColocalizationResult for each marker pair.
    """
    n_markers = len(marker_names)
    n_pairs = n_markers * (n_markers - 1) // 2
    n_locations = X.shape[0]

    if verbose:
        logger.info(
            f"Unified colocalization: {n_markers} markers ({n_pairs} pairs) "
            f"across {n_locations} locations, scoring={scoring_method}"
        )
        logger.info(f"Building {neighbor_k}-NN neighborhood graph...")

    # Build neighborhood weight matrix
    W = _build_neighborhood_matrix(coords, k=neighbor_k)

    if verbose:
        logger.info(f"Running parallel analysis with {n_jobs} workers...")

    # Generate all pairs
    pairs_ij = [(i, j) for i in range(n_markers) for j in range(i + 1, n_markers)]

    # Parallel processing
    results = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_process_pair_unified)(
            i, j, marker_names, X, W, n_permutations, seed, scoring_method
        )
        for i, j in pairs_ij
    )

    if verbose:
        # Summary statistics
        scores = [r.unified_score for r in results]
        pvals = [r.spatial_pattern_pvalue for r in results]
        sig_count = sum(1 for p in pvals if p < 0.05)

        logger.info(
            f"Complete: {sig_count}/{n_pairs} pairs with p < 0.05, "
            f"score range [{min(scores):.3f}, {max(scores):.3f}]"
        )

    return results


def build_colocalization_graph(
    results: List[UnifiedColocalizationResult],
    marker_names: List[str],
    *,
    pvalue_threshold: float = 0.05,
    min_score: Optional[float] = None,
    top_k: int = 5,
    verbose: bool = True,
) -> List[Tuple[str, str, float]]:
    """
    Build filtered colocalization graph from unified results.

    Filtering steps:
    1. P-value threshold on spatial pattern correlation
    2. Optional minimum score threshold
    3. Mutual top-k sparsification

    Args:
        results: Output from analyze_colocalization_unified().
        marker_names: All marker names.
        pvalue_threshold: Maximum p-value to keep edge (default: 0.05).
        min_score: Minimum unified score (optional).
        top_k: Keep only mutual top-k partners per marker.
        verbose: Log progress.

    Returns:
        List of (marker_a, marker_b, score) edges for graph construction.
    """
    # Step 1: P-value filter
    filtered = [r for r in results if r.spatial_pattern_pvalue < pvalue_threshold]

    if verbose:
        logger.info(
            f"P-value filter (p < {pvalue_threshold}): "
            f"{len(filtered)}/{len(results)} pairs"
        )

    # Step 2: Optional score filter
    if min_score is not None:
        filtered = [r for r in filtered if r.unified_score >= min_score]
        if verbose:
            logger.info(f"Score filter (>= {min_score}): {len(filtered)} pairs")

    # Step 3: Mutual top-k
    if top_k > 0 and len(filtered) > 0:
        # Build per-marker rankings
        marker_to_pairs: Dict[str, List[Tuple[str, float]]] = {m: [] for m in marker_names}

        for r in filtered:
            marker_to_pairs[r.marker_a].append((r.marker_b, r.unified_score))
            marker_to_pairs[r.marker_b].append((r.marker_a, r.unified_score))

        # Get top-k for each marker
        marker_topk: Dict[str, set] = {}
        for marker, pairs in marker_to_pairs.items():
            sorted_pairs = sorted(pairs, key=lambda x: x[1], reverse=True)
            marker_topk[marker] = set(p[0] for p in sorted_pairs[:top_k])

        # Keep only mutual top-k edges
        mutual_filtered = []
        for r in filtered:
            a, b = r.marker_a, r.marker_b
            if b in marker_topk.get(a, set()) and a in marker_topk.get(b, set()):
                mutual_filtered.append(r)

        filtered = mutual_filtered
        if verbose:
            logger.info(f"Mutual top-{top_k}: {len(filtered)} pairs")

    # Build edge list
    edges = [
        (r.marker_a, r.marker_b, r.unified_score)
        for r in filtered
    ]

    if verbose:
        logger.info(f"Final graph: {len(edges)} edges")

    return edges


def run_unified_profile_discovery(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    *,
    neighbor_k: int = 15,
    n_permutations: int = 499,
    pvalue_threshold: float = 0.05,
    top_k: int = 5,
    n_jobs: int = -1,
    seed: int = 42,
    verbose: bool = True,
) -> Tuple[List[UnifiedColocalizationResult], List[Tuple[str, str, float]]]:
    """
    Complete unified colocalization pipeline.

    Runs:
    1. Unified colocalization analysis (resolution-agnostic)
    2. Graph construction with filtering

    The resulting edges can be fed to the existing discover_profiles() function
    for hierarchical clustering and profile extraction.

    Args:
        X: Expression matrix (n_locations, n_markers).
        coords: Spatial coordinates (n_locations, 2).
        marker_names: Names for each marker.
        neighbor_k: Neighbors for spatial context.
        n_permutations: Permutations for significance.
        pvalue_threshold: P-value cutoff for edges.
        top_k: Mutual top-k filter.
        n_jobs: Parallel workers.
        seed: Random seed.
        verbose: Log progress.

    Returns:
        Tuple of (all_results, filtered_edges).
    """
    # Run analysis
    results = analyze_colocalization_unified(
        X, coords, marker_names,
        neighbor_k=neighbor_k,
        n_permutations=n_permutations,
        n_jobs=n_jobs,
        seed=seed,
        verbose=verbose,
    )

    # Build graph
    edges = build_colocalization_graph(
        results, marker_names,
        pvalue_threshold=pvalue_threshold,
        top_k=top_k,
        verbose=verbose,
    )

    return results, edges
