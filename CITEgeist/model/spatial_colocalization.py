"""
Spatial colocalization analysis for marker pairs and profile discovery.

Analyzes three types of spatial relationships:
1. Same-spot co-occurrence: Do markers appear together in the same spots?
2. Expression correlation: Are marker intensities correlated?
3. Adjacent-spot enrichment: Are markers enriched in neighboring spots?

Then discovers cell type profiles by:
1. Building a significance-filtered graph from colocalization p-values
2. Finding connected components (separate lineages)
3. Hierarchical clustering within each component
4. Dynamic tree cutting to extract profiles automatically
"""

import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.cluster.hierarchy import linkage, fcluster, inconsistent
from scipy.spatial import cKDTree
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr, spearmanr
from joblib import Parallel, delayed

# Module logger
logger = logging.getLogger(__name__)


@dataclass
class MarkerPairColocalization:
    """Colocalization metrics for a pair of markers."""

    marker_a: str
    marker_b: str

    # Same-spot co-occurrence (binary-based, kept for backwards compatibility)
    jaccard_index: float  # |A & B| / |A | B|
    co_occurrence_spots: int  # Number of spots with both markers
    co_occurrence_fraction: float  # Fraction of spots with both

    # Expression correlation (continuous)
    pearson_r: float  # Pearson correlation of intensities
    spearman_rho: float  # Spearman rank correlation
    correlation_pvalue: float  # P-value for correlation test

    # Continuous similarity metrics (binarization-free)
    cosine_similarity: float  # Scale-invariant pattern similarity
    bivariate_morans_i: float  # Spatial cross-correlation (aggregated if multi-scale)
    bivariate_morans_pvalue: float  # P-value from permutation test

    # Neighborhood analysis (binary-based, kept for backwards compatibility)
    neighbor_enrichment_ab: float  # Enrichment of B in neighbors of A-positive spots
    neighbor_enrichment_ba: float  # Enrichment of A in neighbors of B-positive spots
    neighbor_enrichment_pvalue: float
    mutual_neighbor_enrichment: float  # Mean of A->B and B->A enrichment

    # Combined score
    colocalization_score: float  # Weighted combination of all metrics

    # Multi-scale analysis (optional, for detailed analysis)
    bivariate_morans_per_scale: Optional[Dict[int, float]] = None  # k -> I value
    bivariate_morans_pvalue_per_scale: Optional[Dict[int, float]] = None  # k -> p-value
    bivariate_morans_best_scale: Optional[int] = None  # Scale with highest signal


@dataclass
class ColocalizationResult:
    """Results container for colocalization analysis."""

    pairs: List[MarkerPairColocalization]
    marker_names: List[str]
    n_spots: int
    neighbor_k: int
    marker_specificity: Dict[str, float] = field(default_factory=dict)  # Gini-based specificity scores

    def to_dataframe(self) -> pd.DataFrame:
        """Convert to ranked DataFrame sorted by colocalization_score descending."""
        records = [
            {
                "marker_a": p.marker_a,
                "marker_b": p.marker_b,
                "colocalization_score": p.colocalization_score,
                "jaccard_index": p.jaccard_index,
                "co_occurrence_spots": p.co_occurrence_spots,
                "co_occurrence_fraction": p.co_occurrence_fraction,
                "pearson_r": p.pearson_r,
                "spearman_rho": p.spearman_rho,
                "correlation_pvalue": p.correlation_pvalue,
                "neighbor_enrichment_ab": p.neighbor_enrichment_ab,
                "neighbor_enrichment_ba": p.neighbor_enrichment_ba,
                "neighbor_enrichment_pvalue": p.neighbor_enrichment_pvalue,
                "mutual_neighbor_enrichment": p.mutual_neighbor_enrichment,
            }
            for p in self.pairs
        ]
        df = pd.DataFrame(records)
        return df.sort_values("colocalization_score", ascending=False).reset_index(drop=True)

    def get_pairs_for_marker(self, marker: str) -> List[MarkerPairColocalization]:
        """Get all pairs involving a specific marker."""
        return [p for p in self.pairs if p.marker_a == marker or p.marker_b == marker]

    def top_pairs(self, n: int = 20) -> List[MarkerPairColocalization]:
        """Return top N pairs by colocalization score."""
        sorted_pairs = sorted(self.pairs, key=lambda x: x.colocalization_score, reverse=True)
        return sorted_pairs[:n]


def _build_neighbor_graph(
    coords: NDArray[np.floating],
    k: int,
) -> List[List[int]]:
    """
    Build k-NN neighbor list for each spot.

    Args:
        coords: Spatial coordinates (n_spots, 2).
        k: Number of neighbors per spot.

    Returns:
        List of neighbor indices for each spot.
    """
    n_spots = coords.shape[0]
    tree = cKDTree(coords)

    query_k = min(k + 1, n_spots)
    _, indices = tree.query(coords, k=query_k)

    # Handle 1D case
    if indices.ndim == 1:
        return [[] for _ in range(n_spots)]

    # Remove self (first index)
    neighbors = [list(idx[1:]) for idx in indices]
    return neighbors


def _binarize_markers(
    X: NDArray[np.floating],
    threshold_percentile: float,
) -> NDArray[np.bool_]:
    """
    Convert continuous expression to binary (positive/negative).

    Uses per-marker percentile threshold.

    Args:
        X: Expression matrix (n_spots, n_markers).
        threshold_percentile: Percentile to use as cutoff (e.g., 75).

    Returns:
        Binary matrix (n_spots, n_markers) indicating positive spots.
    """
    thresholds = np.percentile(X, threshold_percentile, axis=0)
    return X > thresholds


def _compute_jaccard(
    binary_a: NDArray[np.bool_],
    binary_b: NDArray[np.bool_],
) -> Tuple[float, int, float]:
    """
    Compute Jaccard index and co-occurrence statistics.

    Jaccard index: |A & B| / |A | B|

    Args:
        binary_a: Binary indicator for marker A (n_spots,).
        binary_b: Binary indicator for marker B (n_spots,).

    Returns:
        Tuple of (jaccard_index, co_occurrence_spots, co_occurrence_fraction).
    """
    intersection = np.sum(binary_a & binary_b)
    union = np.sum(binary_a | binary_b)
    n_spots = len(binary_a)

    jaccard = intersection / union if union > 0 else 0.0
    co_occurrence_fraction = intersection / n_spots if n_spots > 0 else 0.0

    return float(jaccard), int(intersection), float(co_occurrence_fraction)


def _compute_correlation(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
) -> Tuple[float, float, float]:
    """
    Compute Pearson and Spearman correlations.

    Args:
        values_a: Expression values for marker A (n_spots,).
        values_b: Expression values for marker B (n_spots,).

    Returns:
        Tuple of (pearson_r, spearman_rho, p_value).
    """
    # Handle constant arrays
    if np.std(values_a) < 1e-10 or np.std(values_b) < 1e-10:
        return 0.0, 0.0, 1.0

    try:
        pearson_r, p_pearson = pearsonr(values_a, values_b)
        spearman_rho, p_spearman = spearmanr(values_a, values_b)
        # Use the more conservative p-value
        p_value = max(p_pearson, p_spearman)
        return float(pearson_r), float(spearman_rho), float(p_value)
    except Exception:
        return 0.0, 0.0, 1.0


def _compute_cosine_similarity(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
) -> float:
    """
    Compute cosine similarity between two marker expression vectors.

    Cosine similarity is scale-invariant and measures whether markers
    have similar spatial patterns regardless of absolute magnitude.

    cos(A, B) = (A · B) / (||A|| ||B||)

    Args:
        values_a: Expression values for marker A (n_spots,).
        values_b: Expression values for marker B (n_spots,).

    Returns:
        Cosine similarity in range [-1, 1].
    """
    norm_a = np.linalg.norm(values_a)
    norm_b = np.linalg.norm(values_b)

    if norm_a < 1e-10 or norm_b < 1e-10:
        return 0.0

    return float(np.dot(values_a, values_b) / (norm_a * norm_b))


def _compute_bivariate_morans_i(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    neighbors: List[List[int]],
) -> float:
    """
    Compute bivariate Moran's I (spatial cross-correlation).

    This measures whether high values of A co-occur spatially with high values
    of B in neighboring spots. Answers: "Do markers peak in the same areas?"

    I_AB = Σ_i Σ_j w_ij * (A_i - μ_A) * (B_j - μ_B) / (n * σ_A * σ_B)

    Where w_ij = 1 if j is a neighbor of i, 0 otherwise.

    Args:
        values_a: Expression values for marker A (n_spots,).
        values_b: Expression values for marker B (n_spots,).
        neighbors: List of neighbor indices per spot.

    Returns:
        Bivariate Moran's I in range approximately [-1, 1].
    """
    n_spots = len(values_a)

    # Center the values
    mean_a = np.mean(values_a)
    mean_b = np.mean(values_b)
    std_a = np.std(values_a)
    std_b = np.std(values_b)

    if std_a < 1e-10 or std_b < 1e-10:
        return 0.0

    a_centered = (values_a - mean_a) / std_a
    b_centered = (values_b - mean_b) / std_b

    # Compute spatially-lagged B (mean of B in neighborhood of each spot)
    b_lagged = np.zeros(n_spots)
    for i in range(n_spots):
        if len(neighbors[i]) > 0:
            b_lagged[i] = np.mean(b_centered[neighbors[i]])

    # Bivariate Moran's I = correlation between A and spatially-lagged B
    # This captures: "When A is high, are neighboring B values also high?"
    bivariate_i = np.mean(a_centered * b_lagged)

    return float(bivariate_i)


def _compute_bivariate_morans_i_with_weights(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    spatial_weights: NDArray[np.floating],
) -> float:
    """
    Fast bivariate Moran's I using precomputed spatial weights matrix.

    This is optimized for repeated calls during tree building where the
    spatial structure doesn't change but markers do.

    Args:
        values_a: Expression values for marker A (n_spots,).
        values_b: Expression values for marker B (n_spots,).
        spatial_weights: Row-normalized spatial weights matrix (n_spots, n_spots).
            W[i,j] > 0 if j is a neighbor of i, rows should sum to 1.

    Returns:
        Bivariate Moran's I in range approximately [-1, 1].
    """
    # Center and standardize
    mean_a = np.mean(values_a)
    mean_b = np.mean(values_b)
    std_a = np.std(values_a)
    std_b = np.std(values_b)

    if std_a < 1e-10 or std_b < 1e-10:
        return 0.0

    a_centered = (values_a - mean_a) / std_a
    b_centered = (values_b - mean_b) / std_b

    # Spatially-lagged B using matrix multiplication (fast!)
    b_lagged = spatial_weights @ b_centered

    # Bivariate Moran's I
    bivariate_i = np.mean(a_centered * b_lagged)

    return float(bivariate_i)


def _build_spatial_weights_matrix(
    coords: NDArray[np.floating],
    k: int = 6,
) -> NDArray[np.floating]:
    """
    Build row-normalized spatial weights matrix from coordinates.

    Args:
        coords: Spatial coordinates (n_spots, 2).
        k: Number of nearest neighbors.

    Returns:
        Row-normalized spatial weights matrix (n_spots, n_spots).
    """
    from sklearn.neighbors import NearestNeighbors

    n_spots = coords.shape[0]

    # Find k nearest neighbors
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm='ball_tree')
    nn.fit(coords)
    distances, indices = nn.kneighbors(coords)

    # Build sparse-ish weights matrix (dense for simplicity, could optimize later)
    W = np.zeros((n_spots, n_spots), dtype=np.float64)
    for i in range(n_spots):
        # Skip self (index 0), use neighbors (indices 1:k+1)
        neighbor_indices = indices[i, 1:]
        W[i, neighbor_indices] = 1.0

    # Row-normalize
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0  # Avoid division by zero
    W = W / row_sums

    return W


def _compute_bivariate_morans_i_pvalue(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    neighbors: List[List[int]],
    rng: np.random.Generator,
    n_perm: int = 199,
) -> Tuple[float, float]:
    """
    Compute bivariate Moran's I with permutation-based p-value.

    This tests the null hypothesis that there is no spatial cross-correlation
    between markers A and B. Under the null, shuffling B's spatial positions
    should not change the expected I value.

    Uses a smaller number of permutations (default 199) since we only need
    rough p-value estimates for FDR filtering, not precise values.

    Args:
        values_a: Expression values for marker A (n_spots,).
        values_b: Expression values for marker B (n_spots,).
        neighbors: List of neighbor indices per spot.
        rng: Random number generator.
        n_perm: Number of permutations (default: 199 for ~0.5% resolution).

    Returns:
        Tuple of (observed_I, p_value).
    """
    # Compute observed bivariate Moran's I
    observed_i = _compute_bivariate_morans_i(values_a, values_b, neighbors)

    if n_perm <= 0:
        # No permutation test - return placeholder p-value
        return observed_i, 0.5

    # Permutation test: shuffle B's spatial assignment
    null_distribution = np.zeros(n_perm)
    for p in range(n_perm):
        shuffled_b = rng.permutation(values_b)
        null_distribution[p] = _compute_bivariate_morans_i(values_a, shuffled_b, neighbors)

    # One-sided p-value: fraction of null >= observed (for positive association)
    # For negative association, we'd use null <= observed
    # Use absolute value to detect both positive and negative associations
    p_value = (np.sum(np.abs(null_distribution) >= np.abs(observed_i)) + 1) / (n_perm + 1)

    return observed_i, float(p_value)


def _compute_bivariate_morans_i_multiscale(
    values_a: NDArray[np.floating],
    values_b: NDArray[np.floating],
    multi_scale_neighbors: Dict[int, List[List[int]]],
    rng: np.random.Generator,
    n_perm: int = 199,
    aggregation: str = "max",
) -> Tuple[float, float, Dict[int, float], Dict[int, float], int]:
    """
    Compute bivariate Moran's I at multiple spatial scales.

    Larger neighborhoods may capture colocalization signal better in mixed/scattered
    data where cell types are interspersed rather than clustered.

    Args:
        values_a: Expression values for marker A (n_spots,).
        values_b: Expression values for marker B (n_spots,).
        multi_scale_neighbors: Dict mapping k -> neighbor list for that scale.
        rng: Random number generator.
        n_perm: Number of permutations per scale (default: 199).
        aggregation: How to aggregate across scales - "max", "mean", or "weighted".

    Returns:
        Tuple of:
        - aggregated_i: Final bivariate Moran's I (aggregated across scales)
        - aggregated_pvalue: Combined p-value
        - scale_values: Dict mapping k -> I value at that scale
        - scale_pvalues: Dict mapping k -> p-value at that scale
        - best_scale: Scale (k value) with highest signal
    """
    scale_values = {}
    scale_pvalues = {}

    # Reduce permutations per scale to maintain reasonable runtime
    effective_n_perm = max(49, n_perm // len(multi_scale_neighbors))

    for k, neighbors_k in multi_scale_neighbors.items():
        i_val, p_val = _compute_bivariate_morans_i_pvalue(
            values_a, values_b, neighbors_k, rng, n_perm=effective_n_perm
        )
        scale_values[k] = i_val
        scale_pvalues[k] = p_val

    # Aggregation
    scales = list(scale_values.keys())
    values_array = np.array([scale_values[k] for k in scales])
    pvalues_array = np.array([scale_pvalues[k] for k in scales])

    if aggregation == "max":
        # Use the scale with highest bivariate Moran's I
        best_idx = np.argmax(values_array)
        aggregated_i = values_array[best_idx]
        aggregated_pvalue = pvalues_array[best_idx]
        best_scale = scales[best_idx]
    elif aggregation == "mean":
        # Average across scales
        aggregated_i = np.mean(values_array)
        # Fisher's method for combining p-values
        pvalues_clipped = np.clip(pvalues_array, 1e-10, 1.0)
        chi2_stat = -2 * np.sum(np.log(pvalues_clipped))
        from scipy.stats import chi2
        aggregated_pvalue = 1 - chi2.cdf(chi2_stat, df=2 * len(pvalues_array))
        best_scale = scales[np.argmax(values_array)]
    elif aggregation == "weighted":
        # Weight by absolute value (higher |I| = more weight)
        weights = np.abs(values_array) + 0.01
        weights = weights / weights.sum()
        aggregated_i = np.sum(weights * values_array)
        aggregated_pvalue = np.min(pvalues_array)  # Most conservative
        best_scale = scales[np.argmax(values_array)]
    else:
        raise ValueError(f"Invalid aggregation method: {aggregation}")

    return aggregated_i, aggregated_pvalue, scale_values, scale_pvalues, best_scale


def _compute_neighbor_enrichment(
    binary_a: NDArray[np.bool_],
    binary_b: NDArray[np.bool_],
    neighbors: List[List[int]],
    rng: np.random.Generator,
    n_perm: int,
) -> Tuple[float, float, float]:
    """
    Compute enrichment of B in neighbors of A-positive spots.

    Enrichment = (observed B+ neighbors of A+ spots) / (expected by chance)

    Also computes enrichment of A in neighbors of B-positive spots.

    Args:
        binary_a: Binary indicator for marker A (n_spots,).
        binary_b: Binary indicator for marker B (n_spots,).
        neighbors: List of neighbor indices per spot.
        rng: Random number generator.
        n_perm: Number of permutations for p-value.

    Returns:
        Tuple of (enrichment_ab, enrichment_ba, p_value).
    """
    def compute_enrichment_one_way(source_binary, target_binary):
        """Compute enrichment of target in neighbors of source-positive spots."""
        source_positive = np.where(source_binary)[0]
        if len(source_positive) == 0:
            return 0.0, np.array([])

        # Count target-positive neighbors of source-positive spots
        observed = 0
        total_neighbors = 0
        for spot in source_positive:
            spot_neighbors = neighbors[spot]
            if len(spot_neighbors) > 0:
                observed += np.sum(target_binary[spot_neighbors])
                total_neighbors += len(spot_neighbors)

        if total_neighbors == 0:
            return 0.0, np.array([])

        # Expected by chance (background rate of target+)
        expected = np.mean(target_binary) * total_neighbors
        enrichment = observed / expected if expected > 0 else 0.0

        # Permutation null
        null_enrichments = np.zeros(n_perm)
        for i in range(n_perm):
            perm_target = rng.permutation(target_binary)
            perm_observed = 0
            for spot in source_positive:
                spot_neighbors = neighbors[spot]
                if len(spot_neighbors) > 0:
                    perm_observed += np.sum(perm_target[spot_neighbors])
            null_enrichments[i] = perm_observed / expected if expected > 0 else 0.0

        return enrichment, null_enrichments

    # A -> B enrichment
    enrichment_ab, null_ab = compute_enrichment_one_way(binary_a, binary_b)

    # B -> A enrichment
    enrichment_ba, null_ba = compute_enrichment_one_way(binary_b, binary_a)

    # Combined p-value (two-sided)
    if len(null_ab) > 0 and len(null_ba) > 0:
        # Test if observed is more extreme than null
        p_ab = (1 + np.sum(null_ab >= enrichment_ab)) / (1 + len(null_ab))
        p_ba = (1 + np.sum(null_ba >= enrichment_ba)) / (1 + len(null_ba))
        # Use the more conservative p-value
        p_value = max(p_ab, p_ba)
    else:
        p_value = 1.0

    return float(enrichment_ab), float(enrichment_ba), float(p_value)


def _compute_colocalization_score(
    spearman_rho: float,
    cosine_similarity: float,
    bivariate_morans_i: float,
) -> float:
    """
    Combined colocalization score using continuous (binarization-free) metrics.

    All metrics are continuous and scale-invariant, avoiding arbitrary thresholds.

    Weights:
    - Spearman correlation: 0.3 (rank-based correlation, handles non-linear)
    - Cosine similarity: 0.3 (scale-invariant pattern matching)
    - Bivariate Moran's I: 0.4 (spatial cross-correlation, most informative)

    Args:
        spearman_rho: Spearman rank correlation [-1, 1].
        cosine_similarity: Cosine similarity [-1, 1].
        bivariate_morans_i: Spatial cross-correlation [-1, 1].

    Returns:
        Combined colocalization score in [0, 1].
    """
    # Normalize all to [0, 1] range (from [-1, 1])
    norm_spearman = (spearman_rho + 1) / 2
    norm_cosine = (cosine_similarity + 1) / 2
    norm_bivariate = (bivariate_morans_i + 1) / 2

    return 0.3 * norm_spearman + 0.3 * norm_cosine + 0.4 * norm_bivariate


def _process_marker_pair(
    i: int,
    j: int,
    analyze_names: List[str],
    analyze_X: NDArray[np.floating],
    analyze_X_smooth: NDArray[np.floating],
    binary: NDArray[np.bool_],
    neighbors: List[List[int]],
    n_permutations: int,
    seed: int,
    multi_scale_neighbors: Optional[Dict[int, List[List[int]]]] = None,
    multi_scale_aggregation: str = "max",
) -> MarkerPairColocalization:
    """
    Process a single marker pair for colocalization analysis.

    This function is designed to be called in parallel across all pairs.

    Args:
        i, j: Indices of the two markers.
        analyze_names: List of marker names.
        analyze_X: Raw expression matrix for markers being analyzed.
        analyze_X_smooth: Spatially smoothed expression matrix.
        binary: Binarized marker matrix.
        neighbors: Neighbor list per spot.
        n_permutations: Number of permutations for significance testing.
        seed: Random seed (will be adjusted by pair indices for reproducibility).
        multi_scale_neighbors: Optional dict mapping k -> neighbor list for multi-scale analysis.
        multi_scale_aggregation: Aggregation method for multi-scale ("max", "mean", "weighted").

    Returns:
        MarkerPairColocalization object for this pair.
    """
    # Use deterministic seed based on pair indices for reproducibility
    pair_seed = seed + i * 1000 + j
    rng = np.random.default_rng(pair_seed)

    name_a = analyze_names[i]
    name_b = analyze_names[j]

    # Raw values for correlation and score metrics
    values_a = analyze_X[:, i]
    values_b = analyze_X[:, j]
    # Smoothed values for bivariate Moran's I (better for mixed data)
    smooth_a = analyze_X_smooth[:, i]
    smooth_b = analyze_X_smooth[:, j]
    binary_a = binary[:, i]
    binary_b = binary[:, j]

    # Same-spot metrics (binary-based, for backwards compatibility)
    jaccard, co_spots, co_frac = _compute_jaccard(binary_a, binary_b)

    # Correlation metrics (continuous, on raw data)
    pearson_r, spearman_rho, corr_pvalue = _compute_correlation(values_a, values_b)

    # Continuous similarity metrics (binarization-free) - used for score
    cosine_sim = _compute_cosine_similarity(values_a, values_b)

    # Bivariate Moran's I - single-scale or multi-scale
    if multi_scale_neighbors is not None and len(multi_scale_neighbors) > 1:
        # Multi-scale analysis: compute at multiple k values and aggregate
        bivariate_i, bivariate_pval, scale_values, scale_pvalues, best_scale = \
            _compute_bivariate_morans_i_multiscale(
                smooth_a, smooth_b, multi_scale_neighbors, rng,
                n_perm=n_permutations, aggregation=multi_scale_aggregation
            )
    else:
        # Single-scale (backward compatible)
        bivariate_i, bivariate_pval = _compute_bivariate_morans_i_pvalue(
            smooth_a, smooth_b, neighbors, rng, n_perm=n_permutations
        )
        scale_values = None
        scale_pvalues = None
        best_scale = None

    # Neighborhood metrics (binary-based, alternative for high-seg data)
    # Binarization focuses on PEAK expression spots - stricter but works for clear boundaries
    enrich_ab, enrich_ba, enrich_pvalue = _compute_neighbor_enrichment(
        binary_a, binary_b, neighbors, rng, n_permutations
    )
    mutual_enrich = (enrich_ab + enrich_ba) / 2

    # Combined score using continuous metrics only
    score = _compute_colocalization_score(spearman_rho, cosine_sim, bivariate_i)

    return MarkerPairColocalization(
        marker_a=name_a,
        marker_b=name_b,
        jaccard_index=jaccard,
        co_occurrence_spots=co_spots,
        co_occurrence_fraction=co_frac,
        pearson_r=pearson_r,
        spearman_rho=spearman_rho,
        correlation_pvalue=corr_pvalue,
        cosine_similarity=cosine_sim,
        bivariate_morans_i=bivariate_i,
        bivariate_morans_pvalue=bivariate_pval,
        neighbor_enrichment_ab=enrich_ab,
        neighbor_enrichment_ba=enrich_ba,
        neighbor_enrichment_pvalue=enrich_pvalue,
        mutual_neighbor_enrichment=mutual_enrich,
        colocalization_score=score,
        bivariate_morans_per_scale=scale_values,
        bivariate_morans_pvalue_per_scale=scale_pvalues,
        bivariate_morans_best_scale=best_scale,
    )


def _make_names_unique(names: List[str]) -> Tuple[List[str], Dict[str, str]]:
    """
    Make marker names unique by appending suffixes to duplicates.

    Args:
        names: List of marker names (may contain duplicates).

    Returns:
        Tuple of (unique_names, mapping from unique back to original).
    """
    counts: Dict[str, int] = {}
    unique_names = []
    unique_to_original = {}

    for name in names:
        if name not in counts:
            counts[name] = 0
            unique_name = name
        else:
            counts[name] += 1
            unique_name = f"{name}_{counts[name]}"

        unique_names.append(unique_name)
        unique_to_original[unique_name] = name

    return unique_names, unique_to_original


def analyze_marker_colocalization(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    *,
    markers_to_analyze: Optional[List[str]] = None,
    neighbor_k: int = 6,
    smooth_k: int = 6,
    signal_threshold_percentile: float = 75.0,
    n_permutations: int = 999,
    seed: int = 1234,
    verbose: bool = True,
    # Multi-scale neighborhood parameters
    multi_scale_k: Optional[List[int]] = [6, 12, 24, 48, 64],
    multi_scale_aggregation: str = "max",
) -> ColocalizationResult:
    """
    Analyze pairwise spatial colocalization between markers.

    Computes three types of relationships:
    1. Same-spot: Jaccard index and co-occurrence frequency
    2. Correlation: Pearson and Spearman correlation of intensities
    3. Spatial: Bivariate Moran's I (spatial cross-correlation) with permutation test

    The bivariate Moran's I p-value is used for FDR filtering in profile discovery,
    as it directly tests whether markers co-localize in spatial patterns.

    Args:
        X: Expression matrix (n_spots, n_markers).
        coords: Spatial coordinates for each spot (n_spots, 2).
        marker_names: Names for each marker column.
        markers_to_analyze: Subset of markers to analyze. If None, analyze all pairwise.
            Recommended: Pass result.interesting_markers from identify_interesting_markers().
        neighbor_k: Number of nearest neighbors for spatial analysis (default: 6).
        smooth_k: Number of neighbors for spatial smoothing before bivariate Moran's I (default: 6).
            Smoothing reduces noise and improves detection of spatial co-expression patterns,
            especially in mixed tissue contexts where cell types overlap spatially.
        signal_threshold_percentile: Percentile threshold for classifying positive spots (default: 75).
        n_permutations: Number of permutations for bivariate Moran's I p-value (default: 199).
            Higher values give more precise p-values but take longer.
        seed: Random seed for reproducibility (default: 1234).
        verbose: Log progress information (default: True).
        multi_scale_k: List of k values for multi-scale neighborhood analysis
            (default: [6, 12, 24, 48, 64]). Computes bivariate Moran's I at each scale
            and aggregates according to multi_scale_aggregation. Set to None to disable.
            Scales >= n_spots are automatically filtered out.
        multi_scale_aggregation: How to aggregate across scales - "max" (default),
            "mean", or "weighted". "max" selects the best signal at optimal scale.

    Returns:
        ColocalizationResult with DataFrame of marker pairs ranked by colocalization score.

    Example:
        >>> # First identify interesting markers
        >>> interest_result = identify_interesting_markers(X, coords, marker_names)
        >>> # Then analyze colocalization among interesting markers
        >>> coloc_result = analyze_marker_colocalization(
        ...     X, coords, marker_names,
        ...     markers_to_analyze=interest_result.interesting_markers
        ... )
        >>> df = coloc_result.to_dataframe()
        >>> df.head(10)  # Top 10 colocalizing pairs
    """
    X = np.asarray(X, dtype=np.float64)
    coords = np.asarray(coords, dtype=np.float64)

    n_spots, n_markers = X.shape

    if len(marker_names) != n_markers:
        raise ValueError(
            f"Number of marker names ({len(marker_names)}) must match "
            f"number of columns in X ({n_markers})"
        )

    # Handle duplicate marker names by making them unique
    marker_names_unique, unique_to_original = _make_names_unique(list(marker_names))
    has_duplicates = len(set(marker_names)) != len(marker_names)
    if has_duplicates and verbose:
        logging.info(f"Found duplicate marker names, making unique: {[m for m in marker_names_unique if m != unique_to_original[m]]}")

    # Filter to markers of interest
    # Use unique names for internal analysis
    if markers_to_analyze is not None:
        # Match against both original and unique names
        markers_to_analyze_set = set(markers_to_analyze)
        marker_indices = [
            i for i, (orig, uniq) in enumerate(zip(marker_names, marker_names_unique))
            if orig in markers_to_analyze_set or uniq in markers_to_analyze_set
        ]
        if len(marker_indices) == 0:
            raise ValueError("No markers from markers_to_analyze found in marker_names")
        analyze_names = [marker_names_unique[i] for i in marker_indices]
        analyze_X = X[:, marker_indices]
    else:
        marker_indices = list(range(n_markers))
        analyze_names = list(marker_names_unique)
        analyze_X = X

    n_analyze = len(analyze_names)

    if verbose:
        n_pairs = n_analyze * (n_analyze - 1) // 2
        logging.info(
            f"Analyzing colocalization for {n_analyze} markers ({n_pairs} pairs) "
            f"across {n_spots} spots"
        )

    # Build neighbor graph
    if verbose:
        logging.info(f"Building {neighbor_k}-NN neighbor graph...")
    neighbors = _build_neighbor_graph(coords, neighbor_k)

    # Spatial smoothing for bivariate Moran's I (reduces noise, improves mixed data detection)
    if smooth_k > 0:
        if verbose:
            logging.info(f"Applying spatial smoothing (k={smooth_k} neighbors) for bivariate Moran's I...")
        smooth_neighbors = _build_neighbor_graph(coords, smooth_k)
        analyze_X_smooth = np.zeros_like(analyze_X)
        for spot_idx in range(n_spots):
            neighbor_vals = analyze_X[smooth_neighbors[spot_idx], :]
            analyze_X_smooth[spot_idx] = np.mean(neighbor_vals, axis=0)
    else:
        analyze_X_smooth = analyze_X

    # Binarize markers (kept for backwards-compatible metrics)
    if verbose:
        logging.info(f"Binarizing markers at {signal_threshold_percentile}th percentile...")
    binary = _binarize_markers(analyze_X, signal_threshold_percentile)

    # Compute pairwise colocalization
    total_pairs = n_analyze * (n_analyze - 1) // 2

    if verbose:
        logging.info(f"Computing pairwise colocalization (continuous metrics)...")

    # Generate all pair indices
    pair_indices = [(i, j) for i in range(n_analyze) for j in range(i + 1, n_analyze)]

    # Determine number of workers (use environment variable or default to -1 for all CPUs)
    n_jobs = int(os.environ.get('CITEGEIST_N_JOBS', -1))
    if verbose:
        n_cpus = os.cpu_count() or 1
        logging.info(f"Running parallel colocalization analysis with {n_cpus if n_jobs == -1 else n_jobs} workers...")

    # Build multi-scale neighbor graphs if requested
    if multi_scale_k is not None and len(multi_scale_k) > 1:
        # Filter out k values >= n_spots (can't have more neighbors than spots)
        valid_scales = [k for k in multi_scale_k if k < n_spots]
        if len(valid_scales) < len(multi_scale_k) and verbose:
            skipped = [k for k in multi_scale_k if k >= n_spots]
            logging.warning(f"Skipping scales {skipped} (>= n_spots={n_spots})")
        if len(valid_scales) > 1:
            if verbose:
                logging.info(f"Building multi-scale neighbor graphs (k={valid_scales})...")
            multi_scale_neighbors = {
                k: _build_neighbor_graph(coords, k) for k in valid_scales
            }
            if verbose:
                logging.info(f"Multi-scale aggregation method: {multi_scale_aggregation}")
        else:
            # Fall back to single-scale if not enough valid scales
            multi_scale_neighbors = None
            if verbose:
                logging.info("Only one valid scale available, using single-scale mode")
    else:
        multi_scale_neighbors = None

    # Run in parallel using joblib
    pairs = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_process_marker_pair)(
            i, j, analyze_names, analyze_X, analyze_X_smooth, binary,
            neighbors, n_permutations, seed,
            multi_scale_neighbors, multi_scale_aggregation
        )
        for i, j in pair_indices
    )

    if verbose:
        logging.info(f"Processed {len(pairs)}/{total_pairs} pairs")

    # Sort by score
    pairs.sort(key=lambda x: x.colocalization_score, reverse=True)

    # Compute marker specificity (Gini coefficient)
    # Higher = more specific (concentrated in few spots), Lower = more promiscuous
    marker_specificity = _compute_marker_specificity(analyze_X, analyze_names)

    if verbose:
        # Log promiscuous markers (specificity < 0.3)
        promiscuous = [m for m, s in marker_specificity.items() if s < 0.3]
        specific = [m for m, s in marker_specificity.items() if s > 0.6]
        if promiscuous:
            logging.info(f"Promiscuous markers (Gini < 0.3): {promiscuous}")
        if specific:
            logging.info(f"Highly specific markers (Gini > 0.6): {specific[:5]}...")

    result = ColocalizationResult(
        pairs=pairs,
        marker_names=analyze_names,
        n_spots=n_spots,
        neighbor_k=neighbor_k,
        marker_specificity=marker_specificity,
    )

    if verbose:
        top_pair = pairs[0] if pairs else None
        if top_pair:
            logging.info(
                f"Top colocalization: {top_pair.marker_a} <-> {top_pair.marker_b} "
                f"(score={top_pair.colocalization_score:.3f})"
            )

    return result


# =============================================================================
# Profile Discovery
# =============================================================================


@dataclass
class LineageDendrogram:
    """Hierarchical clustering result for a single lineage (connected component)."""

    markers: List[str]  # Markers in this lineage
    linkage_matrix: NDArray[np.floating]  # scipy linkage matrix
    distance_matrix: NDArray[np.floating]  # Pairwise distances used

    def get_newick(self) -> str:
        """Convert to Newick format for visualization."""
        from scipy.cluster.hierarchy import to_tree

        def _to_newick(node, labels):
            if node.is_leaf():
                return labels[node.id]
            else:
                left = _to_newick(node.get_left(), labels)
                right = _to_newick(node.get_right(), labels)
                return f"({left}:{node.dist/2:.3f},{right}:{node.dist/2:.3f})"

        tree = to_tree(self.linkage_matrix)
        return _to_newick(tree, self.markers) + ";"


@dataclass
class ProfileDiscoveryResult:
    """Results from profile discovery."""

    profiles: List[List[str]]  # List of marker sets (each is a cell type profile)
    lineage_dendrograms: Dict[int, LineageDendrogram]  # One dendrogram per lineage
    singletons: List[str]  # Markers that form their own profile (no significant edges)
    modularity: float  # How well profiles explain the colocalization structure
    n_significant_edges: int  # Number of edges that passed significance filter
    alpha: float  # Significance threshold used

    def to_dataframe(self) -> pd.DataFrame:
        """Convert profiles to DataFrame."""
        records = []
        for i, profile in enumerate(self.profiles):
            records.append({
                "profile_id": i,
                "n_markers": len(profile),
                "markers": ", ".join(sorted(profile)),
                "is_singleton": len(profile) == 1,
            })
        return pd.DataFrame(records)

    def get_profile_for_marker(self, marker: str) -> Optional[List[str]]:
        """Get the profile containing a specific marker."""
        for profile in self.profiles:
            if marker in profile:
                return profile
        return None

    def summary(self) -> str:
        """Return a summary string."""
        n_profiles = len(self.profiles)
        n_singletons = len(self.singletons)
        n_multi = n_profiles - n_singletons
        sizes = [len(p) for p in self.profiles]
        return (
            f"Discovered {n_profiles} profiles: {n_multi} multi-marker, {n_singletons} singletons\n"
            f"Profile sizes: {sorted(sizes, reverse=True)}\n"
            f"Modularity: {self.modularity:.3f}\n"
            f"Significant edges: {self.n_significant_edges}"
        )


@dataclass
class ProfileTreeNode:
    """A node in the hierarchical profile tree."""

    node_id: str
    markers: List[str]  # Markers assigned to THIS node (not descendants)
    children: List["ProfileTreeNode"]
    parent_id: Optional[str]
    depth: int
    nmf_weights: Optional[Dict[str, float]] = None  # marker -> weight from NMF

    @property
    def is_leaf(self) -> bool:
        """True if this is a leaf node (no children)."""
        return len(self.children) == 0

    def get_all_markers(self) -> List[str]:
        """Get all markers in this subtree (self + descendants)."""
        markers = list(self.markers)
        for child in self.children:
            markers.extend(child.get_all_markers())
        return markers


@dataclass
class ProfileTree:
    """Hierarchical tree of marker profiles."""

    root: ProfileTreeNode
    n_levels: int

    def get_leaves(self) -> List[ProfileTreeNode]:
        """Return all leaf nodes (final cell type profiles)."""
        leaves: List[ProfileTreeNode] = []
        self._collect_leaves(self.root, leaves)
        return leaves

    def _collect_leaves(self, node: ProfileTreeNode, leaves: List[ProfileTreeNode]) -> None:
        """Recursively collect leaf nodes."""
        if node.is_leaf:
            leaves.append(node)
        else:
            for child in node.children:
                self._collect_leaves(child, leaves)

    def get_depth(self) -> int:
        """Return maximum depth of tree."""
        return self._max_depth(self.root)

    def _max_depth(self, node: ProfileTreeNode) -> int:
        """Recursively compute max depth."""
        if node.is_leaf:
            return node.depth
        return max(self._max_depth(c) for c in node.children)


@dataclass
class HierarchicalProfileResult:
    """Results from hierarchical profile discovery.

    This container holds the complete hierarchical tree structure along with
    flattened profiles that are compatible with downstream Module 3
    deconvolution.

    Attributes:
        tree: Full hierarchical tree for inspection and visualization.
        flat_profiles: Cell type -> list of markers, for Module 3 compatibility.
        depth_per_branch: Branch name -> depth, indicating hierarchy complexity.
        shared_markers: Marker -> list of cell types that share this marker.
        reconstruction_error: Final reconstruction error from profile fitting.
    """

    tree: ProfileTree  # Full hierarchy for inspection
    flat_profiles: Dict[str, List[str]]  # cell_type -> [markers] for Module 3
    depth_per_branch: Dict[str, int]  # branch_name -> depth
    shared_markers: Dict[str, List[str]]  # marker -> [cell_types that share it]
    reconstruction_error: float  # Final reconstruction error

    def to_profile_dict(self) -> Dict[str, List[str]]:
        """Convert to Module 3 compatible profile dictionary.

        Returns:
            Dictionary mapping cell type names to marker lists, suitable
            for use with CitegeistModel.load_cell_profile_dict().
        """
        return dict(self.flat_profiles)

    def summary(self) -> str:
        """Return a human-readable summary string.

        Returns:
            Multi-line string with key statistics about the hierarchical
            profiles.
        """
        n_profiles = len(self.flat_profiles)
        n_shared = len(self.shared_markers)
        max_depth = self.tree.get_depth()
        return (
            f"Hierarchical profiles: {n_profiles} cell types\n"
            f"Tree depth: {max_depth}\n"
            f"Shared markers: {n_shared}\n"
            f"Reconstruction error: {self.reconstruction_error:.4f}"
        )


def _compute_nmf_weights(
    X: NDArray[np.floating],
    marker_names: List[str],
    node: ProfileTreeNode,
) -> Dict[str, Dict[str, float]]:
    """
    Compute NMF weights at an internal node.

    Runs NMF with n_components = number of children to learn
    how markers should be allocated to each child branch.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        node: Internal node with children

    Returns:
        Dict mapping child_id -> {marker: weight}
    """
    from sklearn.decomposition import NMF

    if node.is_leaf or len(node.children) == 0:
        return {}

    # Get all markers in this subtree
    all_markers = node.get_all_markers()
    if len(all_markers) == 0:
        return {c.node_id: {} for c in node.children}

    marker_to_idx = {m: i for i, m in enumerate(marker_names)}
    indices = [marker_to_idx[m] for m in all_markers if m in marker_to_idx]

    if len(indices) == 0:
        return {c.node_id: {} for c in node.children}

    # Extract submatrix
    X_node = X[:, indices]
    X_node = np.maximum(X_node, 0)

    # Number of components = number of children
    n_components = len(node.children)
    n_components = min(n_components, len(indices), X_node.shape[0] - 1)
    n_components = max(1, n_components)

    try:
        nmf = NMF(n_components=n_components, init='nndsvda', max_iter=200, random_state=42)
        W = nmf.fit_transform(X_node)  # (n_spots, n_components)
        H = nmf.components_  # (n_components, n_markers_in_node)

        # Map components to children based on which markers they load on
        # H[k, j] = how much component k contributes to marker j

        # Build list of markers in node order (matching X_node columns)
        markers_in_node = [all_markers[i] for i in range(len(all_markers))
                          if all_markers[i] in marker_to_idx]

        # For each child, find which component best matches its markers
        child_to_component = {}

        for child_idx, child in enumerate(node.children):
            child_markers = child.get_all_markers()
            child_marker_indices = [markers_in_node.index(m) for m in child_markers
                                    if m in markers_in_node]

            if len(child_marker_indices) == 0:
                child_to_component[child.node_id] = child_idx % n_components
                continue

            # Find component with highest average loading on child markers
            best_component = 0
            best_score = -np.inf
            for comp_idx in range(n_components):
                score = np.mean(H[comp_idx, child_marker_indices])
                if score > best_score:
                    best_score = score
                    best_component = comp_idx

            child_to_component[child.node_id] = best_component

        # Build weights dict
        weights: Dict[str, Dict[str, float]] = {}
        for child in node.children:
            comp_idx = child_to_component[child.node_id]
            weights[child.node_id] = {}
            for marker_idx, marker in enumerate(markers_in_node):
                if marker_idx < H.shape[1]:
                    weights[child.node_id][marker] = float(H[comp_idx, marker_idx])

        # Normalize weights per marker (so they sum to 1 across children)
        for marker in markers_in_node:
            total = sum(weights[c.node_id].get(marker, 0) for c in node.children)
            if total > 1e-10:
                for child in node.children:
                    if marker in weights[child.node_id]:
                        weights[child.node_id][marker] /= total

        return weights

    except Exception as e:
        logger.warning(f"NMF failed at node {node.node_id}: {e}")
        # Fallback: equal weights
        weights = {}
        for child in node.children:
            child_markers = child.get_all_markers()
            weights[child.node_id] = {m: 1.0 / len(node.children) for m in child_markers}
        return weights


def _apply_fdr_correction(
    pairs: List[MarkerPairColocalization],
    alpha: float = 0.05,
    fallback_to_raw: bool = True,
    pvalue_source: str = "bivariate_morans",
) -> Tuple[List[MarkerPairColocalization], NDArray[np.bool_], str]:
    """
    Apply Benjamini-Hochberg FDR correction to p-values.

    Args:
        pairs: List of marker pair colocalization results.
        alpha: FDR threshold (default: 0.05 for q < 5%).
        fallback_to_raw: If True and FDR finds nothing but raw p-values
            suggest signal exists, fall back to raw p-value threshold.
        pvalue_source: Which p-value to use for FDR correction:
            - 'bivariate_morans' (default): Tests spatial cross-correlation
              between markers. Works well for both mixed and high-seg data.
            - 'neighbor_enrichment': Binary-based, tests if marker B is enriched
              in neighbors of A-positive spots. Best for high-segregation data.
            - 'correlation': Uses correlation p-value (Pearson/Spearman).

    Returns:
        Tuple of (filtered pairs, boolean mask, method_used).
        method_used is "fdr" or "raw_pvalue".
    """
    if len(pairs) == 0:
        return [], np.array([], dtype=bool), "fdr"

    # Select p-value source based on data characteristics
    if pvalue_source == "bivariate_morans":
        pvalues = np.array([p.bivariate_morans_pvalue for p in pairs])
    elif pvalue_source == "correlation":
        pvalues = np.array([p.correlation_pvalue for p in pairs])
    else:  # neighbor_enrichment
        pvalues = np.array([p.neighbor_enrichment_pvalue for p in pairs])

    # Benjamini-Hochberg procedure
    n = len(pvalues)
    sorted_indices = np.argsort(pvalues)
    sorted_pvalues = pvalues[sorted_indices]

    # BH critical values: (rank / n) * alpha
    bh_critical = (np.arange(1, n + 1) / n) * alpha

    # Find largest k where p_(k) <= (k/n) * alpha
    below_threshold = sorted_pvalues <= bh_critical

    if np.any(below_threshold):
        # FDR found significant pairs
        max_k = np.max(np.where(below_threshold)[0])
        significant_sorted_indices = sorted_indices[:max_k + 1]

        is_significant = np.zeros(n, dtype=bool)
        is_significant[significant_sorted_indices] = True

        filtered_pairs = [p for p, sig in zip(pairs, is_significant) if sig]
        return filtered_pairs, is_significant, "fdr"

    # FDR found nothing - check if this is a permutation resolution issue
    min_pval = sorted_pvalues[0]
    n_below_alpha = np.sum(pvalues < alpha)

    if fallback_to_raw and n_below_alpha > 0:
        # There are pairs with p < alpha, but FDR is too strict
        # This happens when permutation resolution (1/n_perms) > alpha/n_tests
        # Fall back to raw p-value threshold

        # Use the minimum p-value as threshold (all pairs at min_pval are signal)
        # This is conservative: only keep pairs with the best possible p-value
        is_significant = pvalues <= min_pval

        logging.warning(
            f"FDR correction found 0 significant pairs, but {n_below_alpha} have p < {alpha}. "
            f"Permutation resolution (min_p={min_pval:.4f}) may be insufficient for "
            f"{n} tests. Falling back to raw p-value threshold (p <= {min_pval:.4f})."
        )

        filtered_pairs = [p for p, sig in zip(pairs, is_significant) if sig]
        return filtered_pairs, is_significant, "raw_pvalue"

    # No discoveries
    return [], np.zeros(n, dtype=bool), "fdr"


def _apply_mutual_topk(
    pairs: List[MarkerPairColocalization],
    k: int = 3,
    specificity: Optional[Dict[str, float]] = None,
) -> List[MarkerPairColocalization]:
    """
    Apply mutual top-k sparsification to edges.

    For each marker, keep only its top-k partners by colocalization score.
    Then keep an edge only if BOTH markers rank each other in their top-k.

    This prevents the "one giant component" collapse by limiting hub markers.

    Args:
        pairs: List of marker pair colocalization results.
        k: Number of top partners to keep per marker (default: 3).
        specificity: Optional dict mapping marker names to specificity scores [0,1].
            If provided, effective score = score * sqrt(s_a * s_b), which
            downweights edges involving broad/ubiquitous markers.

    Returns:
        Filtered list of pairs that pass mutual top-k criterion.
    """
    if len(pairs) == 0 or k <= 0:
        return pairs

    # Build marker -> [(partner, effective_score, pair)] mapping
    marker_edges: Dict[str, List[Tuple[str, float, MarkerPairColocalization]]] = defaultdict(list)
    for p in pairs:
        score = p.colocalization_score
        # Apply specificity weighting if provided
        if specificity is not None:
            s_a = specificity.get(p.marker_a, 0.5)
            s_b = specificity.get(p.marker_b, 0.5)
            score = score * np.sqrt(s_a * s_b)
        marker_edges[p.marker_a].append((p.marker_b, score, p))
        marker_edges[p.marker_b].append((p.marker_a, score, p))

    # For each marker, keep top-k by score
    topk_partners: Dict[str, Set[str]] = {}
    for marker, edges in marker_edges.items():
        # Sort by score descending, take top k
        sorted_edges = sorted(edges, key=lambda x: x[1], reverse=True)[:k]
        topk_partners[marker] = {e[0] for e in sorted_edges}

    # Keep only mutual edges (both rank each other in top-k)
    kept = []
    seen: Set[Tuple[str, str]] = set()
    for p in pairs:
        key = tuple(sorted([p.marker_a, p.marker_b]))
        if key in seen:
            continue

        # Check mutual: A has B in top-k AND B has A in top-k
        a_has_b = p.marker_b in topk_partners.get(p.marker_a, set())
        b_has_a = p.marker_a in topk_partners.get(p.marker_b, set())

        if a_has_b and b_has_a:
            kept.append(p)
            seen.add(key)

    return kept


def _compute_marker_specificity(
    X: NDArray[np.floating],
    marker_names: List[str],
) -> Dict[str, float]:
    """
    Compute Gini coefficient for each marker (higher = more specific).

    Specific markers (high Gini) are concentrated in few spots.
    Broad markers (low Gini) are spread across many spots.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column

    Returns:
        Dict mapping marker name to specificity score [0, 1]
    """
    specificity = {}
    for i, name in enumerate(marker_names):
        values = X[:, i]
        values = values[values > 0]  # Only positive values
        if len(values) < 2:
            specificity[name] = 0.5  # Default for rare markers
            continue

        # Gini coefficient: 0 = uniform, 1 = concentrated
        sorted_vals = np.sort(values)
        n = len(sorted_vals)
        total = sorted_vals.sum()
        if total == 0:
            specificity[name] = 0.5
            continue

        # Gini = 1 - 2 * (area under Lorenz curve)
        cumsum = np.cumsum(sorted_vals)
        gini = (n + 1 - 2 * np.sum(cumsum) / total) / n
        specificity[name] = max(0, min(1, gini))

    return specificity


def _find_connected_components(
    markers: List[str],
    edges: List[Tuple[str, str, float]],
) -> List[Set[str]]:
    """
    Find connected components in the marker graph.

    Uses union-find algorithm.

    Args:
        markers: List of all marker names.
        edges: List of (marker_a, marker_b, weight) tuples.

    Returns:
        List of sets, each containing markers in one component.
    """
    # Initialize each marker as its own component
    parent = {m: m for m in markers}

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Union markers connected by edges
    for m1, m2, _ in edges:
        union(m1, m2)

    # Group by root
    components = defaultdict(set)
    for m in markers:
        components[find(m)].add(m)

    return list(components.values())


def _build_distance_matrix(
    markers: List[str],
    pairs: List[MarkerPairColocalization],
) -> NDArray[np.floating]:
    """
    Build a distance matrix from colocalization scores.

    Distance = 1 - colocalization_score (so high colocalization = low distance).

    Args:
        markers: List of marker names (defines matrix order).
        pairs: Colocalization pairs with scores.

    Returns:
        Distance matrix (n_markers, n_markers).
    """
    n = len(markers)
    marker_to_idx = {m: i for i, m in enumerate(markers)}

    # Initialize with max distance
    dist = np.ones((n, n))
    np.fill_diagonal(dist, 0)

    # Fill in distances from colocalization scores
    for pair in pairs:
        if pair.marker_a in marker_to_idx and pair.marker_b in marker_to_idx:
            i = marker_to_idx[pair.marker_a]
            j = marker_to_idx[pair.marker_b]
            d = 1.0 - pair.colocalization_score
            dist[i, j] = d
            dist[j, i] = d

    return dist


def _find_gap_threshold_gmm(
    gaps: NDArray[np.floating],
    seed: int = 1234,
) -> Tuple[int, float, NDArray[np.bool_]]:
    """
    Use GMM with BIC to find threshold for "large" gaps in merge distances.

    Args:
        gaps: Array of gap values (differences between consecutive merge distances).
        seed: Random seed.

    Returns:
        Tuple of:
        - n_components: Number of GMM components selected by BIC
        - threshold: Gap value separating "normal" from "large" gaps
        - is_large: Boolean array indicating which gaps are "large"
    """
    from sklearn.mixture import GaussianMixture

    gaps = gaps[np.isfinite(gaps) & (gaps > 0)]
    if len(gaps) < 4:
        return 1, np.inf, np.zeros(len(gaps), dtype=bool)

    gaps_2d = gaps.reshape(-1, 1)

    # Fit GMMs with 1-3 components (don't need many for gaps)
    best_bic = np.inf
    best_gmm = None
    best_k = 1

    for k in range(1, 4):
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type='full',
                random_state=seed,
                n_init=3,
                max_iter=200,
            )
            gmm.fit(gaps_2d)
            bic = gmm.bic(gaps_2d)

            if bic < best_bic:
                best_bic = bic
                best_gmm = gmm
                best_k = k
        except Exception:
            continue

    if best_gmm is None or best_k == 1:
        # No separation found - no large gaps
        return 1, np.inf, np.zeros(len(gaps), dtype=bool)

    # Identify the "large gap" component (highest mean)
    means = best_gmm.means_.flatten()
    large_component = np.argmax(means)

    # Get assignments
    labels = best_gmm.predict(gaps_2d)
    is_large = labels == large_component

    # Threshold is midpoint between largest "normal" mean and "large" mean
    sorted_means = np.sort(means)
    threshold = (sorted_means[-2] + sorted_means[-1]) / 2 if len(sorted_means) > 1 else sorted_means[0]

    return best_k, float(threshold), is_large


def _split_dendrogram_by_gaps(
    linkage_matrix: NDArray[np.floating],
    markers: List[str],
    seed: int = 1234,
    min_gap_ratio: float = 2.0,
) -> List[List[str]]:
    """
    Split a dendrogram into separate lineages based on unusually large gaps.

    Only splits if there's a gap that is significantly larger than others,
    indicating a clear lineage boundary.

    Args:
        linkage_matrix: scipy linkage matrix.
        markers: List of marker names (leaves).
        seed: Random seed.
        min_gap_ratio: Minimum ratio of largest gap to median gap to trigger split.
            Default 2.0 means the largest gap must be at least 2x the median.

    Returns:
        List of marker lists, each representing a separate lineage.
    """
    n_markers = len(markers)

    if n_markers <= 3:
        return [markers]

    # Get merge distances
    distances = linkage_matrix[:, 2]

    # Compute gaps between consecutive merges
    gaps = np.diff(distances)

    if len(gaps) < 3:
        return [markers]

    # Find the largest gap
    max_gap_idx = np.argmax(gaps)
    max_gap = gaps[max_gap_idx]

    # Compare to median gap (excluding the max)
    other_gaps = np.concatenate([gaps[:max_gap_idx], gaps[max_gap_idx+1:]])
    if len(other_gaps) == 0:
        return [markers]

    median_gap = np.median(other_gaps)

    # Only split if the largest gap is significantly larger than typical
    if median_gap <= 0 or max_gap / median_gap < min_gap_ratio:
        # No clear lineage boundary
        return [markers]

    # Also check: the gap should be in the upper portion of the tree
    # (early merges are within-profile, late merges are between-lineage)
    # Only consider gaps in the top 50% of merges
    if max_gap_idx < len(gaps) * 0.3:
        # Gap is too early - probably within-lineage variation
        return [markers]

    # Cut at this gap
    cut_distance = distances[max_gap_idx + 1] - 1e-10

    try:
        clusters = fcluster(linkage_matrix, t=cut_distance, criterion='distance')
    except Exception:
        return [markers]

    # Group markers by cluster
    cluster_to_markers = defaultdict(list)
    for i, c in enumerate(clusters):
        cluster_to_markers[c].append(markers[i])

    lineages = list(cluster_to_markers.values())

    # Only accept the split if we get valid lineages (at least 2 markers each)
    # A pair of markers is a valid lineage - represents co-expressed proteins
    if any(len(lin) < 2 for lin in lineages):
        # Split would create singletons - not valid lineage
        return [markers]

    # Log the split for debugging
    logging.info(f"  Gap-split at merge {max_gap_idx} (gap={max_gap:.3f}, {max_gap/median_gap:.1f}x median)")
    for i, lin in enumerate(lineages):
        logging.info(f"    Lineage {i+1}: {lin}")

    # Recursively check each lineage for further splits (max 1 level deep)
    final_lineages = []
    for lineage in lineages:
        if len(lineage) <= 4:
            final_lineages.append(lineage)
        else:
            # Build sub-dendrogram and check for one more split
            lin_dist = _build_distance_matrix(lineage, [])  # Empty pairs = max distance
            # Actually we need the pairs, but we don't have access here
            # Just keep as single lineage for now
            final_lineages.append(lineage)

    return final_lineages


def _dynamic_tree_cut(
    linkage_matrix: NDArray[np.floating],
    n_markers: int,
) -> NDArray[np.int_]:
    """
    Cut dendrogram using gap-based adaptive method.

    Looks for large jumps in merge distances to identify natural cluster boundaries.
    Conservative for small lineages (3-5 markers) that have already passed
    FDR, top-k, and initial gap-splitting criteria.

    Args:
        linkage_matrix: scipy linkage matrix.
        n_markers: Number of markers (leaves).

    Returns:
        Cluster assignments (n_markers,).
    """
    if n_markers <= 1:
        return np.array([0] * n_markers)

    if n_markers == 2:
        # Two markers - check if they should be one cluster or two
        # Use the merge distance - if very low, same cluster
        if linkage_matrix[0, 2] < 0.5:  # distance < 0.5 means score > 0.5
            return np.array([0, 0])
        else:
            return np.array([0, 1])

    # Get merge distances and compute gaps
    distances = linkage_matrix[:, 2]
    gaps = np.diff(distances)

    if len(gaps) == 0:
        return np.zeros(n_markers, dtype=int)

    # For small lineages (3-5 markers), be very conservative about cutting
    # These lineages already passed FDR, top-k, GMM threshold, and gap-based splitting
    # Only cut if there's an extremely clear gap (>= 3x median)
    if n_markers <= 5:
        max_gap = np.max(gaps)
        median_gap = np.median(gaps)

        # Require gap to be at least 3x median for small lineages
        if median_gap <= 0 or max_gap / median_gap < 3.0:
            # No clear internal structure - keep as single cluster
            return np.zeros(n_markers, dtype=int)

        # Found a clear gap - cut there
        max_gap_idx = np.argmax(gaps)
        cut_distance = distances[max_gap_idx + 1] - 1e-10

        try:
            clusters = fcluster(linkage_matrix, t=cut_distance, criterion='distance')
            # Only accept if we don't create singletons
            cluster_sizes = np.bincount(clusters)
            if np.min(cluster_sizes) >= 2:
                return clusters
            else:
                # Would create singletons - keep as single cluster
                return np.zeros(n_markers, dtype=int)
        except Exception:
            return np.zeros(n_markers, dtype=int)

    # For larger lineages (6+ markers), use adaptive GMM-based method
    n_components, threshold, is_large = _find_gap_threshold_gmm(gaps)

    if n_components > 1 and np.any(is_large):
        # Cut at the largest gap
        large_gap_indices = np.where(is_large)[0]
        gap_values = gaps[large_gap_indices]
        primary_cut_idx = large_gap_indices[np.argmax(gap_values)]

        # Cut distance is just before this merge
        cut_distance = distances[primary_cut_idx + 1] - 1e-10

        try:
            clusters = fcluster(linkage_matrix, t=cut_distance, criterion='distance')
            return clusters
        except Exception:
            pass

    # Fallback for large lineages: use inconsistency-based method
    # But penalize creating many small clusters
    try:
        incons = inconsistent(linkage_matrix, d=2)
        best_t = 1.0
        best_score = -np.inf

        for t in np.linspace(0.5, 2.0, 16):
            try:
                clusters = fcluster(linkage_matrix, t=t, criterion='inconsistent', depth=2)
                n_clusters = len(set(clusters))
                cluster_sizes = np.bincount(clusters)
                min_cluster_size = np.min(cluster_sizes)

                if n_clusters == 1:
                    score = 0.5  # Slight preference for keeping together
                elif n_clusters == n_markers:
                    score = 0  # Penalize all singletons
                elif min_cluster_size < 2:
                    score = 0.1  # Penalize creating singletons
                else:
                    # Prefer fewer, larger clusters
                    score = min_cluster_size / n_markers

                if score > best_score:
                    best_score = score
                    best_t = t
            except Exception:
                continue

        clusters = fcluster(linkage_matrix, t=best_t, criterion='inconsistent', depth=2)
        return clusters
    except Exception:
        # Ultimate fallback - keep as single cluster
        return np.zeros(n_markers, dtype=int)


def _compute_modularity(
    profiles: List[List[str]],
    pairs: List[MarkerPairColocalization],
    significant_pairs: List[MarkerPairColocalization],
) -> float:
    """
    Compute modularity: fraction of significant edges within profiles vs expected.

    Modularity Q = (within-profile edges / total) - (expected by random assignment)

    Args:
        profiles: List of marker sets.
        pairs: All colocalization pairs.
        significant_pairs: Pairs that passed significance filter.

    Returns:
        Modularity score in [-0.5, 1.0], higher is better.
    """
    if len(significant_pairs) == 0:
        return 0.0

    # Build profile lookup
    marker_to_profile = {}
    for i, profile in enumerate(profiles):
        for m in profile:
            marker_to_profile[m] = i

    # Count within-profile edges
    within = 0
    total_weight = 0

    for pair in significant_pairs:
        weight = pair.colocalization_score
        total_weight += weight

        p1 = marker_to_profile.get(pair.marker_a, -1)
        p2 = marker_to_profile.get(pair.marker_b, -1)

        if p1 == p2 and p1 != -1:
            within += weight

    if total_weight == 0:
        return 0.0

    # Expected within-profile edges if edges were random
    # Sum of (profile_size / total)^2 for each profile
    total_markers = len(marker_to_profile)
    if total_markers == 0:
        return 0.0

    expected = sum((len(p) / total_markers) ** 2 for p in profiles)

    # Modularity
    observed = within / total_weight
    modularity = observed - expected

    return float(modularity)


def discover_profiles_continuous(
    colocalization_result: ColocalizationResult,
    *,
    top_k: int = 5,
    distance_metric: str = "colocalization_score",
    seed: int = 1234,
    verbose: bool = True,
) -> ProfileDiscoveryResult:
    """
    Discover cell type profiles using continuous edge weighting (no FDR gate).

    Instead of binary FDR filtering, builds a complete distance matrix from all
    marker pairs and uses hierarchical clustering with gap-based cutting to
    separate lineages and profiles. This eliminates sensitivity to permutation
    count since no p-value thresholds are used.

    Algorithm:
    1. Build complete distance matrix from colocalization scores (all pairs)
    2. Optional: mutual top-k masking (set non-top-k distances to max)
    3. Hierarchical clustering (average linkage) over all markers
    4. Gap-based lineage splitting (large gaps = separate lineages)
    5. Dynamic tree cutting within each lineage to extract profiles
    6. Markers with uniformly high distance become singletons naturally

    Args:
        colocalization_result: Result from analyze_marker_colocalization().
        top_k: Number of top partners to keep per marker in mutual top-k
            distance masking (default: 5). Set to 0 to disable.
            Unlike FDR-based filtering, this only sets non-top-k distances
            to maximum (1.0) rather than removing edges entirely.
        distance_metric: Which metric to use for distance:
            - 'colocalization_score' (default): Combined score (0.3*spearman +
              0.3*cosine + 0.4*bivariate_morans), normalized to [0,1].
              Distance = 1 - score.
            - 'bivariate_morans': Uses bivariate Moran's I only.
              Distance = 1 - normalized_I.
        seed: Random seed for reproducibility (default: 1234).
        verbose: Log progress information (default: True).

    Returns:
        ProfileDiscoveryResult with discovered profiles, dendrograms, and metrics.
    """
    pairs = colocalization_result.pairs
    all_markers = colocalization_result.marker_names
    marker_specificity = colocalization_result.marker_specificity
    n_markers = len(all_markers)

    if verbose:
        logging.info(f"Discovering profiles (continuous weighting) from {len(pairs)} marker pairs...")
        logging.info(f"Distance metric: {distance_metric}, top_k: {top_k}")

    if n_markers < 2:
        return ProfileDiscoveryResult(
            profiles=[[m] for m in all_markers],
            lineage_dendrograms={},
            singletons=list(all_markers),
            modularity=0.0,
            n_significant_edges=0,
            alpha=0.0,
        )

    # Step 1: Build complete distance matrix from ALL pairs
    marker_to_idx = {m: i for i, m in enumerate(all_markers)}

    if distance_metric == "bivariate_morans":
        # Use bivariate Moran's I with min-max normalization
        morans_values = [p.bivariate_morans_i for p in pairs]
        if len(morans_values) > 0:
            min_i = min(morans_values)
            max_i = max(morans_values)
            range_i = max_i - min_i if max_i > min_i else 1.0
        else:
            min_i, range_i = 0.0, 1.0

        dist_matrix = np.ones((n_markers, n_markers))
        np.fill_diagonal(dist_matrix, 0.0)

        for pair in pairs:
            i = marker_to_idx[pair.marker_a]
            j = marker_to_idx[pair.marker_b]
            normalized = (pair.bivariate_morans_i - min_i) / range_i
            d = 1.0 - normalized
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d
    else:
        # Use colocalization_score (already in [0, 1])
        dist_matrix = np.ones((n_markers, n_markers))
        np.fill_diagonal(dist_matrix, 0.0)

        for pair in pairs:
            i = marker_to_idx[pair.marker_a]
            j = marker_to_idx[pair.marker_b]
            d = 1.0 - pair.colocalization_score
            dist_matrix[i, j] = d
            dist_matrix[j, i] = d

    if verbose:
        upper_tri = dist_matrix[np.triu_indices(n_markers, k=1)]
        logging.info(
            f"Distance matrix: min={upper_tri.min():.3f}, "
            f"median={np.median(upper_tri):.3f}, max={upper_tri.max():.3f}"
        )

    # Step 2: Optional mutual top-k distance masking
    if top_k > 0 and len(pairs) > 0:
        top_k_mask = np.zeros((n_markers, n_markers), dtype=bool)

        for i in range(n_markers):
            dists = dist_matrix[i, :].copy()
            dists[i] = np.inf  # Exclude self

            # Apply specificity weighting if available
            if marker_specificity:
                marker_name = all_markers[i]
                s_i = marker_specificity.get(marker_name, 1.0)
                for j in range(n_markers):
                    if j != i:
                        s_j = marker_specificity.get(all_markers[j], 1.0)
                        dists[j] = dists[j] / (np.sqrt(s_i * s_j) + 1e-10)

            k_actual = min(top_k, n_markers - 1)
            top_k_indices = np.argsort(dists)[:k_actual]
            for j in top_k_indices:
                top_k_mask[i, j] = True

        mutual_mask = top_k_mask & top_k_mask.T

        n_masked = 0
        for i in range(n_markers):
            for j in range(i + 1, n_markers):
                if not mutual_mask[i, j]:
                    dist_matrix[i, j] = 1.0
                    dist_matrix[j, i] = 1.0
                    n_masked += 1

        n_total_pairs = n_markers * (n_markers - 1) // 2
        n_kept = n_total_pairs - n_masked
        if verbose:
            logging.info(
                f"Mutual top-{top_k} masking: {n_kept}/{n_total_pairs} pairs kept "
                f"(specificity-weighted: {marker_specificity is not None})"
            )

    # Step 3: Hierarchical clustering over all markers
    condensed_dist = squareform(dist_matrix)
    Z = linkage(condensed_dist, method='average')

    if verbose:
        merge_distances = Z[:, 2]
        logging.info(
            f"Linkage: {len(Z)} merges, "
            f"distances [{merge_distances.min():.3f} ... {merge_distances.max():.3f}]"
        )

    # Step 4: Gap-based lineage splitting
    lineages = _split_dendrogram_by_gaps(Z, all_markers, seed=seed)

    if verbose:
        logging.info(f"Gap analysis: {len(lineages)} lineages")
        for i, lin in enumerate(lineages):
            logging.info(f"  Lineage {i + 1} ({len(lin)} markers): {lin}")

    # Step 5: Dynamic tree cutting within each lineage
    profiles = []
    lineage_dendrograms = {}
    lineage_idx = 0

    for lineage_markers in lineages:
        n_lineage = len(lineage_markers)

        if n_lineage == 1:
            profiles.append(lineage_markers)
            continue

        if n_lineage == 2:
            i = marker_to_idx[lineage_markers[0]]
            j = marker_to_idx[lineage_markers[1]]
            if dist_matrix[i, j] < 0.7:
                profiles.append(lineage_markers)
            else:
                profiles.append([lineage_markers[0]])
                profiles.append([lineage_markers[1]])
            continue

        # Build lineage-specific dendrogram
        lin_marker_indices = [marker_to_idx[m] for m in lineage_markers]
        lin_dist = dist_matrix[np.ix_(lin_marker_indices, lin_marker_indices)]
        lin_condensed = squareform(lin_dist)
        lin_Z = linkage(lin_condensed, method='average')

        lineage_dendrograms[lineage_idx] = LineageDendrogram(
            markers=lineage_markers,
            linkage_matrix=lin_Z,
            distance_matrix=lin_dist,
        )
        lineage_idx += 1

        cluster_labels = _dynamic_tree_cut(lin_Z, n_lineage)

        cluster_to_markers = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            cluster_to_markers[label].append(lineage_markers[i])

        for cluster_markers in cluster_to_markers.values():
            profiles.append(cluster_markers)

    # Identify singletons
    singletons = [p[0] for p in profiles if len(p) == 1]

    # Compute modularity
    significant_pairs = [p for p in pairs if p.colocalization_score > 0.5]
    modularity = _compute_modularity(profiles, pairs, significant_pairs)

    n_meaningful_edges = sum(1 for p in pairs if p.colocalization_score > 0.5)

    if verbose:
        logging.info(f"Discovered {len(profiles)} profiles, modularity = {modularity:.3f}")

    result = ProfileDiscoveryResult(
        profiles=profiles,
        lineage_dendrograms=lineage_dendrograms,
        singletons=singletons,
        modularity=modularity,
        n_significant_edges=n_meaningful_edges,
        alpha=0.0,
    )

    if verbose:
        logging.info(result.summary())

    return result


def _compute_profile_activity_map(
    X: NDArray[np.floating],
    marker_names: List[str],
    profile: List[str],
) -> NDArray[np.floating]:
    """
    Compute activity map for a profile (mean expression across profile markers).

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        profile: List of marker names in this profile

    Returns:
        Activity map (n_spots,)
    """
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}
    indices = [marker_to_idx[m] for m in profile if m in marker_to_idx]

    if not indices:
        return np.zeros(X.shape[0])

    return X[:, indices].mean(axis=1)


# =============================================================================
# MODULE 2c: PROFILE SELECTION (SPATIAL VARIANCE-BASED)
# =============================================================================


@dataclass
class ProfileSelectionResult:
    """Results from Module 2c profile selection."""

    selected_profiles: List[List[str]]  # Ordered list of selected profiles
    optimal_n: int  # Number of profiles selected
    variance_explained: NDArray[np.floating]  # Cumulative variance explained at each step
    proportion_smoothness: NDArray[np.floating]  # Smoothness of proportion estimates
    stopping_reason: str  # Why selection stopped
    all_profiles: List[List[str]]  # All candidate profiles (for reference)
    marginal_gains: NDArray[np.floating]  # Marginal variance gain at each step

    def summary(self) -> str:
        """Return a summary string."""
        total_ve = self.variance_explained[-1] if len(self.variance_explained) > 0 else 0.0
        return (
            f"Selected {self.optimal_n} of {len(self.all_profiles)} profiles\n"
            f"Total spatial variance explained: {total_ve:.1%}\n"
            f"Stopping reason: {self.stopping_reason}"
        )


def _compute_protein_variance_explained(
    X: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
) -> float:
    """
    Compute fraction of protein signal variance explained by profile activities.

    This measures how much of the total variance in protein expression data
    is captured by the markers included in the selected profiles.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        profiles: List of profiles (each is a list of marker names)

    Returns:
        Fraction of total variance explained (0 to 1)
    """
    if len(profiles) == 0:
        return 0.0

    # Get all markers covered by profiles
    covered_markers = set()
    for profile in profiles:
        covered_markers.update(profile)

    # Find column indices for covered markers
    marker_to_idx = {m: i for i, m in enumerate(marker_names)}
    covered_idx = [marker_to_idx[m] for m in covered_markers if m in marker_to_idx]

    if len(covered_idx) == 0:
        return 0.0

    # Compute total variance across all markers
    total_var = np.var(X, axis=0).sum()
    if total_var == 0:
        return 1.0

    # Compute variance explained by covered markers
    covered_var = np.var(X[:, covered_idx], axis=0).sum()

    return float(covered_var / total_var)


def _compute_spatial_variance_explained(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    neighbor_k: int = 8,
) -> float:
    """
    Compute fraction of spatial variance explained by profile activities.

    Uses Moran's I weighted sum: sum of (variance_i * |Moran_i|) across profiles,
    normalized by total spatial variance in the data.

    Args:
        X: Expression matrix (n_spots, n_markers)
        coords: Spatial coordinates (n_spots, 2)
        marker_names: Names for each column
        profiles: List of profiles (each is a list of marker names)
        neighbor_k: K for spatial neighbor graph

    Returns:
        Fraction of spatial variance explained (0-1)
    """
    if len(profiles) == 0:
        return 0.0

    try:
        from libpysal.weights import KNN as LibPySAL_KNN
        from esda.moran import Moran
    except ImportError:
        logger.warning("libpysal/esda not available - returning 0")
        return 0.0

    # Build spatial weights
    try:
        W = LibPySAL_KNN.from_array(coords, k=neighbor_k)
        W.transform = 'r'  # Row-standardize
    except Exception as e:
        logger.warning(f"Could not build spatial weights: {e}")
        return 0.0

    # Compute total spatial variance in the data (baseline)
    total_spatial_var = 0.0
    for j in range(X.shape[1]):
        col = X[:, j]
        if np.var(col) < 1e-10:
            continue
        try:
            mi = Moran(col, W)
            total_spatial_var += np.var(col) * abs(mi.I)
        except Exception:
            continue

    if total_spatial_var < 1e-10:
        return 1.0  # No spatial variance to explain

    # Compute spatial variance explained by profiles
    explained_var = 0.0
    for profile in profiles:
        activity = _compute_profile_activity_map(X, marker_names, profile)
        if np.var(activity) < 1e-10:
            continue
        try:
            mi = Moran(activity, W)
            explained_var += np.var(activity) * abs(mi.I)
        except Exception:
            continue

    return min(1.0, explained_var / total_spatial_var)


def _compute_proportion_smoothness(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    neighbor_k: int = 8,
) -> float:
    """
    Compute smoothness of profile-based proportion estimates.

    Uses average Moran's I of profile activities as a proxy for how
    spatially coherent the proportions would be.

    Returns:
        Mean Moran's I across profile activities (higher = smoother)
    """
    if len(profiles) == 0:
        return 0.0

    try:
        from libpysal.weights import KNN as LibPySAL_KNN
        from esda.moran import Moran
    except ImportError:
        return 0.0

    try:
        W = LibPySAL_KNN.from_array(coords, k=neighbor_k)
        W.transform = 'r'
    except Exception:
        return 0.0

    morans = []
    for profile in profiles:
        activity = _compute_profile_activity_map(X, marker_names, profile)
        if np.var(activity) < 1e-10:
            continue
        try:
            mi = Moran(activity, W)
            morans.append(mi.I)
        except Exception:
            continue

    return float(np.mean(morans)) if morans else 0.0


def _score_profile_spatial_contribution(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    profile: List[str],
    selected_profiles: List[List[str]],
    neighbor_k: int = 8,
) -> float:
    """
    Score a candidate profile by its marginal spatial variance contribution.

    Returns the increase in spatial variance explained when adding this profile.
    """
    # Current variance explained
    current_ve = _compute_spatial_variance_explained(
        X, coords, marker_names, selected_profiles, neighbor_k
    )

    # Variance explained with new profile
    new_profiles = selected_profiles + [profile]
    new_ve = _compute_spatial_variance_explained(
        X, coords, marker_names, new_profiles, neighbor_k
    )

    return new_ve - current_ve


def rescue_singletons(
    profiles: List[List[str]],
    signal_masks: NDArray[np.bool_],
    signal_mask_marker_names: List[str],
    min_unique_coverage: float = 0.3,
    min_signal_fraction: float = 0.05,
    verbose: bool = False,
) -> List[List[str]]:
    """
    Filter singletons by unique spatial coverage.

    Keeps singletons that cover spatial territory not explained by multi-marker
    profiles. Drops noise singletons that overlap heavily with multi-marker
    profile territory.

    Algorithm:
    1. Separate profiles into singletons (1 marker) and multi-marker (2+)
    2. Compute "explained territory" = spots where at least 2 markers from
       the same multi-marker profile co-express (both in signal component).
       This prevents broadly-expressed markers from inflating territory.
    3. For each singleton, compute:
       - unique_coverage: fraction of its signal spots outside explained territory
       - signal_fraction: fraction of all spots classified as signal
    4. Keep singleton if unique_coverage >= min_unique_coverage AND
       signal_fraction >= min_signal_fraction

    Args:
        profiles: List of profiles (each a list of marker names).
        signal_masks: Boolean array (n_spots, n_markers) from Module 1 GMM.
            True = spot is in signal component for that marker.
        signal_mask_marker_names: Marker names corresponding to columns of signal_masks.
        min_unique_coverage: Minimum fraction of signal spots that must be
            outside explained territory (default: 0.3).
        min_signal_fraction: Minimum fraction of all spots that must be
            classified as signal (default: 0.05).
        verbose: Log rescue decisions (default: False).

    Returns:
        Filtered list of profiles with noise singletons removed.
    """
    n_spots = signal_masks.shape[0]

    # Build marker name -> column index mapping
    marker_to_col = {name: i for i, name in enumerate(signal_mask_marker_names)}

    # Separate singletons from multi-marker profiles
    singletons = []
    multi_marker = []
    for profile in profiles:
        if len(profile) == 1:
            singletons.append(profile)
        else:
            multi_marker.append(profile)

    if not singletons:
        if verbose:
            logging.info("No singletons to rescue/filter")
        return profiles

    # Compute explained territory: spots with co-expression in multi-marker profiles
    # A spot is "explained" if at least 2 markers from the SAME profile have signal there.
    # This avoids broadly-expressed markers (CD45, CD3E) inflating the territory.
    explained = np.zeros(n_spots, dtype=bool)
    for profile in multi_marker:
        # Count how many markers in this profile have signal per spot
        profile_cols = [marker_to_col[m] for m in profile if m in marker_to_col]
        if len(profile_cols) < 2:
            # Single-marker match in multi-marker profile: use that marker's signal
            for col in profile_cols:
                explained |= signal_masks[:, col]
        else:
            # Co-expression: spots where >= 2 markers from this profile have signal
            coexpr_count = np.zeros(n_spots, dtype=int)
            for col in profile_cols:
                coexpr_count += signal_masks[:, col].astype(int)
            explained |= (coexpr_count >= 2)

    n_explained = explained.sum()
    if verbose:
        logging.info(
            f"Explained territory (co-expression): {n_explained}/{n_spots} spots "
            f"({n_explained / n_spots:.1%}) covered by {len(multi_marker)} multi-marker profiles"
        )

    # Evaluate each singleton
    rescued = []
    dropped = []
    for profile in singletons:
        marker = profile[0]
        if marker not in marker_to_col:
            if verbose:
                logging.info(f"  Singleton [{marker}]: not in signal masks, dropping")
            dropped.append(marker)
            continue

        col = marker_to_col[marker]
        signal_spots = signal_masks[:, col]
        n_signal = signal_spots.sum()
        signal_fraction = n_signal / n_spots

        if n_signal == 0:
            if verbose:
                logging.info(f"  Singleton [{marker}]: no signal spots, dropping")
            dropped.append(marker)
            continue

        # Unique coverage: signal spots NOT in explained territory
        unique_spots = signal_spots & ~explained
        unique_coverage = unique_spots.sum() / n_signal

        keep = unique_coverage >= min_unique_coverage and signal_fraction >= min_signal_fraction

        if verbose:
            status = "KEEP" if keep else "DROP"
            logging.info(
                f"  Singleton [{marker}]: signal_fraction={signal_fraction:.3f}, "
                f"unique_coverage={unique_coverage:.3f} -> {status}"
            )

        if keep:
            rescued.append(profile)
        else:
            dropped.append(marker)

    result = multi_marker + rescued

    if verbose:
        logging.info(
            f"Singleton rescue: {len(rescued)} kept, {len(dropped)} dropped "
            f"(total profiles: {len(result)})"
        )

    return result


def select_profiles(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    interesting_markers: List[str],
    colocalization_result: "ColocalizationResult",
    min_spatial_explained: float = 0.90,
    min_protein_explained: float = 0.90,
    max_profiles: int = 15,
    min_profiles: int = 5,
    neighbor_k: int = 8,
    verbose: bool = False,
) -> ProfileSelectionResult:
    """
    Module 2c: Select profiles by dual variance contribution.

    Greedily selects profiles that maximize spatial variance explained,
    stopping when BOTH variance targets are reached.

    Uses dual variance checkpoints:
    1. Protein signal variance - fraction of total marker variance covered
    2. Spatial coverage - Moran's I weighted spatial variance explained

    Selection continues by rank (best contributing profile first) until
    BOTH targets reach 90% or all profiles are exhausted.

    Args:
        X: Expression matrix (n_spots, n_markers), should be RAW (not CLR-transformed)
        coords: Spatial coordinates (n_spots, 2)
        marker_names: Names for each column
        profiles: Candidate profiles from Module 2b (discover_profiles)
        interesting_markers: Markers identified as spatially interesting (Module 1)
        colocalization_result: Result from Module 2a (for reference)
        min_spatial_explained: Target fraction of spatial variance (default: 0.90)
        min_protein_explained: Target fraction of protein signal variance (default: 0.90)
        max_profiles: Maximum number of profiles to select (default: 15)
        min_profiles: Minimum profiles to select before checking marginal gain (default: 5)
        neighbor_k: K for spatial neighbor graph (default: 8)
        verbose: Whether to print progress

    Returns:
        ProfileSelectionResult with selected profiles and metrics
    """
    if verbose:
        logger.info(f"Module 2c: Selecting from {len(profiles)} candidate profiles")
        logger.info(f"  Target spatial variance: {min_spatial_explained:.0%}")
        logger.info(f"  Target protein variance: {min_protein_explained:.0%}")
        logger.info(f"  Selection: by rank until BOTH targets met")

    if len(profiles) == 0:
        return ProfileSelectionResult(
            selected_profiles=[],
            optimal_n=0,
            variance_explained=np.array([]),
            proportion_smoothness=np.array([]),
            stopping_reason="No candidate profiles",
            all_profiles=[],
            marginal_gains=np.array([]),
        )

    selected: List[List[str]] = []
    remaining = [list(p) if isinstance(p, (set, frozenset)) else p for p in profiles]
    variance_explained_list = []
    smoothness_list = []
    marginal_gains_list = []
    stopping_reason = "max_profiles reached"

    for i in range(min(max_profiles, len(remaining))):
        # Score all remaining profiles by marginal contribution
        best_profile = None
        best_gain = -np.inf
        best_idx = -1

        for idx, profile in enumerate(remaining):
            gain = _score_profile_spatial_contribution(
                X, coords, marker_names, profile, selected, neighbor_k
            )
            if gain > best_gain:
                best_gain = gain
                best_profile = profile
                best_idx = idx

        # Check stopping conditions
        if best_profile is None:
            stopping_reason = "No remaining profiles"
            break

        # Add best profile (no marginal gain check - select until variance targets met)
        selected.append(best_profile)
        remaining.pop(best_idx)

        # Compute metrics - dual variance checkpoints
        current_spatial_ve = _compute_spatial_variance_explained(
            X, coords, marker_names, selected, neighbor_k
        )
        current_protein_ve = _compute_protein_variance_explained(
            X, marker_names, selected
        )
        current_smooth = _compute_proportion_smoothness(
            X, coords, marker_names, selected, neighbor_k
        )

        variance_explained_list.append(current_spatial_ve)
        smoothness_list.append(current_smooth)
        marginal_gains_list.append(best_gain)

        if verbose:
            print(f"  [{i+1}] Added {best_profile} (gain: {best_gain:.3%}, "
                  f"spatial: {current_spatial_ve:.1%}, protein: {current_protein_ve:.1%})")

        # Check if BOTH variance targets reached (dual checkpoint)
        spatial_met = current_spatial_ve >= min_spatial_explained
        protein_met = current_protein_ve >= min_protein_explained

        if spatial_met and protein_met:
            stopping_reason = (f"Both targets reached "
                             f"(spatial: {current_spatial_ve:.0%}, protein: {current_protein_ve:.0%})")
            break

    return ProfileSelectionResult(
        selected_profiles=selected,
        optimal_n=len(selected),
        variance_explained=np.array(variance_explained_list),
        proportion_smoothness=np.array(smoothness_list),
        stopping_reason=stopping_reason,
        all_profiles=[list(p) if isinstance(p, (set, frozenset)) else p for p in profiles],
        marginal_gains=np.array(marginal_gains_list),
    )


# =============================================================================
# HIERARCHICAL NMF-BASED PROFILE DISCOVERY (Main Entry Point)
# =============================================================================


def _detect_adaptive_threshold(
    topk_pairs: List[MarkerPairColocalization],
    fdr_pairs: List[MarkerPairColocalization],
    fallback: float = 0.85,
    verbose: bool = True,
) -> float:
    """
    Automatically detect the optimal adaptive ratio threshold for this sample.

    Approach: Compute edge ratios (score / min(best_a, best_b)) and find the
    largest gap in the distribution. Within-celltype edges cluster near 1.0,
    cross-celltype bridges have lower ratios. The gap separates them.

    Args:
        topk_pairs: Pairs after top-k filtering
        fdr_pairs: Pairs after FDR filtering (for computing best scores)
        fallback: Fallback threshold if no clear gap is found
        verbose: Log detection details

    Returns:
        Detected threshold ratio (between 0.75 and 0.95)
    """
    if len(topk_pairs) < 3:
        return fallback

    # Build per-marker best scores
    marker_best: Dict[str, float] = defaultdict(float)
    for p in fdr_pairs:
        marker_best[p.marker_a] = max(marker_best[p.marker_a], p.colocalization_score)
        marker_best[p.marker_b] = max(marker_best[p.marker_b], p.colocalization_score)

    # Compute ratio for each edge
    edge_ratios = []
    for p in topk_pairs:
        best_a = marker_best.get(p.marker_a, 0)
        best_b = marker_best.get(p.marker_b, 0)
        if min(best_a, best_b) > 0:
            ratio = p.colocalization_score / min(best_a, best_b)
            edge_ratios.append(ratio)

    if len(edge_ratios) < 3:
        return fallback

    edge_ratios = np.array(sorted(edge_ratios))

    # Find largest gap in the distribution
    # Only consider gaps in the range [0.75, 0.95] - outside this range is unlikely to be useful
    gaps = []
    for i in range(len(edge_ratios) - 1):
        lower = edge_ratios[i]
        upper = edge_ratios[i + 1]
        midpoint = (lower + upper) / 2
        gap_size = upper - lower

        # Only consider gaps where the midpoint is in [0.75, 0.95]
        if 0.75 <= midpoint <= 0.95 and gap_size > 0.02:  # Min gap size of 2%
            gaps.append((gap_size, midpoint, lower, upper))

    if gaps:
        # Use the largest gap
        gaps.sort(reverse=True)
        best_gap = gaps[0]
        detected_threshold = best_gap[1]  # Use midpoint of gap

        if verbose:
            logger.info(f"Auto-detected threshold: {detected_threshold:.3f} "
                       f"(gap of {best_gap[0]:.3f} between {best_gap[2]:.3f} and {best_gap[3]:.3f})")

        return detected_threshold

    # No clear gap found - use elbow detection
    # Find where removing edges causes the biggest change in edge count
    # This is a simpler heuristic: use the ratio at the 25th percentile
    # (keeping top 75% of edges by ratio)
    if len(edge_ratios) >= 4:
        percentile_threshold = np.percentile(edge_ratios, 25)
        if 0.75 <= percentile_threshold <= 0.95:
            if verbose:
                logger.info(f"Auto-detected threshold (25th percentile): {percentile_threshold:.3f}")
            return percentile_threshold

    if verbose:
        logger.info(f"No clear gap found, using fallback threshold: {fallback}")
    return fallback


def discover_hierarchical_profiles_continuous(
    coloc_result: ColocalizationResult,
    antibody_expression: NDArray[np.floating],
    marker_names: List[str],
    coords: NDArray[np.floating],
    improvement_threshold: float = 0.05,
    sharing_ratio: float = 0.5,
    sharing_min_I: float = 0.2,
    max_depth: int = 5,
    neighbor_k: int = 6,
    top_k: int = 3,
    distance_metric: str = "colocalization_score",
    verbose: bool = True,
) -> HierarchicalProfileResult:
    """
    Discover hierarchical profiles using flat-first, hierarchy-second approach.

    Stage 1: Runs discover_profiles_continuous() to get high-quality flat
    profiles (proven 86-100% GT coverage on Xenium benchmarks).

    Stage 2: Builds a post-hoc hierarchy over the flat profiles by clustering
    them by Jaccard similarity of marker sets. The tree is purely
    organizational — it does not modify the profiles.

    The function signature is unchanged from the tree-first version so that
    callers (run_benchmark.py) require no modification.

    Args:
        coloc_result: Result from analyze_marker_colocalization()
        antibody_expression: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column in expression matrix
        coords: Spatial coordinates (n_spots, 2) [unused, kept for API compat]
        improvement_threshold: Unused, kept for API compatibility.
        sharing_ratio: Unused, kept for API compatibility.
        sharing_min_I: Unused, kept for API compatibility.
        max_depth: Unused, kept for API compatibility.
        neighbor_k: Unused, kept for API compatibility.
        top_k: Number of top partners for mutual top-k masking (default 3).
        distance_metric: Which metric for distance (default 'colocalization_score').
        verbose: Log progress (default True)

    Returns:
        HierarchicalProfileResult with flat profiles and post-hoc tree.
    """
    if verbose:
        logger.info("Hierarchical profile discovery (flat-first, hierarchy post-hoc)...")

    # ── Stage 1: Flat profile discovery ───────────────────────────────────
    flat_result = discover_profiles_continuous(
        colocalization_result=coloc_result,
        top_k=top_k,
        distance_metric=distance_metric,
        verbose=verbose,
    )

    # Convert to dict format expected by HierarchicalProfileResult
    flat_profiles: Dict[str, List[str]] = {}
    for i, profile_markers in enumerate(flat_result.profiles):
        flat_profiles[f"profile_{i}"] = list(profile_markers)
    for i, singleton in enumerate(flat_result.singletons):
        flat_profiles[f"singleton_{i}"] = [singleton]

    if verbose:
        logger.info(f"Stage 1: {len(flat_result.profiles)} multi-marker profiles, "
                     f"{len(flat_result.singletons)} singletons")

    # ── Stage 2: Post-hoc hierarchy ───────────────────────────────────────
    tree, shared_markers = _build_posthoc_hierarchy(flat_profiles)

    if verbose:
        logger.info(f"Stage 2: tree depth {tree.n_levels}, "
                     f"{len(shared_markers)} shared markers")

    # Compute reconstruction error
    final_error = _compute_final_reconstruction_error(
        antibody_expression, marker_names, list(flat_profiles.values())
    )

    if verbose:
        logger.info(f"Discovered {len(flat_profiles)} profiles, "
                     f"reconstruction error: {final_error:.4f}")

    return HierarchicalProfileResult(
        tree=tree,
        flat_profiles=flat_profiles,
        depth_per_branch={},
        shared_markers=shared_markers,
        reconstruction_error=final_error,
    )


def _filter_coloc_result_for_markers(
    coloc_result: ColocalizationResult,
    markers: List[str],
) -> ColocalizationResult:
    """Filter colocalization result to only include specified markers."""
    marker_set = set(markers)
    filtered_pairs = [
        p for p in coloc_result.pairs
        if p.marker_a in marker_set and p.marker_b in marker_set
    ]
    return ColocalizationResult(
        pairs=filtered_pairs,
        marker_names=markers,
        n_spots=coloc_result.n_spots,
        neighbor_k=coloc_result.neighbor_k,
        marker_specificity={m: coloc_result.marker_specificity.get(m, 0.5)
                          for m in markers} if coloc_result.marker_specificity else None,
    )


def _compute_all_nmf_weights(
    node: ProfileTreeNode,
    X: NDArray[np.floating],
    marker_names: List[str],
    all_weights: Dict[str, Dict[str, Dict[str, float]]],
) -> None:
    """Recursively compute NMF weights for all internal nodes."""
    if node.is_leaf:
        return

    # Compute weights for this node
    weights = _compute_nmf_weights(X, marker_names, node)
    all_weights[node.node_id] = weights

    # Recurse on children
    for child in node.children:
        _compute_all_nmf_weights(child, X, marker_names, all_weights)


def _compute_final_reconstruction_error(
    X: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
) -> float:
    """Compute overall reconstruction error for final profiles."""
    from scipy.optimize import nnls

    if len(profiles) == 0:
        return 1.0

    marker_to_idx = {m: i for i, m in enumerate(marker_names)}

    # Build profile matrix: each column is mean expression of profile markers
    profile_vectors = []
    for profile in profiles:
        indices = [marker_to_idx[m] for m in profile if m in marker_to_idx]
        if len(indices) > 0:
            profile_expr = X[:, indices].mean(axis=1)
            profile_vectors.append(profile_expr)

    if len(profile_vectors) == 0:
        return 1.0

    P = np.column_stack(profile_vectors)

    # For each marker, compute reconstruction error
    total_error = 0.0
    all_markers = set(m for p in profiles for m in p)

    for marker in all_markers:
        if marker not in marker_to_idx:
            continue
        y = X[:, marker_to_idx[marker]]
        try:
            coeffs, _ = nnls(P, y)
            y_pred = P @ coeffs
            error = np.mean((y - y_pred) ** 2)
            total_error += error
        except Exception:
            total_error += 1.0

    return total_error / max(len(all_markers), 1)


def _build_posthoc_hierarchy(
    flat_profiles: Dict[str, List[str]],
) -> Tuple[ProfileTree, Dict[str, List[str]]]:
    """
    Build a post-hoc hierarchy over pre-computed flat profiles.

    Clusters profiles by Jaccard similarity of their marker sets, producing
    a ProfileTree for visualization/reporting.  The flat profiles themselves
    are NOT modified — this is purely organizational.

    Args:
        flat_profiles: Dict mapping profile_name -> list of marker names.

    Returns:
        Tuple of (ProfileTree, shared_markers_dict).
        shared_markers_dict maps marker -> list of profile names that contain it.
    """
    from scipy.cluster.hierarchy import linkage, to_tree
    from scipy.spatial.distance import squareform

    profile_names = list(flat_profiles.keys())
    profile_sets = [set(markers) for markers in flat_profiles.values()]
    n = len(profile_names)

    # Identify shared markers (appear in 2+ profiles)
    marker_to_profiles: Dict[str, List[str]] = defaultdict(list)
    for name, markers in flat_profiles.items():
        for m in markers:
            marker_to_profiles[m].append(name)
    shared_markers = {m: ps for m, ps in marker_to_profiles.items() if len(ps) > 1}

    # Trivial cases
    if n <= 1:
        if n == 1:
            node = ProfileTreeNode(
                node_id=profile_names[0],
                markers=flat_profiles[profile_names[0]],
                children=[], parent_id=None, depth=0,
            )
        else:
            node = ProfileTreeNode(
                node_id="root", markers=[], children=[],
                parent_id=None, depth=0,
            )
        return ProfileTree(root=node, n_levels=1), shared_markers

    # Pairwise Jaccard distance
    D = np.ones((n, n))
    np.fill_diagonal(D, 0.0)
    for i in range(n):
        for j in range(i + 1, n):
            intersection = len(profile_sets[i] & profile_sets[j])
            union = len(profile_sets[i] | profile_sets[j])
            jaccard = intersection / union if union > 0 else 0.0
            D[i, j] = 1.0 - jaccard
            D[j, i] = 1.0 - jaccard

    # Hierarchical clustering
    condensed = squareform(D)
    Z = linkage(condensed, method='average')
    scipy_root = to_tree(Z)

    # Convert scipy tree to ProfileTree
    def _convert(node, depth=0, parent_id=None):
        if node.is_leaf():
            idx = node.id
            return ProfileTreeNode(
                node_id=profile_names[idx],
                markers=flat_profiles[profile_names[idx]],
                children=[], parent_id=parent_id, depth=depth,
            )
        node_id = f"internal_{node.id}"
        left = _convert(node.get_left(), depth + 1, node_id)
        right = _convert(node.get_right(), depth + 1, node_id)
        return ProfileTreeNode(
            node_id=node_id, markers=[], children=[left, right],
            parent_id=parent_id, depth=depth,
        )

    root = _convert(scipy_root)
    n_levels = _compute_tree_depth(root)
    return ProfileTree(root=root, n_levels=n_levels), shared_markers

