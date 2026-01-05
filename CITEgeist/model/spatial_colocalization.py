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
    bivariate_morans_i: float  # Spatial cross-correlation
    bivariate_morans_pvalue: float  # P-value from permutation test

    # Neighborhood analysis (binary-based, kept for backwards compatibility)
    neighbor_enrichment_ab: float  # Enrichment of B in neighbors of A-positive spots
    neighbor_enrichment_ba: float  # Enrichment of A in neighbors of B-positive spots
    neighbor_enrichment_pvalue: float
    mutual_neighbor_enrichment: float  # Mean of A->B and B->A enrichment

    # Combined score
    colocalization_score: float  # Weighted combination of all metrics


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
    # Bivariate Moran's I with permutation test - default for FDR correction
    # Uses smoothed data to reduce noise and improve detection in mixed tissue
    bivariate_i, bivariate_pval = _compute_bivariate_morans_i_pvalue(
        smooth_a, smooth_b, neighbors, rng, n_perm=n_permutations
    )

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
    n_permutations: int = 199,
    seed: int = 1234,
    verbose: bool = True,
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

    # Run in parallel using joblib
    pairs = Parallel(n_jobs=n_jobs, verbose=0)(
        delayed(_process_marker_pair)(
            i, j, analyze_names, analyze_X, analyze_X_smooth, binary,
            neighbors, n_permutations, seed
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


def _triangle_assembly(
    markers: List[str],
    edges: List[Tuple[str, str, float]],
    verbose: bool = True,
) -> List[List[str]]:
    """
    Assemble profiles using triangle-first strategy with two-support attachment.

    Algorithm:
    1. Find all triangles (3 nodes with all 3 edges) - strongest evidence
    2. Triangles become triplet profiles (assigned strongest-first by mean weight)
    3. Remaining edges become pair profiles
    4. TWO-SUPPORT ATTACHMENT: For each pair {a,b}, attach singletons that have
       edges to BOTH a and b (forming a triplet with structural evidence)
    5. Remaining nodes (no two-support) become singletons

    Args:
        markers: List of marker names
        edges: List of (marker_a, marker_b, weight) tuples
        verbose: Whether to log assembly details

    Returns:
        List of profiles (each is a list of marker names)
    """
    if len(edges) == 0:
        # No edges - all singletons
        return [[m] for m in markers]

    # Build adjacency set for fast lookup
    adj: Dict[str, Set[str]] = defaultdict(set)
    edge_weights: Dict[Tuple[str, str], float] = {}
    for a, b, w in edges:
        adj[a].add(b)
        adj[b].add(a)
        edge_key = tuple(sorted([a, b]))
        edge_weights[edge_key] = w

    # Step 1: Find all triangles
    triangles: List[Tuple[Tuple[str, str, str], float]] = []
    seen_triangles: Set[Tuple[str, str, str]] = set()

    for a in adj:
        for b in adj[a]:
            if b <= a:
                continue
            # Find common neighbors (potential triangle completions)
            common = adj[a] & adj[b]
            for c in common:
                if c <= b:
                    continue
                tri = tuple(sorted([a, b, c]))
                if tri not in seen_triangles:
                    # Compute mean edge weight for this triangle
                    w_ab = edge_weights[tuple(sorted([a, b]))]
                    w_ac = edge_weights[tuple(sorted([a, c]))]
                    w_bc = edge_weights[tuple(sorted([b, c]))]
                    mean_weight = (w_ab + w_ac + w_bc) / 3
                    triangles.append((tri, mean_weight))
                    seen_triangles.add(tri)

    # Sort triangles by mean weight (strongest first)
    triangles.sort(key=lambda x: x[1], reverse=True)

    if verbose and triangles:
        logging.info(f"  Found {len(triangles)} triangles in graph")

    # Step 2: Greedily assign triangles (no node reuse)
    used_nodes: Set[str] = set()
    profiles: List[List[str]] = []

    for tri, weight in triangles:
        if any(node in used_nodes for node in tri):
            continue  # Skip if any node already used
        profiles.append(list(tri))
        used_nodes.update(tri)
        if verbose:
            logging.info(f"    Triangle profile: {list(tri)} (mean_weight={weight:.3f})")

    # Step 3: Remaining edges become pairs (strongest first)
    sorted_edges = sorted(edges, key=lambda x: x[2], reverse=True)
    pair_profiles: List[Tuple[List[str], float]] = []  # Track pairs for two-support
    for a, b, w in sorted_edges:
        if a in used_nodes or b in used_nodes:
            continue
        pair_profiles.append(([a, b], w))
        used_nodes.update([a, b])
        if verbose:
            logging.info(f"    Pair profile: [{a}, {b}] (weight={w:.3f})")

    # Step 4: TWO-SUPPORT ATTACHMENT
    # For each pair, check if any remaining singleton has edges to BOTH members
    remaining_singletons = set(markers) - used_nodes
    attached_count = 0

    for pair, pair_weight in pair_profiles:
        a, b = pair
        # Find singletons with edges to both a and b
        candidates = []
        for u in remaining_singletons:
            has_edge_a = u in adj.get(a, set())
            has_edge_b = u in adj.get(b, set())
            if has_edge_a and has_edge_b:
                # Compute attachment strength (mean of u-a and u-b edge weights)
                w_ua = edge_weights.get(tuple(sorted([u, a])), 0)
                w_ub = edge_weights.get(tuple(sorted([u, b])), 0)
                attach_strength = (w_ua + w_ub) / 2
                candidates.append((u, attach_strength))

        # Attach the strongest candidate (if any)
        if candidates:
            candidates.sort(key=lambda x: x[1], reverse=True)
            best_u, best_strength = candidates[0]
            pair.append(best_u)
            remaining_singletons.remove(best_u)
            used_nodes.add(best_u)
            attached_count += 1
            if verbose:
                logging.info(
                    f"    Two-support attachment: {best_u} -> {pair[:2]} "
                    f"(strength={best_strength:.3f})"
                )

    # Add pairs (now possibly triplets) to profiles
    for pair, _ in pair_profiles:
        profiles.append(pair)

    if verbose and attached_count > 0:
        logging.info(f"  Two-support attachment: {attached_count} singletons promoted to triplets")

    # Step 5: Remaining nodes become singletons
    for marker in markers:
        if marker not in used_nodes:
            profiles.append([marker])

    return profiles


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

    # Find large gaps using adaptive method
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

    # Fallback: use inconsistency-based method
    try:
        incons = inconsistent(linkage_matrix, d=2)
        best_t = 1.0
        best_score = -np.inf

        for t in np.linspace(0.5, 2.0, 16):
            try:
                clusters = fcluster(linkage_matrix, t=t, criterion='inconsistent', depth=2)
                n_clusters = len(set(clusters))

                if n_clusters == 1:
                    score = 0
                elif n_clusters == n_markers:
                    score = 0
                else:
                    score = n_clusters * (n_markers - n_clusters) / n_markers

                if score > best_score:
                    best_score = score
                    best_t = t
            except Exception:
                continue

        clusters = fcluster(linkage_matrix, t=best_t, criterion='inconsistent', depth=2)
        return clusters
    except Exception:
        # Ultimate fallback
        return np.arange(n_markers)


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


def _fit_score_gmm_with_bic(
    scores: NDArray[np.floating],
    max_components: int = 5,
    seed: int = 1234,
) -> Tuple[int, float, NDArray[np.bool_]]:
    """
    Fit GMM to colocalization scores with BIC-based model selection.

    Tries 1 to max_components GMM fits and selects the best by BIC.
    Then identifies which edges are in the "high colocalization" component(s)
    vs the "noise/spatial proximity" component.

    Args:
        scores: Array of colocalization scores.
        max_components: Maximum number of GMM components to try.
        seed: Random seed.

    Returns:
        Tuple of:
        - n_components: Selected number of components
        - score_threshold: Threshold separating noise from signal
        - is_signal: Boolean array indicating which scores are signal (not noise)
    """
    from sklearn.mixture import GaussianMixture

    scores = scores[np.isfinite(scores)]
    if len(scores) < 10:
        # Not enough data, keep everything
        return 1, 0.0, np.ones(len(scores), dtype=bool)

    scores_2d = scores.reshape(-1, 1)

    # Fit GMMs with 1 to max_components
    best_bic = np.inf
    best_gmm = None
    best_k = 1

    for k in range(1, max_components + 1):
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type='full',
                random_state=seed,
                n_init=3,
                max_iter=200,
            )
            gmm.fit(scores_2d)
            bic = gmm.bic(scores_2d)

            if bic < best_bic:
                best_bic = bic
                best_gmm = gmm
                best_k = k
        except Exception:
            continue

    if best_gmm is None:
        return 1, 0.0, np.ones(len(scores), dtype=bool)

    # Identify components by mean
    means = best_gmm.means_.flatten()
    sorted_indices = np.argsort(means)

    # The lowest-mean component is "noise/spatial proximity"
    # All other components are "signal"
    noise_component = sorted_indices[0]

    # Get component assignments
    labels = best_gmm.predict(scores_2d)

    # Signal = not in the noise component
    is_signal = labels != noise_component

    # Threshold is the boundary between noise and lowest signal component
    if best_k == 1:
        # Only one component - use median as threshold
        threshold = np.median(scores)
        is_signal = scores > threshold
    else:
        # Threshold = midpoint between noise mean and next component mean
        noise_mean = means[noise_component]
        signal_means = [means[i] for i in sorted_indices[1:]]
        lowest_signal_mean = min(signal_means) if signal_means else noise_mean
        threshold = (noise_mean + lowest_signal_mean) / 2

    return best_k, float(threshold), is_signal


def discover_profiles(
    colocalization_result: ColocalizationResult,
    *,
    fdr_alpha: float = 0.05,
    top_k: int = 5,
    min_score: Optional[float] = None,
    use_triangle_assembly: bool = False,
    pvalue_source: str = "bivariate_morans",
    seed: int = 1234,
    verbose: bool = True,
) -> ProfileDiscoveryResult:
    """
    Discover cell type profiles from colocalization analysis.

    Algorithm:
    1. Apply BH-FDR correction to p-values (controls false discovery rate)
    2. Apply mutual top-k sparsification (prevents hub marker collapse)
    3. Optionally apply score threshold (GMM-adaptive or fixed)
    4. Either:
       a) Triangle-first assembly (if use_triangle_assembly=True):
          - Find all triangles (3 nodes with all 3 edges) as triplet profiles
          - Remaining edges become pair profiles
          - Remaining nodes become singletons
       b) Hierarchical clustering (default):
          - Find connected components
          - Within each component, perform hierarchical clustering
          - Dynamic tree cutting to extract profiles
    5. Compute modularity to measure how well profiles explain the data

    Args:
        colocalization_result: Result from analyze_marker_colocalization().
        fdr_alpha: FDR threshold for Benjamini-Hochberg correction (default: 0.05).
            Use stricter threshold (e.g., 0.01) to reduce false edges if needed.
        top_k: Number of top partners to keep per marker in mutual top-k
            sparsification (default: 3). Lower values = more aggressive pruning.
            Set to 0 to disable.
        min_score: Minimum colocalization score. If None (default), uses
            adaptive GMM-based threshold with BIC model selection.
        use_triangle_assembly: If True, use triangle-first assembly instead
            of hierarchical clustering (default: False). Triangle assembly
            is more robust for detecting triplets but doesn't build dendrograms.
        pvalue_source: Which p-value to use for FDR correction:
            - 'bivariate_morans' (default): Tests spatial cross-correlation
              between markers. Works well for both mixed and high-seg data.
            - 'neighbor_enrichment': Binary-based. Best for high-segregation
              data with clear cell type boundaries.
            - 'correlation': Pearson/Spearman p-value.
        seed: Random seed for reproducibility (default: 1234).
        verbose: Log progress information (default: True).

    Returns:
        ProfileDiscoveryResult with discovered profiles, dendrograms, and metrics.

    Example:
        >>> # Run full Module 2 pipeline
        >>> coloc_result = analyze_marker_colocalization(X, coords, markers)
        >>> profiles = discover_profiles(coloc_result)
        >>> print(profiles.summary())
        >>> for profile in profiles.profiles:
        ...     print(f"Profile: {profile}")
    """
    pairs = colocalization_result.pairs
    all_markers = colocalization_result.marker_names
    marker_specificity = colocalization_result.marker_specificity

    if verbose:
        logging.info(f"Discovering profiles from {len(pairs)} marker pairs...")

    # Step 1: Apply BH-FDR correction to p-values
    # Default: bivariate_morans (spatial cross-correlation), works for mixed and high-seg data
    fdr_pairs, _, method_used = _apply_fdr_correction(
        pairs, alpha=fdr_alpha, pvalue_source=pvalue_source
    )

    if verbose:
        if method_used == "fdr":
            logging.info(
                f"FDR correction [{pvalue_source}] (q < {fdr_alpha}): {len(fdr_pairs)}/{len(pairs)} pairs"
            )
        else:
            logging.info(
                f"Raw p-value filter [{pvalue_source}] (fallback): {len(fdr_pairs)}/{len(pairs)} pairs"
            )

    # Step 2: Apply mutual top-k sparsification WITH specificity weighting
    # This downweights edges involving promiscuous markers (low Gini, like CD45)
    if top_k > 0:
        # Use specificity if available, otherwise None (no weighting)
        specificity_for_topk = marker_specificity if marker_specificity else None
        topk_pairs = _apply_mutual_topk(fdr_pairs, k=top_k, specificity=specificity_for_topk)
        if verbose:
            logging.info(
                f"Mutual top-{top_k} sparsification (specificity-weighted): {len(topk_pairs)}/{len(fdr_pairs)} pairs"
            )
    else:
        topk_pairs = fdr_pairs

    # Step 3: Apply optional score threshold (GMM-adaptive or fixed)
    if min_score is None and len(topk_pairs) > 10:
        # Only use GMM if we have enough pairs and user didn't disable
        scores = np.array([p.colocalization_score for p in topk_pairs])

        n_components, score_threshold, is_signal = _fit_score_gmm_with_bic(
            scores, max_components=5, seed=seed
        )

        if verbose:
            logging.info(
                f"GMM-BIC selected {n_components} components, "
                f"score threshold = {score_threshold:.3f}"
            )

        # Keep only signal pairs
        significant_pairs = [
            p for p, sig in zip(topk_pairs, is_signal) if sig
        ]
    elif min_score is not None:
        # Fixed threshold
        score_threshold = min_score
        significant_pairs = [
            p for p in topk_pairs
            if p.colocalization_score >= min_score
        ]

        if verbose:
            logging.info(f"Fixed score filter (>= {min_score}): {len(significant_pairs)} pairs")
    else:
        # No additional filtering - FDR + top-k is sufficient
        significant_pairs = topk_pairs
        score_threshold = 0.0

    significant_edges = [
        (p.marker_a, p.marker_b, p.colocalization_score)
        for p in significant_pairs
    ]

    if verbose:
        logging.info(
            f"Final: {len(significant_edges)} edges for graph construction"
        )

    # Branch: Triangle assembly vs hierarchical clustering
    if use_triangle_assembly:
        # Triangle-first assembly approach
        if verbose:
            logging.info("Using triangle-first assembly...")

        profiles = _triangle_assembly(all_markers, significant_edges, verbose=verbose)

        # Identify singletons for result tracking
        singletons = [p[0] for p in profiles if len(p) == 1]
        lineage_dendrograms = {}  # No dendrograms in triangle mode

        # Compute modularity
        modularity = _compute_modularity(profiles, pairs, significant_pairs)

        if verbose:
            logging.info(f"Discovered {len(profiles)} profiles, modularity = {modularity:.3f}")

        result = ProfileDiscoveryResult(
            profiles=profiles,
            lineage_dendrograms=lineage_dendrograms,
            singletons=singletons,
            modularity=modularity,
            n_significant_edges=len(significant_edges),
            alpha=fdr_alpha,
        )

        if verbose:
            logging.info(result.summary())

        return result

    # Step 2: Find connected components (hierarchical approach)
    components = _find_connected_components(all_markers, significant_edges)

    # Separate singletons (components of size 1 with no significant edges)
    singletons = []
    multi_marker_components = []

    for comp in components:
        if len(comp) == 1:
            marker = list(comp)[0]
            # Check if this marker has ANY significant edge
            has_edge = any(
                (p.marker_a == marker or p.marker_b == marker)
                for p in significant_pairs
            )
            if not has_edge:
                singletons.append(marker)
            else:
                multi_marker_components.append(comp)
        else:
            multi_marker_components.append(comp)

    if verbose:
        logging.info(
            f"Found {len(multi_marker_components)} connected components, "
            f"{len(singletons)} singletons"
        )

    # Step 3: Hierarchical clustering and gap-based lineage splitting
    profiles = []
    lineage_dendrograms = {}
    lineage_idx = 0

    for comp_idx, component in enumerate(multi_marker_components):
        comp_markers = sorted(list(component))
        n_comp = len(comp_markers)

        if verbose:
            logging.info(f"Processing component {comp_idx + 1} ({n_comp} markers): {comp_markers}")

        if n_comp == 1:
            profiles.append(comp_markers)
            continue

        if n_comp == 2:
            # 2-node component already passed FDR + mutual top-k evidence criteria
            # No additional filtering needed - return as pair profile
            profiles.append(comp_markers)
            continue

        # Build distance matrix and dendrogram for the full component
        dist_matrix = _build_distance_matrix(comp_markers, pairs)
        condensed_dist = squareform(dist_matrix)
        Z = linkage(condensed_dist, method='average')

        # Step 3a: Split into separate lineages based on large gaps
        lineages = _split_dendrogram_by_gaps(Z, comp_markers, seed=seed)

        if verbose and len(lineages) > 1:
            logging.info(f"  Gap analysis split into {len(lineages)} lineages")

        # Step 3b: Process each lineage separately
        for lineage_markers in lineages:
            n_lineage = len(lineage_markers)

            if n_lineage == 1:
                profiles.append(lineage_markers)
                continue

            if n_lineage == 2:
                # 2-node lineage already passed FDR + mutual top-k evidence criteria
                # No additional filtering needed - return as pair profile
                profiles.append(lineage_markers)
                continue

            # Build lineage-specific dendrogram
            lin_dist_matrix = _build_distance_matrix(lineage_markers, pairs)
            lin_condensed = squareform(lin_dist_matrix)
            lin_Z = linkage(lin_condensed, method='average')

            # Store dendrogram for this lineage
            lineage_dendrograms[lineage_idx] = LineageDendrogram(
                markers=lineage_markers,
                linkage_matrix=lin_Z,
                distance_matrix=lin_dist_matrix,
            )
            lineage_idx += 1

            if verbose:
                logging.info(f"  Lineage {lineage_idx}: {lineage_markers}")

            # Dynamic tree cutting within lineage
            cluster_labels = _dynamic_tree_cut(lin_Z, n_lineage)

            # Extract profiles from clusters
            cluster_to_markers = defaultdict(list)
            for i, label in enumerate(cluster_labels):
                cluster_to_markers[label].append(lineage_markers[i])

            for cluster_markers in cluster_to_markers.values():
                profiles.append(cluster_markers)

    # Add singletons as individual profiles
    for singleton in singletons:
        profiles.append([singleton])

    # Step 5: Compute modularity
    modularity = _compute_modularity(profiles, pairs, significant_pairs)

    if verbose:
        logging.info(f"Discovered {len(profiles)} profiles, modularity = {modularity:.3f}")

    result = ProfileDiscoveryResult(
        profiles=profiles,
        lineage_dendrograms=lineage_dendrograms,
        singletons=singletons,
        modularity=modularity,
        n_significant_edges=len(significant_edges),
        alpha=fdr_alpha,
    )

    if verbose:
        logging.info(result.summary())

    return result


def _compute_profile_confidence(
    profile: List[str],
    pairs: List[MarkerPairColocalization],
) -> float:
    """
    Compute profile confidence as mean colocalization score of internal edges.

    Args:
        profile: List of marker names in this profile
        pairs: All colocalization pairs

    Returns:
        Confidence score in [0, 1] (0.5 for singletons)
    """
    if len(profile) <= 1:
        return 0.5  # Default for singletons

    profile_set = set(profile)
    internal_scores = []

    for p in pairs:
        if p.marker_a in profile_set and p.marker_b in profile_set:
            internal_scores.append(p.colocalization_score)

    if not internal_scores:
        return 0.5

    return float(np.mean(internal_scores))


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

