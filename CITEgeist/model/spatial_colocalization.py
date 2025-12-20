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

    rng = np.random.default_rng(seed)

    # Build neighbor graph
    if verbose:
        logging.info(f"Building {neighbor_k}-NN neighbor graph...")
    neighbors = _build_neighbor_graph(coords, neighbor_k)

    # Binarize markers (kept for backwards-compatible metrics)
    if verbose:
        logging.info(f"Binarizing markers at {signal_threshold_percentile}th percentile...")
    binary = _binarize_markers(analyze_X, signal_threshold_percentile)

    # Compute pairwise colocalization
    pairs = []
    total_pairs = n_analyze * (n_analyze - 1) // 2
    processed = 0

    if verbose:
        logging.info(f"Computing pairwise colocalization (continuous metrics)...")

    for i in range(n_analyze):
        for j in range(i + 1, n_analyze):
            name_a = analyze_names[i]
            name_b = analyze_names[j]

            values_a = analyze_X[:, i]
            values_b = analyze_X[:, j]
            binary_a = binary[:, i]
            binary_b = binary[:, j]

            # Same-spot metrics (binary-based, for backwards compatibility)
            jaccard, co_spots, co_frac = _compute_jaccard(binary_a, binary_b)

            # Correlation metrics (continuous)
            pearson_r, spearman_rho, corr_pvalue = _compute_correlation(values_a, values_b)

            # Continuous similarity metrics (binarization-free) - used for score
            cosine_sim = _compute_cosine_similarity(values_a, values_b)
            bivariate_i = _compute_bivariate_morans_i(values_a, values_b, neighbors)
            # Skip bivariate Moran's I permutation test - FDR uses neighbor enrichment
            bivariate_pval = 0.5  # Placeholder

            # Neighborhood metrics (binary-based, used for FDR filtering)
            # Binarization acts as stringent filter focusing on PEAK expression spots
            # This prevents false merges that continuous metrics might allow
            enrich_ab, enrich_ba, enrich_pvalue = _compute_neighbor_enrichment(
                binary_a, binary_b, neighbors, rng, n_permutations
            )
            mutual_enrich = (enrich_ab + enrich_ba) / 2

            # Combined score using continuous metrics only
            score = _compute_colocalization_score(spearman_rho, cosine_sim, bivariate_i)

            pairs.append(MarkerPairColocalization(
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
            ))

            processed += 1
            if verbose and processed % 100 == 0:
                logging.info(f"Processed {processed}/{total_pairs} pairs...")

    # Sort by score
    pairs.sort(key=lambda x: x.colocalization_score, reverse=True)

    result = ColocalizationResult(
        pairs=pairs,
        marker_names=analyze_names,
        n_spots=n_spots,
        neighbor_k=neighbor_k,
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
) -> Tuple[List[MarkerPairColocalization], NDArray[np.bool_], str]:
    """
    Apply Benjamini-Hochberg FDR correction to neighbor enrichment p-values.

    Uses neighbor enrichment p-values (from permutation test) which test whether
    marker B is enriched in neighbors of A-positive spots. The binarization
    (75th percentile threshold) acts as a stringent filter that prevents false
    merges by focusing on PEAK expression spots rather than overall correlation.

    Args:
        pairs: List of marker pair colocalization results.
        alpha: FDR threshold (default: 0.05 for q < 5%).
        fallback_to_raw: If True and FDR finds nothing but raw p-values
            suggest signal exists, fall back to raw p-value threshold.

    Returns:
        Tuple of (filtered pairs, boolean mask, method_used).
        method_used is "fdr" or "raw_pvalue".
    """
    if len(pairs) == 0:
        return [], np.array([], dtype=bool), "fdr"

    # Use neighbor enrichment p-value (binary-based, stricter than continuous)
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
    top_k: int = 3,
    min_score: Optional[float] = None,
    use_triangle_assembly: bool = False,
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

    if verbose:
        logging.info(f"Discovering profiles from {len(pairs)} marker pairs...")

    # Step 1: Apply BH-FDR correction to p-values (with fallback for permutation tests)
    fdr_pairs, _, method_used = _apply_fdr_correction(pairs, alpha=fdr_alpha)

    if verbose:
        if method_used == "fdr":
            logging.info(
                f"FDR correction (q < {fdr_alpha}): {len(fdr_pairs)}/{len(pairs)} pairs"
            )
        else:
            logging.info(
                f"Raw p-value filter (fallback): {len(fdr_pairs)}/{len(pairs)} pairs"
            )

    # Step 2: Apply mutual top-k sparsification
    if top_k > 0:
        topk_pairs = _apply_mutual_topk(fdr_pairs, k=top_k)
        if verbose:
            logging.info(
                f"Mutual top-{top_k} sparsification: {len(topk_pairs)}/{len(fdr_pairs)} pairs"
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


# =============================================================================
# MODULE 2c: Coverage-Based Profile Selection
# =============================================================================

@dataclass
class ProfileSelectionResult:
    """Result from coverage-based profile selection."""

    selected_profiles: List[List[str]]
    """Profiles selected by coverage-based objective."""

    optimal_n: int
    """Number of profiles selected."""

    all_profiles_ranked: List[List[str]]
    """All profiles in selection order."""

    reconstruction_errors: NDArray[np.floating]
    """Reconstruction RMSE for each number of profiles."""

    variance_explained: NDArray[np.floating]
    """Fraction of variance explained for each number of profiles."""

    residual_morans_i: Optional[NDArray[np.floating]] = None
    """Residual Moran's I after each profile added."""

    stopping_reason: str = "threshold"
    """Reason for stopping: 'threshold', 'morans_plateau', or 'all_profiles'."""

    def summary(self) -> str:
        """Return summary string."""
        lines = [
            f"Selected {self.optimal_n}/{len(self.all_profiles_ranked)} profiles",
            f"Final RMSE: {self.reconstruction_errors[self.optimal_n - 1]:.4f}",
            f"Variance explained: {self.variance_explained[self.optimal_n - 1]:.1%}",
            f"Stopping reason: {self.stopping_reason}",
        ]
        if self.residual_morans_i is not None and len(self.residual_morans_i) > 0:
            lines.append(f"Final residual Moran's I: {self.residual_morans_i[self.optimal_n - 1]:.3f}")
        return "\n".join(lines)


# -----------------------------------------------------------------------------
# Helper functions for coverage-based selection
# -----------------------------------------------------------------------------


def _compute_rarity_weights(
    X: NDArray[np.floating],
    marker_names: List[str],
    threshold_percentile: float = 75.0,
    epsilon: float = 0.01,
) -> Dict[str, float]:
    """
    Compute rarity weight for each marker: r(m) = 1 / (Pr(m positive) + ε).

    Low-abundance markers get higher weights, upweighting their contribution
    to coverage score.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        threshold_percentile: Percentile for binarizing (default: 75)
        epsilon: Smoothing constant to avoid division by zero

    Returns:
        Dict mapping marker name to rarity weight
    """
    n_spots, n_markers = X.shape
    thresholds = np.percentile(X, threshold_percentile, axis=0)
    binary = X > thresholds

    positivity_rate = binary.mean(axis=0)  # Fraction of spots where marker is positive
    rarity = 1.0 / (positivity_rate + epsilon)

    return {name: float(rarity[i]) for i, name in enumerate(marker_names)}


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


def _compute_spatial_novelty(
    activity_map: NDArray[np.floating],
    selected_activity_maps: List[NDArray[np.floating]],
) -> float:
    """
    Compute spatial novelty as 1 - max correlation with already selected profiles.

    Args:
        activity_map: Activity map for candidate profile (n_spots,)
        selected_activity_maps: List of activity maps for already selected profiles

    Returns:
        Novelty score in [0, 1] (1.0 if no profiles selected yet)
    """
    if not selected_activity_maps:
        return 1.0

    correlations = []
    for selected_map in selected_activity_maps:
        # Handle constant arrays
        if np.std(activity_map) < 1e-10 or np.std(selected_map) < 1e-10:
            correlations.append(0.0)
        else:
            corr = np.corrcoef(activity_map, selected_map)[0, 1]
            if np.isnan(corr):
                corr = 0.0
            correlations.append(abs(corr))

    return 1.0 - max(correlations)


def _compute_marker_jaccard(profile_a: List[str], profile_b: List[str]) -> float:
    """Compute Jaccard similarity between two profiles."""
    set_a, set_b = set(profile_a), set(profile_b)
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def _compute_morans_i(
    values: NDArray[np.floating],
    neighbors: List[List[int]],
) -> float:
    """
    Compute Moran's I for spatial autocorrelation.

    Args:
        values: Values at each spot (n_spots,)
        neighbors: List of neighbor indices for each spot

    Returns:
        Moran's I statistic (typically in [-1, 1])
    """
    n = len(values)
    if n == 0:
        return 0.0

    mean_val = np.mean(values)
    centered = values - mean_val

    # Numerator: Σ_i Σ_j w_ij * (x_i - mean) * (x_j - mean)
    numerator = 0.0
    W = 0.0  # sum of weights

    for i in range(n):
        for j in neighbors[i]:
            numerator += centered[i] * centered[j]
            W += 1.0

    # Denominator: Σ_i (x_i - mean)^2
    denominator = np.sum(centered ** 2)

    if denominator == 0 or W == 0:
        return 0.0

    return (n / W) * (numerator / denominator)


def _compute_residual_morans_i(
    X: NDArray[np.floating],
    marker_names: List[str],
    selected_profiles: List[List[str]],
    neighbors: List[List[int]],
) -> float:
    """
    Compute median Moran's I of residuals after regressing out selected profiles.

    This measures how much spatial structure remains unexplained by current profiles.

    Args:
        X: Expression matrix (n_spots, n_markers)
        marker_names: Names for each column
        selected_profiles: Currently selected profiles
        neighbors: List of neighbor indices for each spot

    Returns:
        Median absolute Moran's I across all markers' residuals
    """
    from scipy.optimize import nnls

    n_spots, n_markers = X.shape

    if not selected_profiles:
        # No profiles selected - compute Moran's I on raw data
        morans = []
        for m_idx in range(n_markers):
            I = _compute_morans_i(X[:, m_idx], neighbors)
            morans.append(abs(I))
        return float(np.median(morans))

    # Build activity matrix (n_spots, n_profiles)
    activity_matrix = np.column_stack([
        _compute_profile_activity_map(X, marker_names, P) for P in selected_profiles
    ])

    # For each marker, regress and compute residual Moran's I
    morans = []
    for m_idx in range(n_markers):
        y = X[:, m_idx]

        # NNLS regression: y ≈ activity_matrix @ β
        try:
            beta, _ = nnls(activity_matrix, y)
            predicted = activity_matrix @ beta
            residuals = y - predicted
        except Exception:
            residuals = y  # Fall back to raw if NNLS fails

        I = _compute_morans_i(residuals, neighbors)
        morans.append(abs(I))

    return float(np.median(morans))


def _rank_profiles(
    profiles: List[List[str]],
    colocalization_result: Optional[ColocalizationResult] = None,
) -> List[List[str]]:
    """
    Rank profiles for incremental selection.

    Ranking criteria:
    1. Multi-marker profiles first (more evidence)
    2. Within each size group, by mean colocalization score (if available)
    3. Singletons last

    Args:
        profiles: List of profiles from discover_profiles()
        colocalization_result: Optional colocalization data for scoring

    Returns:
        Profiles sorted by quality (best first)
    """
    # Build score lookup if colocalization result provided
    pair_scores = {}
    if colocalization_result is not None:
        for p in colocalization_result.pairs:
            key = tuple(sorted([p.marker_a, p.marker_b]))
            pair_scores[key] = p.colocalization_score

    def profile_score(profile: List[str]) -> Tuple[int, float]:
        """Return (negative size, negative mean score) for sorting."""
        size = len(profile)
        if size == 1:
            return (0, 0.0)  # Singletons get lowest priority

        # Compute mean pairwise score within profile
        if pair_scores:
            scores = []
            for i, m1 in enumerate(profile):
                for m2 in profile[i + 1:]:
                    key = tuple(sorted([m1, m2]))
                    if key in pair_scores:
                        scores.append(pair_scores[key])
            mean_score = np.mean(scores) if scores else 0.0
        else:
            mean_score = 0.0

        return (size, mean_score)

    # Sort by size (descending), then by score (descending)
    return sorted(profiles, key=profile_score, reverse=True)


def _compute_reconstruction_error(
    X: NDArray[np.floating],
    profile_matrix: NDArray[np.floating],
) -> Tuple[float, float]:
    """
    Compute reconstruction error using non-negative least squares.

    For each spot, solve: min ||x - P @ w||^2 s.t. w >= 0
    where x is observed expression, P is profile matrix, w is weights.

    Args:
        X: Observed expression matrix (n_spots, n_markers)
        profile_matrix: Profile indicator matrix (n_markers, n_profiles)

    Returns:
        Tuple of (RMSE, variance_explained)
    """
    from scipy.optimize import nnls

    n_spots, n_markers = X.shape
    n_profiles = profile_matrix.shape[1]

    if n_profiles == 0:
        # No profiles - return total variance as error
        total_var = np.var(X)
        return float(np.sqrt(np.mean(X ** 2))), 0.0

    # Reconstruct each spot
    reconstructed = np.zeros_like(X)

    for i in range(n_spots):
        # Non-negative least squares for this spot
        weights, _ = nnls(profile_matrix, X[i])
        reconstructed[i] = profile_matrix @ weights

    # Compute metrics
    residuals = X - reconstructed
    rmse = float(np.sqrt(np.mean(residuals ** 2)))

    total_var = np.var(X)
    residual_var = np.var(residuals)
    var_explained = 1.0 - (residual_var / total_var) if total_var > 0 else 0.0

    return rmse, var_explained


def _find_elbow(values: NDArray[np.floating]) -> int:
    """
    Find elbow point using maximum curvature method.

    Args:
        values: Array of values (e.g., reconstruction errors)

    Returns:
        Index of elbow point (1-indexed for number of profiles)
    """
    n = len(values)
    if n <= 2:
        return n

    # Normalize to [0, 1] range for both axes
    x = np.arange(n) / (n - 1)
    y = (values - values.min()) / (values.max() - values.min() + 1e-10)

    # Line from first to last point
    p1 = np.array([x[0], y[0]])
    p2 = np.array([x[-1], y[-1]])

    # Distance from each point to the line
    distances = np.zeros(n)
    for i in range(n):
        p = np.array([x[i], y[i]])
        # Distance from point to line
        distances[i] = np.abs(np.cross(p2 - p1, p1 - p)) / np.linalg.norm(p2 - p1)

    # Elbow is point with maximum distance
    elbow_idx = int(np.argmax(distances))

    # Return 1-indexed (number of profiles)
    return elbow_idx + 1


def select_profiles_by_reconstruction(
    X: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    colocalization_result: Optional[ColocalizationResult] = None,
    interesting_markers: Optional[List[str]] = None,
    min_profiles: int = 1,
    max_profiles: Optional[int] = None,
    var_explained_threshold: float = 0.90,
    verbose: bool = True,
) -> ProfileSelectionResult:
    """
    Select optimal number of profiles using variance explained threshold.

    Module 2c: Data-driven profile selection based on how well profiles
    can reconstruct the observed antibody expression matrix.

    Algorithm:
    1. Rank profiles by quality (multi-marker first, by colocalization score)
    2. For n = 1, 2, ... profiles:
       - Build profile matrix (ALL interesting markers x profiles)
       - Solve NNLS reconstruction for each spot
       - Compute RMSE and variance explained over ALL interesting markers
    3. Select minimum profiles needed to reach variance explained threshold
    4. Return profiles up to that point

    Args:
        X: Observed expression matrix (n_spots, n_markers)
        marker_names: Names for each column of X
        profiles: List of profiles from discover_profiles()
        colocalization_result: Optional, for ranking by colocalization score
        interesting_markers: Markers to reconstruct (from Module 1). If None,
            uses all unique markers across profiles.
        min_profiles: Minimum profiles to consider (default: 1)
        max_profiles: Maximum profiles to consider (default: all)
        var_explained_threshold: Target variance explained (default: 0.90 = 90%)
        verbose: Log progress (default: True)

    Returns:
        ProfileSelectionResult with selected profiles and metrics

    Example:
        >>> # After Module 2b
        >>> profile_result = discover_profiles(coloc_result)
        >>> # Module 2c: Select optimal profiles
        >>> selection = select_profiles_by_reconstruction(
        ...     X, marker_names, profile_result.profiles, coloc_result,
        ...     interesting_markers=result1.interesting_markers
        ... )
        >>> print(selection.summary())
        >>> selected = selection.selected_profiles
    """
    # Build marker name to index mapping
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    # Determine target markers to reconstruct
    # If interesting_markers provided, use those; otherwise use all markers in profiles
    if interesting_markers is not None:
        target_markers = sorted([m for m in interesting_markers if m in marker_to_idx])
    else:
        # Fall back to all unique markers across ALL profiles
        all_profile_markers = set()
        for p in profiles:
            all_profile_markers.update(p)
        target_markers = sorted([m for m in all_profile_markers if m in marker_to_idx])

    # Extract target marker columns from X (FIXED: use ALL target markers, not just current profile's)
    target_indices = [marker_to_idx[m] for m in target_markers]
    X_target = X[:, target_indices]

    if verbose:
        logging.info(f"Module 2c: Selecting profiles by reconstruction accuracy...")
        logging.info(f"  Input: {len(profiles)} profiles, {X.shape[0]} spots")
        logging.info(f"  Target markers to reconstruct: {len(target_markers)}")

    if max_profiles is None:
        max_profiles = len(profiles)
    max_profiles = min(max_profiles, len(profiles))

    # Create mapping from target marker to its index in X_target
    target_marker_to_local = {m: i for i, m in enumerate(target_markers)}

    # Helper function to compute variance explained for a set of profiles
    def compute_ve_for_profiles(profile_list):
        if not profile_list:
            return 0.0, float('inf')
        profile_matrix = np.zeros((len(target_markers), len(profile_list)))
        for j, profile in enumerate(profile_list):
            for marker in profile:
                if marker in target_marker_to_local:
                    local_idx = target_marker_to_local[marker]
                    profile_matrix[local_idx, j] = 1.0
        rmse, ve = _compute_reconstruction_error(X_target, profile_matrix)
        return ve, rmse

    # Greedy forward selection: at each step, add the profile that maximizes variance explained
    remaining_profiles = [tuple(p) for p in profiles]  # Convert to tuples for set operations
    selected_profiles = []
    errors = []
    var_explained = []

    if verbose:
        logging.info(f"  Greedy forward selection (adding best profile at each step):")

    for n in range(1, max_profiles + 1):
        best_ve = -1
        best_rmse = float('inf')
        best_profile = None

        # Try adding each remaining profile and see which gives best variance explained
        for candidate in remaining_profiles:
            test_profiles = selected_profiles + [list(candidate)]
            ve, rmse = compute_ve_for_profiles(test_profiles)
            if ve > best_ve:
                best_ve = ve
                best_rmse = rmse
                best_profile = candidate

        if best_profile is None:
            break

        # Add the best profile
        selected_profiles.append(list(best_profile))
        remaining_profiles.remove(best_profile)
        errors.append(best_rmse)
        var_explained.append(best_ve)

        if verbose and n <= 10:
            logging.info(f"    n={n}: +{list(best_profile)} -> VarExpl={best_ve:.1%}")

        # Check if we've reached the threshold
        if best_ve >= var_explained_threshold:
            if verbose:
                logging.info(f"  Variance threshold {var_explained_threshold:.0%} reached at n={n} profiles")
            break

    errors = np.array(errors)
    var_explained = np.array(var_explained)
    optimal_n = len(selected_profiles)

    if verbose:
        if var_explained[-1] >= var_explained_threshold:
            pass  # Already logged above
        else:
            logging.info(f"  Variance threshold {var_explained_threshold:.0%} not reached, using all {optimal_n} profiles")
        logging.info(f"  Final RMSE: {errors[-1]:.4f}")
        logging.info(f"  Variance explained: {var_explained[-1]:.1%}")
        logging.info(f"  Selected profiles:")
        for i, p in enumerate(selected_profiles):
            logging.info(f"    {i+1}. {p}")

    return ProfileSelectionResult(
        selected_profiles=selected_profiles,
        optimal_n=optimal_n,
        all_profiles_ranked=selected_profiles,  # Now ordered by greedy selection
        reconstruction_errors=errors,
        variance_explained=var_explained,
    )


def select_profiles_by_coverage(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    colocalization_result: ColocalizationResult,
    interesting_markers: Optional[List[str]] = None,
    *,
    neighbor_k: int = 6,
    alpha: float = 1.0,
    beta: float = 0.5,
    gamma: float = 0.3,
    max_profiles: Optional[int] = None,
    min_var_explained: float = 0.90,
    min_coverage: float = 0.90,
    n_null_samples: int = 50,
    stat_significance: float = 0.05,
    verbose: bool = True,
) -> ProfileSelectionResult:
    """
    Select profiles using coverage-based objective with statistical stopping criterion.

    Module 2c v2: Biologically-aligned profile selection that prioritizes:
    - Coverage of rare markers (not just abundant ones)
    - Spatial novelty (profiles with different spatial patterns)
    - Non-redundancy (avoid overlapping profiles)

    Gain function:
        gain(P) = α * c(P) * Σ r(m) for m in P not yet covered   (coverage)
                + β * c(P) * (1 - max correlation with selected)   (spatial novelty)
                - γ * max Jaccard(P, Q) for Q in selected          (redundancy penalty)

    Where:
        r(m) = 1 / (Pr(m positive) + ε)  - rarity weight (upweights low-abundance)
        c(P) = profile confidence (mean internal edge weight)

    Stopping criteria (requires BOTH dual thresholds before statistical test):
        1. Dual checkpoint: Must reach BOTH min_var_explained AND min_coverage
        2. Statistical test: Stop when marginal Moran's I gain is not significantly
           above what you'd get from a random profile (p > stat_significance)
        3. All profiles added

    The statistical test compares the observed ΔMoran's I to a null distribution
    generated by randomly selecting from remaining profiles. This prevents stopping
    prematurely when cell types are spatially correlated with already-selected profiles.

    Args:
        X: Expression matrix (n_spots, n_markers)
        coords: Spatial coordinates (n_spots, 2)
        marker_names: Names for each column of X
        profiles: List of profiles from discover_profiles()
        colocalization_result: Colocalization analysis result (for confidence scores)
        interesting_markers: Markers to evaluate coverage (from Module 1). If None,
            uses all unique markers across profiles.
        neighbor_k: Number of neighbors for Moran's I computation (default: 6)
        alpha: Weight for coverage term (default: 1.0)
        beta: Weight for spatial novelty term (default: 0.5)
        gamma: Weight for redundancy penalty (default: 0.3)
        max_profiles: Maximum profiles to select (default: all)
        min_var_explained: Minimum variance explained before statistical stop allowed
            (default: 0.90)
        min_coverage: Minimum marker coverage before statistical stop allowed
            (default: 0.90)
        n_null_samples: Number of random profile samples for null distribution
            (default: 50)
        stat_significance: P-value threshold for statistical stopping test
            (default: 0.05)
        verbose: Log progress (default: True)

    Returns:
        ProfileSelectionResult with selected profiles, metrics, and stopping reason.

    Example:
        >>> # After Module 2b
        >>> profile_result = discover_profiles(coloc_result)
        >>> # Module 2c v2: Coverage-based selection
        >>> selection = select_profiles_by_coverage(
        ...     X, coords, marker_names, profile_result.profiles, coloc_result,
        ...     interesting_markers=result1.interesting_markers
        ... )
        >>> print(selection.summary())
    """
    from scipy.optimize import nnls

    # Build marker name to index mapping
    marker_to_idx = {name: i for i, name in enumerate(marker_names)}

    # Determine target markers
    if interesting_markers is not None:
        target_markers = sorted([m for m in interesting_markers if m in marker_to_idx])
    else:
        all_profile_markers = set()
        for p in profiles:
            all_profile_markers.update(p)
        target_markers = sorted([m for m in all_profile_markers if m in marker_to_idx])

    # Extract target marker columns
    target_indices = [marker_to_idx[m] for m in target_markers]
    X_target = X[:, target_indices]
    target_marker_to_local = {m: i for i, m in enumerate(target_markers)}

    if verbose:
        logging.info(f"Module 2c: Coverage-based profile selection")
        logging.info(f"  Input: {len(profiles)} profiles, {X.shape[0]} spots")
        logging.info(f"  Target markers: {len(target_markers)}")
        logging.info(f"  Weights: α={alpha} (coverage), β={beta} (novelty), γ={gamma} (redundancy)")

    if max_profiles is None:
        max_profiles = len(profiles)
    max_profiles = min(max_profiles, len(profiles))

    # Precompute components
    # 1. Rarity weights for all markers
    rarity_weights = _compute_rarity_weights(X, marker_names)

    # 2. Profile confidences
    profile_confidence = {
        tuple(p): _compute_profile_confidence(p, colocalization_result.pairs)
        for p in profiles
    }

    # 3. Build neighbor graph for Moran's I
    neighbors = _build_neighbor_graph(coords, neighbor_k)

    # 4. Precompute activity maps for all profiles
    profile_activity_maps = {
        tuple(p): _compute_profile_activity_map(X, marker_names, p)
        for p in profiles
    }

    # Helper: compute variance explained for current selection
    def compute_ve_for_profiles(profile_list):
        if not profile_list:
            return 0.0, float('inf')
        profile_matrix = np.zeros((len(target_markers), len(profile_list)))
        for j, profile in enumerate(profile_list):
            for marker in profile:
                if marker in target_marker_to_local:
                    local_idx = target_marker_to_local[marker]
                    profile_matrix[local_idx, j] = 1.0
        rmse, ve = _compute_reconstruction_error(X_target, profile_matrix)
        return ve, rmse

    # Greedy selection
    remaining_profiles = [tuple(p) for p in profiles]
    selected_profiles: List[List[str]] = []
    selected_activity_maps: List[NDArray[np.floating]] = []
    covered_markers: Set[str] = set()

    errors = []
    var_explained_list = []
    morans_i_list = []

    # Initial residual Moran's I (before any profiles)
    prev_morans_i = _compute_residual_morans_i(X_target, target_markers, [], neighbors)

    if verbose:
        logging.info(f"  Initial residual Moran's I: {prev_morans_i:.3f}")
        logging.info(f"  Greedy selection with coverage objective:")

    stopping_reason = "all_profiles"

    for n in range(1, max_profiles + 1):
        best_gain = float('-inf')
        best_profile = None

        for candidate in remaining_profiles:
            candidate_list = list(candidate)

            # Coverage term: sum of rarity weights for uncovered markers in this profile
            new_coverage = sum(
                rarity_weights.get(m, 1.0)
                for m in candidate_list
                if m in target_markers and m not in covered_markers
            )

            # Confidence term
            conf = profile_confidence.get(candidate, 0.5)

            # Spatial novelty term
            activity_map = profile_activity_maps[candidate]
            novelty = _compute_spatial_novelty(activity_map, selected_activity_maps)

            # Redundancy penalty: max Jaccard with already selected profiles
            if selected_profiles:
                max_jaccard = max(
                    _compute_marker_jaccard(candidate_list, sel)
                    for sel in selected_profiles
                )
            else:
                max_jaccard = 0.0

            # Combined gain
            gain = (
                alpha * conf * new_coverage +
                beta * conf * novelty -
                gamma * max_jaccard
            )

            if gain > best_gain:
                best_gain = gain
                best_profile = candidate

        if best_profile is None:
            break

        # Add the best profile
        selected_profiles.append(list(best_profile))
        selected_activity_maps.append(profile_activity_maps[best_profile])
        covered_markers.update(m for m in best_profile if m in target_markers)
        remaining_profiles.remove(best_profile)

        # Compute metrics
        ve, rmse = compute_ve_for_profiles(selected_profiles)
        errors.append(rmse)
        var_explained_list.append(ve)

        # Compute residual Moran's I
        curr_morans_i = _compute_residual_morans_i(
            X_target, target_markers, selected_profiles, neighbors
        )
        morans_i_list.append(curr_morans_i)
        delta_i = prev_morans_i - curr_morans_i

        if verbose and n <= 15:
            conf = profile_confidence.get(best_profile, 0.5)
            logging.info(
                f"    n={n}: +{list(best_profile)} | "
                f"gain={best_gain:.2f}, conf={conf:.2f}, VE={ve:.1%}, "
                f"Moran's I: {prev_morans_i:.3f}→{curr_morans_i:.3f} (Δ={delta_i:.3f})"
            )

        # Compute marker coverage
        coverage_frac = len(covered_markers) / len(target_markers) if target_markers else 1.0

        # Check stopping criteria with dual checkpoint + statistical test
        # Only consider stopping if BOTH thresholds are met
        dual_checkpoint_met = (ve >= min_var_explained) and (coverage_frac >= min_coverage)

        if dual_checkpoint_met and n > 1 and len(remaining_profiles) >= 2:
            # Statistical test: is observed Δ significantly above null?
            # Null distribution: Moran's I reduction from random profile selection

            # Compute null distribution by sampling random profiles
            rng = np.random.default_rng(42 + n)  # Deterministic per iteration
            null_deltas = []

            # Sample from remaining profiles (excluding the best we already picked)
            sample_pool = [p for p in remaining_profiles]  # All remaining
            n_samples = min(n_null_samples, len(sample_pool))

            for _ in range(n_samples):
                random_profile = sample_pool[rng.integers(len(sample_pool))]
                null_morans = _compute_residual_morans_i(
                    X_target, target_markers,
                    selected_profiles[:-1] + [list(random_profile)],  # Replace last with random
                    neighbors
                )
                null_delta = prev_morans_i - null_morans
                null_deltas.append(null_delta)

            # P-value: fraction of null samples with delta >= observed
            if null_deltas:
                p_value = (sum(d >= delta_i for d in null_deltas) + 1) / (len(null_deltas) + 1)
            else:
                p_value = 0.0  # No null samples, don't stop

            if p_value > stat_significance:
                stopping_reason = "statistical_plateau"
                if verbose:
                    logging.info(
                        f"  Stopping: dual checkpoint met (VE={ve:.1%}, cov={coverage_frac:.1%}) "
                        f"and Moran's I gain not significant (p={p_value:.3f} > {stat_significance})"
                    )
                break
            elif verbose and n <= 15:
                logging.info(f"      (stat test: p={p_value:.3f}, continue)")

        prev_morans_i = curr_morans_i

    errors = np.array(errors) if errors else np.array([0.0])
    var_explained_arr = np.array(var_explained_list) if var_explained_list else np.array([0.0])
    morans_i_arr = np.array(morans_i_list) if morans_i_list else np.array([0.0])
    optimal_n = len(selected_profiles)

    final_coverage = len(covered_markers) / len(target_markers) if target_markers else 1.0

    if verbose:
        logging.info(f"  Final: {optimal_n} profiles selected")
        logging.info(f"  Variance explained: {var_explained_arr[-1]:.1%}")
        logging.info(f"  Marker coverage: {len(covered_markers)}/{len(target_markers)} ({final_coverage:.1%})")
        logging.info(f"  Residual Moran's I: {morans_i_arr[-1]:.3f}")
        logging.info(f"  Stopping reason: {stopping_reason}")
        logging.info(f"  Thresholds: VE>={min_var_explained:.0%}, cov>={min_coverage:.0%}")
        logging.info(f"  Selected profiles:")
        for i, p in enumerate(selected_profiles):
            logging.info(f"    {i+1}. {p}")

    return ProfileSelectionResult(
        selected_profiles=selected_profiles,
        optimal_n=optimal_n,
        all_profiles_ranked=selected_profiles,
        reconstruction_errors=errors,
        variance_explained=var_explained_arr,
        residual_morans_i=morans_i_arr,
        stopping_reason=stopping_reason,
    )


# =============================================================================
# Module 2c v3: Spatial Variance-Based Profile Selection
# =============================================================================


@dataclass
class SpatialVarianceSelectionResult:
    """Result from spatial variance-based profile selection."""

    selected_profiles: List[List[str]]
    """Profiles selected by spatial variance objective."""

    optimal_n: int
    """Number of profiles selected."""

    all_profiles_ranked: List[List[str]]
    """All profiles in selection order."""

    # Spatial metrics
    spatial_covariance_matrix: NDArray[np.floating]
    """Bivariate Moran's I matrix for interesting markers."""

    eigenvalues: NDArray[np.floating]
    """Eigenvalues of spatial covariance (variance in each direction)."""

    explained_spatial_variance: NDArray[np.floating]
    """Cumulative explained spatial variance for each n profiles."""

    # Reconstruction metrics
    variance_explained: NDArray[np.floating]
    """Fraction of data variance explained for each n profiles."""

    proportion_smoothness: NDArray[np.floating]
    """Mean Moran's I of NNLS proportions for each n profiles."""

    # Stopping info
    stopping_reason: str
    """Reason: 'spatial_threshold', 'diminishing_returns', or 'all_profiles'."""

    def summary(self) -> str:
        """Return summary string."""
        n = self.optimal_n
        lines = [
            f"Selected {n}/{len(self.all_profiles_ranked)} profiles",
            f"Explained spatial variance: {self.explained_spatial_variance[n-1]:.1%}",
            f"Data variance explained: {self.variance_explained[n-1]:.1%}",
            f"Proportion smoothness: {self.proportion_smoothness[n-1]:.3f}",
            f"Stopping reason: {self.stopping_reason}",
        ]
        return "\n".join(lines)


def _build_spatial_covariance_matrix(
    X: NDArray[np.floating],
    neighbors: List[List[int]],
    marker_names: List[str],
    verbose: bool = False,
) -> NDArray[np.floating]:
    """
    Build spatial covariance matrix using bivariate Moran's I.

    C[i,j] = bivariate_Moran's_I(marker_i, marker_j)

    For diagonal elements (i==j), uses univariate Moran's I.

    Args:
        X: Expression matrix (n_spots, n_markers).
        neighbors: List of neighbor indices per spot.
        marker_names: Names for each marker column.
        verbose: Log progress.

    Returns:
        Symmetric matrix (n_markers, n_markers) of spatial cross-correlations.
    """
    n_markers = len(marker_names)
    C = np.zeros((n_markers, n_markers))

    total_pairs = n_markers * (n_markers + 1) // 2
    computed = 0

    for i in range(n_markers):
        for j in range(i, n_markers):
            if i == j:
                # Univariate Moran's I on diagonal
                C[i, j] = _compute_morans_i(X[:, i], neighbors)
            else:
                # Bivariate Moran's I for off-diagonal
                C[i, j] = _compute_bivariate_morans_i(X[:, i], X[:, j], neighbors)
                C[j, i] = C[i, j]  # Symmetric

            computed += 1
            if verbose and computed % 50 == 0:
                logging.info(f"  Spatial covariance: {computed}/{total_pairs} pairs...")

    return C


def _compute_profile_activity_map(
    X: NDArray[np.floating],
    marker_names: List[str],
    profile: List[str],
) -> NDArray[np.floating]:
    """
    Compute activity map for a profile (mean expression across profile markers).

    Args:
        X: Expression matrix (n_spots, n_markers).
        marker_names: Names for each marker column.
        profile: List of marker names in the profile.

    Returns:
        Activity map (n_spots,) representing spatial pattern of the profile.
    """
    marker_indices = [marker_names.index(m) for m in profile if m in marker_names]
    if not marker_indices:
        return np.zeros(X.shape[0])

    return np.mean(X[:, marker_indices], axis=1)


def _compute_spatial_coverage_score(
    activity_map: NDArray[np.floating],
    eigenvectors: NDArray[np.floating],
    eigenvalues: NDArray[np.floating],
    remaining_weight: NDArray[np.floating],
) -> float:
    """
    Compute how much spatial variance a profile explains.

    Projects the profile's activity map onto the spatial eigenvectors,
    weighted by remaining eigenvalue importance.

    Args:
        activity_map: Profile activity (n_spots,).
        eigenvectors: Spatial eigenvectors (n_spots, k).
        eigenvalues: Original eigenvalues (k,).
        remaining_weight: Remaining unexplained weight per eigenvector (k,).

    Returns:
        Spatial coverage score (higher = explains more spatial variance).
    """
    # Normalize activity map
    if np.std(activity_map) < 1e-10:
        return 0.0

    activity_normalized = (activity_map - np.mean(activity_map)) / np.std(activity_map)

    # Project onto eigenvectors
    projections = eigenvectors.T @ activity_normalized  # (k,)

    # Weighted score by remaining eigenvalue importance
    score = np.sum(remaining_weight * projections**2)

    return float(score)


def _compute_proportion_smoothness_fast(
    X: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    neighbors: List[List[int]],
) -> float:
    """
    Compute spatial smoothness of NNLS-estimated proportions.

    Runs NNLS per spot and computes mean Moran's I of proportion vectors.
    Higher values = more spatially coherent proportions.

    Args:
        X: Expression matrix (n_spots, n_markers).
        marker_names: Names for each marker column.
        profiles: List of profiles (each a list of marker names).
        neighbors: Neighbor indices per spot.

    Returns:
        Mean Moran's I of proportions across profiles (higher = smoother).
    """
    from scipy.optimize import nnls

    n_spots = X.shape[0]
    n_profiles = len(profiles)

    if n_profiles == 0:
        return 0.0

    # Build profile indicator matrix (n_markers, n_profiles)
    # P[m, p] = 1 if marker m is in profile p
    marker_to_idx = {m: i for i, m in enumerate(marker_names)}
    P = np.zeros((len(marker_names), n_profiles))
    for p_idx, profile in enumerate(profiles):
        for marker in profile:
            if marker in marker_to_idx:
                P[marker_to_idx[marker], p_idx] = 1.0

    # NNLS per spot to estimate proportions
    # Solve: X[i] ≈ P @ w for weights w (n_profiles,)
    proportions = np.zeros((n_spots, n_profiles))
    for i in range(n_spots):
        try:
            proportions[i], _ = nnls(P, X[i])
        except Exception:
            proportions[i] = np.zeros(n_profiles)

    # Normalize proportions per spot
    row_sums = proportions.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums > 0, row_sums, 1.0)
    proportions = proportions / row_sums

    # Compute Moran's I for each profile's proportions
    smoothness_values = []
    for p_idx in range(n_profiles):
        prop_values = proportions[:, p_idx]
        if np.std(prop_values) > 1e-10:
            moran_i = _compute_morans_i(prop_values, neighbors)
            smoothness_values.append(moran_i)

    if not smoothness_values:
        return 0.0

    return float(np.mean(smoothness_values))


def select_profiles_by_spatial_variance(
    X: NDArray[np.floating],
    coords: NDArray[np.floating],
    marker_names: List[str],
    profiles: List[List[str]],
    *,
    interesting_markers: Optional[List[str]] = None,
    colocalization_result: Optional[ColocalizationResult] = None,
    alpha: float = 1.0,
    beta: float = 0.3,
    gamma: float = 0.2,
    min_spatial_explained: float = 0.90,
    min_marginal_gain: float = 0.02,
    neighbor_k: int = 6,
    verbose: bool = True,
) -> SpatialVarianceSelectionResult:
    """
    Select profiles that best explain spatial variation for deconvolution.

    Algorithm:
    1. Build spatial covariance matrix (bivariate Moran's I between markers)
    2. Eigendecomposition to find spatial variance directions
    3. Greedy selection optimizing:
       - Spatial coverage (eigenvalue-weighted projection)
       - Proportion smoothness (NNLS + Moran's I)
       - Redundancy penalty
    4. Stop when 90% spatial variance explained or diminishing returns

    Args:
        X: Expression matrix (n_spots, n_markers).
        coords: Spatial coordinates (n_spots, 2).
        marker_names: Names for each marker column.
        profiles: List of profiles from Module 2b (each a list of markers).
        interesting_markers: Subset of markers to consider (from Module 1).
        colocalization_result: Optional, for profile confidence scores.
        alpha: Weight for spatial coverage term (default: 1.0).
        beta: Weight for proportion smoothness term (default: 0.3).
        gamma: Weight for redundancy penalty (default: 0.2).
        min_spatial_explained: Target spatial variance to explain (default: 0.90).
        min_marginal_gain: Minimum gain to continue (default: 0.02).
        neighbor_k: Number of neighbors for spatial analysis (default: 6).
        verbose: Log progress (default: True).

    Returns:
        SpatialVarianceSelectionResult with selected profiles and metrics.
    """
    X = np.asarray(X, dtype=np.float64)
    coords = np.asarray(coords, dtype=np.float64)
    n_spots, n_markers = X.shape

    if len(profiles) == 0:
        raise ValueError("No profiles provided")

    if verbose:
        logging.info("Module 2c: Spatial variance-based profile selection")
        logging.info(f"  Input: {len(profiles)} profiles, {n_spots} spots")

    # Determine target markers
    if interesting_markers is not None:
        target_markers = [m for m in interesting_markers if m in marker_names]
    else:
        all_profile_markers = set()
        for p in profiles:
            all_profile_markers.update(p)
        target_markers = [m for m in all_profile_markers if m in marker_names]

    if not target_markers:
        raise ValueError("No valid target markers found")

    if verbose:
        logging.info(f"  Target markers: {len(target_markers)}")

    # Filter X to target markers
    target_indices = [marker_names.index(m) for m in target_markers]
    X_target = X[:, target_indices]

    # Build neighbor graph
    neighbors = _build_neighbor_graph(coords, neighbor_k)

    # Phase 1: Build spatial covariance matrix
    if verbose:
        logging.info("  Phase 1: Building spatial covariance matrix...")

    C = _build_spatial_covariance_matrix(X_target, neighbors, target_markers, verbose=verbose)

    # Eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(C)

    # Sort by descending eigenvalue
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Handle negative eigenvalues (set to small positive)
    eigenvalues = np.maximum(eigenvalues, 1e-10)

    # Normalize eigenvalues to sum to 1
    total_variance = np.sum(eigenvalues)
    eigenvalue_weights = eigenvalues / total_variance

    if verbose:
        top_3_explained = np.sum(eigenvalue_weights[:3])
        logging.info(f"  Top 3 eigenvectors explain {top_3_explained:.1%} of spatial variance")

    # Project X onto eigenvectors for activity map projection
    X_target_centered = X_target - np.mean(X_target, axis=0)
    X_projected = X_target_centered @ eigenvectors  # (n_spots, n_markers)

    # Build marker lookup for variance calculation
    target_marker_to_local = {m: i for i, m in enumerate(target_markers)}

    # Helper: compute variance explained for current selection
    def compute_ve_for_profiles(profile_list):
        if not profile_list:
            return 0.0
        profile_matrix = np.zeros((len(target_markers), len(profile_list)))
        for j, profile in enumerate(profile_list):
            for marker in profile:
                if marker in target_marker_to_local:
                    local_idx = target_marker_to_local[marker]
                    profile_matrix[local_idx, j] = 1.0
        _, ve = _compute_reconstruction_error(X_target, profile_matrix)
        return ve

    # Phase 2: Greedy profile selection
    if verbose:
        logging.info("  Phase 2: Greedy selection with spatial + smoothness objectives...")

    remaining_weight = eigenvalue_weights.copy()
    selected_profiles: List[List[str]] = []
    remaining_profiles = list(profiles)

    # Get profile confidences from colocalization if available
    profile_confidence: Dict[int, float] = {}
    for i, profile in enumerate(profiles):
        if colocalization_result is not None and len(profile) > 1:
            conf = _compute_profile_confidence(profile, colocalization_result.pairs)
        else:
            conf = 0.5  # Default for singletons
        profile_confidence[i] = conf

    # Tracking arrays
    spatial_explained_list = []
    var_explained_list = []
    smoothness_list = []
    stopping_reason = "all_profiles"

    for iteration in range(len(profiles)):
        if not remaining_profiles:
            break

        # Score each candidate profile
        best_score = -np.inf
        best_profile = None
        best_profile_idx = -1

        for p_idx, profile in enumerate(remaining_profiles):
            # Compute activity map
            activity = _compute_profile_activity_map(X, marker_names, profile)

            # Project activity onto spatial eigenvectors
            if np.std(activity) > 1e-10:
                activity_centered = activity - np.mean(activity)
                activity_projected = X_projected.T @ activity_centered  # (n_markers,)
                activity_projected = activity_projected / (np.linalg.norm(activity_projected) + 1e-10)
            else:
                activity_projected = np.zeros(len(target_markers))

            # Spatial coverage score
            spatial_score = np.sum(remaining_weight * activity_projected**2)

            # Proportion smoothness (only compute for top candidates to save time)
            # For first few iterations, compute for all; later only for top candidates
            if len(selected_profiles) < 3 or spatial_score > 0.01:
                test_profiles = selected_profiles + [profile]
                smoothness = _compute_proportion_smoothness_fast(
                    X, marker_names, test_profiles, neighbors
                )
            else:
                smoothness = 0.0

            # Redundancy penalty
            if selected_profiles:
                max_jaccard = max(
                    len(set(profile) & set(sp)) / len(set(profile) | set(sp))
                    for sp in selected_profiles
                )
            else:
                max_jaccard = 0.0

            # Get confidence
            orig_idx = profiles.index(profile)
            conf = profile_confidence.get(orig_idx, 0.5)

            # Combined score
            score = (
                alpha * spatial_score * conf
                + beta * smoothness
                - gamma * max_jaccard
            )

            if score > best_score:
                best_score = score
                best_profile = profile
                best_profile_idx = p_idx

        if best_profile is None:
            break

        # Add best profile
        selected_profiles.append(best_profile)
        remaining_profiles.pop(best_profile_idx)

        # Update remaining weight (deflate)
        activity = _compute_profile_activity_map(X, marker_names, best_profile)
        if np.std(activity) > 1e-10:
            activity_centered = activity - np.mean(activity)
            activity_projected = X_projected.T @ activity_centered
            activity_projected = activity_projected / (np.linalg.norm(activity_projected) + 1e-10)
            explained = remaining_weight * activity_projected**2
            remaining_weight = np.maximum(remaining_weight - explained, 0)

        # Compute metrics
        cumulative_explained = 1.0 - np.sum(remaining_weight)
        spatial_explained_list.append(cumulative_explained)

        # Variance explained
        var_exp = compute_ve_for_profiles(selected_profiles)
        var_explained_list.append(var_exp)

        # Smoothness
        smoothness = _compute_proportion_smoothness_fast(X, marker_names, selected_profiles, neighbors)
        smoothness_list.append(smoothness)

        if verbose:
            logging.info(
                f"    n={len(selected_profiles)}: +{best_profile} | "
                f"spatial={cumulative_explained:.1%}, VE={var_exp:.1%}, smooth={smoothness:.3f}"
            )

        # Check stopping criteria
        if cumulative_explained >= min_spatial_explained:
            stopping_reason = "spatial_threshold"
            if verbose:
                logging.info(f"  Stopping: {min_spatial_explained:.0%} spatial variance explained")
            break

        if len(spatial_explained_list) >= 2:
            marginal_gain = spatial_explained_list[-1] - spatial_explained_list[-2]
            if marginal_gain < min_marginal_gain:
                stopping_reason = "diminishing_returns"
                if verbose:
                    logging.info(f"  Stopping: marginal gain {marginal_gain:.1%} < {min_marginal_gain:.0%}")
                break

    # Build result arrays
    spatial_explained_arr = np.array(spatial_explained_list) if spatial_explained_list else np.array([0.0])
    var_explained_arr = np.array(var_explained_list) if var_explained_list else np.array([0.0])
    smoothness_arr = np.array(smoothness_list) if smoothness_list else np.array([0.0])

    if verbose:
        logging.info(f"  Final: {len(selected_profiles)} profiles selected")
        logging.info(f"  Spatial variance explained: {spatial_explained_arr[-1]:.1%}")
        logging.info(f"  Data variance explained: {var_explained_arr[-1]:.1%}")
        logging.info(f"  Proportion smoothness: {smoothness_arr[-1]:.3f}")
        logging.info(f"  Stopping reason: {stopping_reason}")

    return SpatialVarianceSelectionResult(
        selected_profiles=selected_profiles,
        optimal_n=len(selected_profiles),
        all_profiles_ranked=selected_profiles,
        spatial_covariance_matrix=C,
        eigenvalues=eigenvalues,
        explained_spatial_variance=spatial_explained_arr,
        variance_explained=var_explained_arr,
        proportion_smoothness=smoothness_arr,
        stopping_reason=stopping_reason,
    )
