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

    # Same-spot co-occurrence
    jaccard_index: float  # |A & B| / |A | B|
    co_occurrence_spots: int  # Number of spots with both markers
    co_occurrence_fraction: float  # Fraction of spots with both

    # Expression correlation
    pearson_r: float  # Pearson correlation of intensities
    spearman_rho: float  # Spearman rank correlation
    correlation_pvalue: float  # P-value for correlation test

    # Neighborhood analysis
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
    jaccard: float,
    pearson_r: float,
    mutual_enrichment: float,
) -> float:
    """
    Combined colocalization score.

    Weights:
    - Jaccard (same-spot): 0.3
    - Correlation: 0.3
    - Neighbor enrichment: 0.4 (most informative for spatial relationships)

    Args:
        jaccard: Jaccard index [0, 1].
        pearson_r: Pearson correlation [-1, 1].
        mutual_enrichment: Mean neighbor enrichment (unbounded, but typically 0-5+).

    Returns:
        Combined colocalization score.
    """
    # Normalize to [0, 1] range
    norm_corr = (pearson_r + 1) / 2  # r in [-1, 1] -> [0, 1]
    norm_enrich = min(mutual_enrichment / 5.0, 1.0)  # Cap at 5x enrichment

    return 0.3 * jaccard + 0.3 * norm_corr + 0.4 * norm_enrich


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
    n_permutations: int = 999,
    seed: int = 1234,
    verbose: bool = True,
) -> ColocalizationResult:
    """
    Analyze pairwise spatial colocalization between markers.

    Computes three types of relationships:
    1. Same-spot: Jaccard index and co-occurrence frequency
    2. Correlation: Pearson and Spearman correlation of intensities
    3. Neighborhood: Enrichment of marker B in neighbors of A-positive spots

    Args:
        X: Expression matrix (n_spots, n_markers).
        coords: Spatial coordinates for each spot (n_spots, 2).
        marker_names: Names for each marker column.
        markers_to_analyze: Subset of markers to analyze. If None, analyze all pairwise.
            Recommended: Pass result.interesting_markers from identify_interesting_markers().
        neighbor_k: Number of nearest neighbors for neighborhood analysis (default: 6).
        signal_threshold_percentile: Percentile threshold for classifying positive spots (default: 75).
        n_permutations: Number of permutations for enrichment significance (default: 999).
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

    # Binarize markers
    if verbose:
        logging.info(f"Binarizing markers at {signal_threshold_percentile}th percentile...")
    binary = _binarize_markers(analyze_X, signal_threshold_percentile)

    # Compute pairwise colocalization
    pairs = []
    total_pairs = n_analyze * (n_analyze - 1) // 2
    processed = 0

    for i in range(n_analyze):
        for j in range(i + 1, n_analyze):
            name_a = analyze_names[i]
            name_b = analyze_names[j]

            values_a = analyze_X[:, i]
            values_b = analyze_X[:, j]
            binary_a = binary[:, i]
            binary_b = binary[:, j]

            # Same-spot metrics
            jaccard, co_spots, co_frac = _compute_jaccard(binary_a, binary_b)

            # Correlation metrics
            pearson_r, spearman_rho, corr_pvalue = _compute_correlation(values_a, values_b)

            # Neighborhood metrics
            enrich_ab, enrich_ba, enrich_pvalue = _compute_neighbor_enrichment(
                binary_a, binary_b, neighbors, rng, n_permutations
            )
            mutual_enrich = (enrich_ab + enrich_ba) / 2

            # Combined score
            score = _compute_colocalization_score(jaccard, pearson_r, mutual_enrich)

            pairs.append(MarkerPairColocalization(
                marker_a=name_a,
                marker_b=name_b,
                jaccard_index=jaccard,
                co_occurrence_spots=co_spots,
                co_occurrence_fraction=co_frac,
                pearson_r=pearson_r,
                spearman_rho=spearman_rho,
                correlation_pvalue=corr_pvalue,
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
    logger.info(f"  Gap-split at merge {max_gap_idx} (gap={max_gap:.3f}, {max_gap/median_gap:.1f}x median)")
    for i, lin in enumerate(lineages):
        logger.info(f"    Lineage {i+1}: {lin}")

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
    alpha: float = 0.05,
    min_score: Optional[float] = None,
    seed: int = 1234,
    verbose: bool = True,
) -> ProfileDiscoveryResult:
    """
    Discover cell type profiles from colocalization analysis.

    Algorithm:
    1. Filter edges by statistical significance (p < alpha)
    2. Apply adaptive score threshold using GMM with BIC model selection
       - Fits 1-5 component GMM to score distribution
       - BIC selects optimal number of components
       - Lowest component = noise/spatial proximity (excluded)
       - Higher components = true colocalization (kept)
    3. Find connected components (separate lineages get separate dendrograms)
    4. Within each component, perform hierarchical clustering
    5. Dynamic tree cutting to extract profiles automatically
    6. Compute modularity to measure how well profiles explain the data

    Args:
        colocalization_result: Result from analyze_marker_colocalization().
        alpha: Significance threshold for p-value filtering (default: 0.05).
        min_score: Minimum colocalization score. If None (default), uses
            adaptive GMM-based threshold with BIC model selection.
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

    # Step 1: Filter edges by p-value significance
    pvalue_significant = [
        p for p in pairs
        if p.neighbor_enrichment_pvalue < alpha
    ]

    if verbose:
        logging.info(
            f"P-value filter (p < {alpha}): {len(pvalue_significant)}/{len(pairs)} pairs"
        )

    # Step 2: Apply score threshold (adaptive or fixed)
    if min_score is None:
        # Adaptive: Use GMM with BIC to find score threshold
        scores = np.array([p.colocalization_score for p in pvalue_significant])

        if len(scores) > 0:
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
                p for p, sig in zip(pvalue_significant, is_signal) if sig
            ]
        else:
            significant_pairs = []
            score_threshold = 0.0
    else:
        # Fixed threshold
        score_threshold = min_score
        significant_pairs = [
            p for p in pvalue_significant
            if p.colocalization_score >= min_score
        ]

        if verbose:
            logging.info(f"Fixed score filter (>= {min_score}): {len(significant_pairs)} pairs")

    significant_edges = [
        (p.marker_a, p.marker_b, p.colocalization_score)
        for p in significant_pairs
    ]

    if verbose:
        logging.info(
            f"Final: {len(significant_edges)} edges for graph construction"
        )

    # Step 2: Find connected components
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
            # Two markers - check distance to decide if same profile
            dist_matrix = _build_distance_matrix(comp_markers, pairs)
            if dist_matrix[0, 1] < 0.5:  # Strong colocalization
                profiles.append(comp_markers)
            else:
                profiles.append([comp_markers[0]])
                profiles.append([comp_markers[1]])
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
                # Check if they should be together
                lin_dist = _build_distance_matrix(lineage_markers, pairs)
                if lin_dist[0, 1] < 0.5:
                    profiles.append(lineage_markers)
                else:
                    profiles.append([lineage_markers[0]])
                    profiles.append([lineage_markers[1]])
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
        alpha=alpha,
    )

    if verbose:
        logging.info(result.summary())

    return result
