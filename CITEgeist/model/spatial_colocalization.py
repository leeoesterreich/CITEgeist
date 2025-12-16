"""
Spatial colocalization analysis for marker pairs.

Analyzes three types of spatial relationships:
1. Same-spot co-occurrence: Do markers appear together in the same spots?
2. Expression correlation: Are marker intensities correlated?
3. Adjacent-spot enrichment: Are markers enriched in neighboring spots?

This module produces informational reports only - no automatic profile selection.
"""

import logging
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.spatial import cKDTree
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

    # Filter to markers of interest
    if markers_to_analyze is not None:
        marker_indices = [
            i for i, m in enumerate(marker_names)
            if m in markers_to_analyze
        ]
        if len(marker_indices) == 0:
            raise ValueError("No markers from markers_to_analyze found in marker_names")
        analyze_names = [marker_names[i] for i in marker_indices]
        analyze_X = X[:, marker_indices]
    else:
        marker_indices = list(range(n_markers))
        analyze_names = marker_names
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
