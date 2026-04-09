"""
Spot-level enhancement for CITEgeist profile discovery.

This module provides functionality for:
1. Per-spot reconstruction error tracking
2. Tissue region discovery (without explicit segmentation)
3. Region-specific profile scoring

These enhancements help improve profile selection by considering spatial
heterogeneity in the tissue and identifying region-specific expression patterns.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform


logger = logging.getLogger(__name__)


@dataclass
class TissueRegion:
    """
    Represents a discovered tissue region.

    Attributes:
        region_id: Unique identifier for this region.
        spot_indices: List of spot indices belonging to this region.
        centroid_coords: Spatial centroid (x, y) of the region.
        characteristic_profiles: Profile indices that are characteristic of this region.
        mean_reconstruction_error: Average reconstruction error in this region.
        region_weight: Weight for this region in optimization (based on size/importance).
        profile_scores: Scores for each profile in this region (n_profiles,).
    """
    region_id: int
    spot_indices: List[int]
    centroid_coords: Tuple[float, float]
    characteristic_profiles: List[int] = field(default_factory=list)
    mean_reconstruction_error: float = 0.0
    region_weight: float = 1.0
    profile_scores: Optional[np.ndarray] = None


def compute_spot_reconstruction_errors(
    Z: np.ndarray,
    Y_proportions: np.ndarray,
    candidate_profiles: List[set],
    selected_indices: List[int],
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """
    Compute per-spot reconstruction error (SSE).

    error[s] = sum_m (Z[s,m] - predicted[s,m])^2

    where predicted[s,m] = sum_p (Y[s,p] * profile_has_marker[p,m])

    Args:
        Z: Standardized expression matrix (n_spots, n_markers).
        Y_proportions: Proportion matrix (n_spots, n_selected).
        candidate_profiles: List of all candidate profiles.
        selected_indices: Indices of selected profiles in candidate_profiles.

    Returns:
        spot_errors: Per-spot SSE (n_spots,).
        diagnostics: Dict with diagnostic info (total_error, mean_error, etc.).
    """
    n_spots, n_markers = Z.shape
    n_selected = len(selected_indices)

    if n_selected == 0:
        return np.zeros(n_spots), {"status": "no_profiles", "total_error": 0.0}

    # Build profile-marker matrix for selected profiles
    # A[p, m] = 1 if profile p contains marker m
    A = np.zeros((n_selected, n_markers))
    for i, p_idx in enumerate(selected_indices):
        profile = candidate_profiles[p_idx]
        for m in profile:
            if m < n_markers:
                A[i, m] = 1.0

    # Compute predicted expression: predicted[s, m] = sum_p Y[s, p] * A[p, m]
    # This is a matrix multiplication: predicted = Y @ A
    predicted = Y_proportions @ A

    # Compute per-spot SSE
    residuals = Z - predicted
    spot_errors = np.sum(residuals ** 2, axis=1)

    # Diagnostics
    total_error = np.sum(spot_errors)
    mean_error = np.mean(spot_errors)
    median_error = np.median(spot_errors)
    max_error = np.max(spot_errors)
    high_error_spots = np.where(spot_errors > np.percentile(spot_errors, 90))[0]

    diagnostics = {
        "status": "computed",
        "total_error": float(total_error),
        "mean_error": float(mean_error),
        "median_error": float(median_error),
        "max_error": float(max_error),
        "n_high_error_spots": len(high_error_spots),
        "high_error_spot_indices": high_error_spots.tolist(),
        "r2": float(1.0 - total_error / (np.sum(Z ** 2) + 1e-9)),
    }

    logger.info(
        f"Reconstruction: R²={diagnostics['r2']:.4f}, "
        f"mean_SSE={mean_error:.4f}, max_SSE={max_error:.4f}"
    )

    return spot_errors, diagnostics


def discover_tissue_regions(
    coords: np.ndarray,
    Y_proportions: np.ndarray,
    spot_errors: Optional[np.ndarray] = None,
    n_regions: Optional[int] = None,
    spatial_weight: float = 0.5,
    min_region_size: int = 5,
    max_regions: int = 20,
) -> List[TissueRegion]:
    """
    Discover tissue regions via hierarchical clustering.

    Clustering is performed on a weighted combination of:
    - Spatial proximity (coordinates)
    - Profile proportion similarity

    Args:
        coords: Spatial coordinates (n_spots, 2).
        Y_proportions: Profile proportions (n_spots, n_profiles).
        spot_errors: Optional per-spot errors (n_spots,) for weighting.
        n_regions: Number of regions to discover. If None, auto-determined.
        spatial_weight: Weight for spatial vs profile similarity (0-1).
            0 = purely profile-based, 1 = purely spatial.
        min_region_size: Minimum number of spots per region.
        max_regions: Maximum number of regions.

    Returns:
        List of TissueRegion objects.
    """
    n_spots = coords.shape[0]

    if n_spots < min_region_size:
        # Single region if too few spots
        return [TissueRegion(
            region_id=0,
            spot_indices=list(range(n_spots)),
            centroid_coords=(float(np.mean(coords[:, 0])), float(np.mean(coords[:, 1]))),
            region_weight=1.0,
        )]

    # Normalize coordinates for distance computation
    coords_norm = (coords - coords.min(axis=0)) / (coords.max(axis=0) - coords.min(axis=0) + 1e-9)

    # Normalize proportions
    Y_norm = Y_proportions / (Y_proportions.max(axis=0, keepdims=True) + 1e-9)

    # Combine spatial and profile features
    profile_weight = 1.0 - spatial_weight
    combined_features = np.hstack([
        coords_norm * spatial_weight,
        Y_norm * profile_weight
    ])

    # Compute pairwise distances
    distances = pdist(combined_features, metric='euclidean')

    # Hierarchical clustering
    linkage_matrix = linkage(distances, method='ward')

    # Determine number of clusters
    if n_regions is None:
        # Auto-determine based on silhouette or gap statistic
        # Simple heuristic: sqrt(n_spots) with bounds
        n_regions = max(2, min(max_regions, int(np.sqrt(n_spots) / 2)))

    n_regions = min(n_regions, n_spots // min_region_size)
    n_regions = max(2, n_regions)

    # Cut tree to get cluster assignments
    cluster_labels = fcluster(linkage_matrix, n_regions, criterion='maxclust')

    # Build TissueRegion objects
    regions = []
    for region_id in range(1, n_regions + 1):
        spot_mask = cluster_labels == region_id
        spot_indices = np.where(spot_mask)[0].tolist()

        if len(spot_indices) < min_region_size:
            continue

        # Compute region centroid
        region_coords = coords[spot_mask]
        centroid = (float(np.mean(region_coords[:, 0])), float(np.mean(region_coords[:, 1])))

        # Compute region weight (based on size)
        region_weight = len(spot_indices) / n_spots

        # Compute mean reconstruction error if available
        mean_error = 0.0
        if spot_errors is not None:
            mean_error = float(np.mean(spot_errors[spot_mask]))

        region = TissueRegion(
            region_id=region_id,
            spot_indices=spot_indices,
            centroid_coords=centroid,
            mean_reconstruction_error=mean_error,
            region_weight=region_weight,
        )
        regions.append(region)

    logger.info(f"Discovered {len(regions)} tissue regions")
    for region in regions:
        logger.info(
            f"  Region {region.region_id}: {len(region.spot_indices)} spots, "
            f"error={region.mean_reconstruction_error:.4f}"
        )

    return regions


def compute_region_specific_profile_scores(
    Z: np.ndarray,
    candidate_profiles: List[set],
    regions: List[TissueRegion],
) -> Dict[int, np.ndarray]:
    """
    For each region, score how well each profile explains that region.

    Profile score = variance explained by the profile's markers in that region.

    Args:
        Z: Standardized expression matrix (n_spots, n_markers).
        candidate_profiles: List of candidate profiles.
        regions: List of TissueRegion objects.

    Returns:
        Dict mapping region_id to profile_scores array (n_candidates,).
    """
    n_candidates = len(candidate_profiles)
    region_scores = {}

    for region in regions:
        spot_indices = region.spot_indices
        if len(spot_indices) == 0:
            region_scores[region.region_id] = np.zeros(n_candidates)
            continue

        Z_region = Z[spot_indices, :]
        scores = np.zeros(n_candidates)

        for p, profile in enumerate(candidate_profiles):
            markers = list(profile)
            if len(markers) == 0:
                continue

            # Profile score = variance explained in this region
            marker_data = Z_region[:, markers]
            profile_expr = marker_data.mean(axis=1)
            predicted = np.outer(profile_expr, np.ones(len(markers)))

            before_error = np.sum(marker_data ** 2)
            after_error = np.sum((marker_data - predicted) ** 2)

            if before_error > 1e-9:
                scores[p] = max(0.0, (before_error - after_error) / before_error)

        region_scores[region.region_id] = scores
        region.profile_scores = scores

        # Identify characteristic profiles for this region
        top_profiles = np.argsort(scores)[-5:][::-1]  # Top 5
        region.characteristic_profiles = top_profiles[scores[top_profiles] > 0.1].tolist()

    return region_scores


def compute_global_profile_region_fit(
    region_scores: Dict[int, np.ndarray],
    regions: List[TissueRegion],
    aggregation: str = "max",
) -> np.ndarray:
    """
    Compute global profile-region fit scores.

    Each profile gets a score based on how well it fits across regions.

    Args:
        region_scores: Dict from compute_region_specific_profile_scores.
        regions: List of TissueRegion objects.
        aggregation: How to combine region scores - "max", "mean", "weighted".
            - "max": Best fit across any region
            - "mean": Average fit across regions
            - "weighted": Weighted average by region size

    Returns:
        profile_region_fit: (n_candidates,) array of fit scores in [0, 1].
    """
    if not regions or not region_scores:
        return np.array([])

    n_candidates = len(list(region_scores.values())[0])
    scores_matrix = np.zeros((len(regions), n_candidates))
    weights = np.zeros(len(regions))

    for i, region in enumerate(regions):
        scores_matrix[i] = region_scores[region.region_id]
        weights[i] = region.region_weight

    if aggregation == "max":
        return np.max(scores_matrix, axis=0)
    elif aggregation == "mean":
        return np.mean(scores_matrix, axis=0)
    elif aggregation == "weighted":
        weights = weights / (weights.sum() + 1e-9)
        return np.average(scores_matrix, axis=0, weights=weights)
    else:
        raise ValueError(f"Unknown aggregation method: {aggregation}")


def compute_spot_weights_from_errors(
    spot_errors: np.ndarray,
    method: str = "inverse",
    temperature: float = 1.0,
) -> np.ndarray:
    """
    Compute spot weights based on reconstruction errors.

    High-error spots may need more attention in profile selection.

    Args:
        spot_errors: Per-spot reconstruction errors (n_spots,).
        method: Weighting method:
            - "inverse": 1 / (1 + error) - downweight high-error spots
            - "proportional": error / sum(error) - upweight high-error spots
            - "softmax": softmax(error * temperature) - upweight high-error spots
        temperature: Temperature for softmax method.

    Returns:
        spot_weights: (n_spots,) array of weights summing to 1.
    """
    if method == "inverse":
        weights = 1.0 / (1.0 + spot_errors)
    elif method == "proportional":
        weights = spot_errors / (spot_errors.sum() + 1e-9)
    elif method == "softmax":
        scaled = spot_errors * temperature
        scaled = scaled - scaled.max()  # Numerical stability
        exp_scaled = np.exp(scaled)
        weights = exp_scaled / (exp_scaled.sum() + 1e-9)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Normalize to sum to 1
    weights = weights / (weights.sum() + 1e-9)
    return weights


def identify_poorly_explained_regions(
    regions: List[TissueRegion],
    error_percentile: float = 75,
) -> List[TissueRegion]:
    """
    Identify regions with high reconstruction error.

    These regions may benefit from additional or different profiles.

    Args:
        regions: List of TissueRegion objects.
        error_percentile: Percentile threshold for "high" error.

    Returns:
        List of regions with error above threshold.
    """
    if not regions:
        return []

    errors = np.array([r.mean_reconstruction_error for r in regions])
    threshold = np.percentile(errors, error_percentile)

    poorly_explained = [r for r in regions if r.mean_reconstruction_error >= threshold]

    logger.info(
        f"Identified {len(poorly_explained)} poorly explained regions "
        f"(error >= {threshold:.4f})"
    )

    return poorly_explained


def suggest_profiles_for_region(
    region: TissueRegion,
    candidate_profiles: List[set],
    current_selected: List[int],
    region_scores: np.ndarray,
    top_k: int = 3,
) -> List[Tuple[int, float]]:
    """
    Suggest additional profiles that might help explain a region.

    Args:
        region: TissueRegion to analyze.
        candidate_profiles: List of all candidate profiles.
        current_selected: Indices of currently selected profiles.
        region_scores: Profile scores for this region.
        top_k: Number of suggestions.

    Returns:
        List of (profile_index, score) tuples for suggested profiles.
    """
    # Find unselected profiles with high scores for this region
    unselected_mask = np.ones(len(candidate_profiles), dtype=bool)
    for idx in current_selected:
        unselected_mask[idx] = False

    # Rank unselected profiles by region score
    unselected_scores = region_scores * unselected_mask
    top_indices = np.argsort(unselected_scores)[-top_k:][::-1]

    suggestions = [
        (int(idx), float(unselected_scores[idx]))
        for idx in top_indices
        if unselected_scores[idx] > 0.01
    ]

    return suggestions


def run_spot_level_enhancement(
    Z: np.ndarray,
    coords: np.ndarray,
    candidate_profiles: List[set],
    Y_proportions: np.ndarray,
    selected_indices: List[int],
    n_regions: Optional[int] = None,
    spatial_weight: float = 0.5,
) -> Dict[str, Any]:
    """
    Run full spot-level enhancement pipeline.

    This is a convenience function that runs all spot-level analyses.

    Args:
        Z: Standardized expression matrix (n_spots, n_markers).
        coords: Spatial coordinates (n_spots, 2).
        candidate_profiles: List of candidate profiles.
        Y_proportions: Proportion matrix (n_spots, n_selected).
        selected_indices: Indices of selected profiles.
        n_regions: Number of regions (auto if None).
        spatial_weight: Weight for spatial clustering.

    Returns:
        Dict with all spot-level enhancement results.
    """
    results = {}

    # 1. Compute reconstruction errors
    spot_errors, error_diagnostics = compute_spot_reconstruction_errors(
        Z, Y_proportions, candidate_profiles, selected_indices
    )
    results["spot_errors"] = spot_errors
    results["error_diagnostics"] = error_diagnostics

    # 2. Discover tissue regions
    regions = discover_tissue_regions(
        coords, Y_proportions, spot_errors,
        n_regions=n_regions, spatial_weight=spatial_weight
    )
    results["regions"] = regions

    # 3. Compute region-specific profile scores
    region_scores = compute_region_specific_profile_scores(
        Z, candidate_profiles, regions
    )
    results["region_scores"] = region_scores

    # 4. Compute global profile-region fit
    profile_region_fit = compute_global_profile_region_fit(
        region_scores, regions, aggregation="max"
    )
    results["profile_region_fit"] = profile_region_fit

    # 5. Identify poorly explained regions
    poorly_explained = identify_poorly_explained_regions(regions)
    results["poorly_explained_regions"] = poorly_explained

    # 6. Suggest profiles for poorly explained regions
    suggestions = {}
    for region in poorly_explained:
        if region.profile_scores is not None:
            suggestions[region.region_id] = suggest_profiles_for_region(
                region, candidate_profiles, selected_indices,
                region.profile_scores, top_k=3
            )
    results["profile_suggestions"] = suggestions

    # 7. Compute spot weights
    spot_weights = compute_spot_weights_from_errors(spot_errors, method="inverse")
    results["spot_weights"] = spot_weights

    logger.info("Spot-level enhancement complete")
    logger.info(f"  Global R²: {error_diagnostics['r2']:.4f}")
    logger.info(f"  Regions discovered: {len(regions)}")
    logger.info(f"  Poorly explained regions: {len(poorly_explained)}")

    return results
