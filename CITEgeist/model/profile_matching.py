"""
Profile matching utilities for comparing discovered profiles to ground truth.

This module provides functions to match auto-discovered cell type profiles
to known ground truth profiles by exact marker composition, and to remap
proportion DataFrames for proper benchmarking comparison.

Author: A. Chang
License: BSD 3-Clause
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Set

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ProfileMatchResult:
    """Container for profile matching results."""

    matched_profiles: Dict[str, str]
    """Mapping from discovered profile name to matched ground truth cell type."""

    unmatched_discovered: List[str]
    """Discovered profiles that don't match any ground truth."""

    missing_ground_truth: List[str]
    """Ground truth cell types with no matching discovered profile."""

    profile_recovery_rate: float
    """Fraction of GT cell types with exact marker match (n_matched / n_gt)."""

    marker_precision: float
    """Of all discovered markers, fraction that exist in ground truth."""

    marker_recall: float
    """Of all ground truth markers, fraction that were discovered."""

    false_discovery_rate: float
    """Fraction of discovered profiles with no ground truth match."""


def extract_markers_from_profile(profile_dict: Dict[str, List[str]]) -> Set[str]:
    """
    Extract the set of markers from a profile dictionary.

    Args:
        profile_dict: Profile definition, e.g., {"Major": ["Marker1", "Marker2"]}

    Returns:
        Set of marker names from the "Major" key.
    """
    return set(profile_dict.get("Major", []))


def match_profiles_to_ground_truth(
    discovered_profiles: Dict[str, Dict[str, List[str]]],
    ground_truth_profiles: Dict[str, Dict[str, List[str]]],
    allow_partial: bool = True,
) -> ProfileMatchResult:
    """
    Match discovered profiles to ground truth by marker composition.

    Matching strategy (in order of priority):
    1. EXACT match: discovered markers == ground truth markers
    2. PARTIAL match (if allow_partial=True): discovered markers ⊆ ground truth markers
       (i.e., all discovered markers exist in ground truth, but not all GT markers discovered)

    For partial matches, if multiple GT profiles could match, prefers:
    - Highest overlap ratio (discovered ∩ GT) / (GT markers)
    - In case of tie, first alphabetically

    Args:
        discovered_profiles: Auto-discovered profiles in CITEgeist format.
            Example: {"Marker1_Marker2": {"Major": ["Marker1", "Marker2"]}}
        ground_truth_profiles: Known cell type profiles.
            Example: {"B-cells": {"Major": ["B-cells_Protein_1", "B-cells_Protein_2"]}}
        allow_partial: If True, allow singleton/partial profiles to match GT.
            A partial match requires ALL discovered markers to be in the GT profile.
            Default: True.

    Returns:
        ProfileMatchResult with matching information and accuracy metrics.

    Examples:
        >>> discovered = {"M1": {"Major": ["M1"]}}  # singleton
        >>> ground_truth = {"CellType": {"Major": ["M1", "M2"]}}
        >>> result = match_profiles_to_ground_truth(discovered, ground_truth)
        >>> result.matched_profiles
        {"M1": "CellType"}  # partial match allowed
    """
    matched_profiles: Dict[str, str] = {}
    unmatched_discovered: List[str] = []
    matched_gt_types: Set[str] = set()

    # Build marker set lookup for ground truth
    gt_marker_sets: Dict[str, Set[str]] = {
        ct: extract_markers_from_profile(profile)
        for ct, profile in ground_truth_profiles.items()
    }

    # For each discovered profile, find best matching ground truth
    for disc_name, disc_profile in discovered_profiles.items():
        disc_markers = extract_markers_from_profile(disc_profile)

        best_match = None
        best_overlap_ratio = 0.0
        is_exact = False

        for gt_name, gt_markers in gt_marker_sets.items():
            # Check for exact match first
            if disc_markers == gt_markers:
                best_match = gt_name
                best_overlap_ratio = 1.0
                is_exact = True
                break  # Exact match is best possible

            # Check for partial match: all discovered markers must be in GT
            if allow_partial and disc_markers and disc_markers.issubset(gt_markers):
                overlap_ratio = len(disc_markers) / len(gt_markers)
                if overlap_ratio > best_overlap_ratio:
                    best_match = gt_name
                    best_overlap_ratio = overlap_ratio

        if best_match is not None:
            matched_profiles[disc_name] = best_match
            matched_gt_types.add(best_match)
            match_type = "exact" if is_exact else f"partial ({best_overlap_ratio:.0%})"
            logger.info(f"Matched profile '{disc_name}' -> '{best_match}' [{match_type}]")
        else:
            unmatched_discovered.append(disc_name)
            logger.warning(
                f"No ground truth match for discovered profile: {disc_name} "
                f"with markers {disc_markers}"
            )

    # Find ground truth types with no matching discovered profile
    missing_ground_truth = [
        gt_name
        for gt_name in ground_truth_profiles.keys()
        if gt_name not in matched_gt_types
    ]

    if missing_ground_truth:
        logger.warning(f"Ground truth cell types not discovered: {missing_ground_truth}")

    # Calculate metrics
    n_gt = len(ground_truth_profiles)
    n_discovered = len(discovered_profiles)
    n_matched_gt = len(matched_gt_types)  # Unique GT types matched

    # Recovery rate: fraction of GT cell types that have at least one matching profile
    profile_recovery_rate = n_matched_gt / n_gt if n_gt > 0 else 0.0
    # FDR: fraction of discovered profiles with no GT match
    false_discovery_rate = (
        len(unmatched_discovered) / n_discovered if n_discovered > 0 else 0.0
    )

    # Marker-level metrics
    all_discovered_markers: Set[str] = set()
    for profile in discovered_profiles.values():
        all_discovered_markers.update(extract_markers_from_profile(profile))

    all_gt_markers: Set[str] = set()
    for profile in ground_truth_profiles.values():
        all_gt_markers.update(extract_markers_from_profile(profile))

    correct_markers = all_discovered_markers & all_gt_markers

    marker_precision = (
        len(correct_markers) / len(all_discovered_markers)
        if all_discovered_markers
        else 0.0
    )
    marker_recall = (
        len(correct_markers) / len(all_gt_markers) if all_gt_markers else 0.0
    )

    return ProfileMatchResult(
        matched_profiles=matched_profiles,
        unmatched_discovered=unmatched_discovered,
        missing_ground_truth=missing_ground_truth,
        profile_recovery_rate=profile_recovery_rate,
        marker_precision=marker_precision,
        marker_recall=marker_recall,
        false_discovery_rate=false_discovery_rate,
    )


def create_remapped_proportions(
    original_proportions: pd.DataFrame,
    match_result: ProfileMatchResult,
    ground_truth_cell_types: List[str],
) -> pd.DataFrame:
    """
    Remap discovered profile names to matched ground truth cell type names.

    For unmatched discovered profiles, aggregate into "Unknown_Discovered" category.
    For missing ground truth types, add columns with zeros.

    Args:
        original_proportions: DataFrame with discovered profile names as columns.
            Shape: (n_spots, n_discovered_profiles)
        match_result: Profile matching result from match_profiles_to_ground_truth.
        ground_truth_cell_types: List of ground truth cell type names in desired order.

    Returns:
        DataFrame with columns renamed to ground truth cell types.
        Unmatched profiles aggregated into "Unknown_Discovered".
        Missing GT types have zero columns.

    Examples:
        >>> props = pd.DataFrame({"M1_M2": [0.5, 0.3], "M3": [0.5, 0.7]})
        >>> match_result.matched_profiles = {"M1_M2": "B-cells"}
        >>> match_result.unmatched_discovered = ["M3"]
        >>> result = create_remapped_proportions(props, match_result, ["B-cells", "T-cells"])
        >>> result.columns.tolist()
        ["B-cells", "T-cells", "Unknown_Discovered"]
    """
    result = pd.DataFrame(index=original_proportions.index)

    # Add matched columns with renamed names
    for disc_name, gt_name in match_result.matched_profiles.items():
        if disc_name in original_proportions.columns:
            result[gt_name] = original_proportions[disc_name]

    # Aggregate unmatched discovered profiles into "Unknown_Discovered"
    unmatched_cols = [
        col
        for col in match_result.unmatched_discovered
        if col in original_proportions.columns
    ]
    if unmatched_cols:
        result["Unknown_Discovered"] = original_proportions[unmatched_cols].sum(axis=1)
        logger.info(
            f"Aggregated {len(unmatched_cols)} unmatched profiles into 'Unknown_Discovered'"
        )

    # Add zeros for missing ground truth types
    for gt_name in match_result.missing_ground_truth:
        if gt_name not in result.columns:
            result[gt_name] = 0.0
            logger.info(f"Added zero column for missing GT type: {gt_name}")

    # Ensure all ground truth columns exist (even if not in missing list)
    for gt_name in ground_truth_cell_types:
        if gt_name not in result.columns:
            result[gt_name] = 0.0

    # Reorder to match ground truth order, then add Unknown_Discovered at end
    final_columns = [c for c in ground_truth_cell_types if c in result.columns]
    if "Unknown_Discovered" in result.columns:
        final_columns.append("Unknown_Discovered")

    return result[final_columns]


def benchmark_profile_discovery(
    discovered_profiles: Dict[str, Dict[str, List[str]]],
    ground_truth_profiles: Dict[str, Dict[str, List[str]]],
) -> Dict[str, Any]:
    """
    Calculate comprehensive metrics for profile discovery accuracy.

    This is a convenience function that wraps match_profiles_to_ground_truth
    and returns a flat dictionary suitable for DataFrame conversion and CSV export.

    Args:
        discovered_profiles: Auto-discovered profiles from discover_profiles().
        ground_truth_profiles: Known/expected cell type profiles.

    Returns:
        Dictionary with discovery accuracy metrics:
        - profile_recovery_rate: Fraction of GT cell types exactly recovered
        - marker_precision: Of all discovered markers, fraction in GT
        - marker_recall: Of all GT markers, fraction discovered
        - false_discovery_rate: Fraction of discovered profiles with no GT match
        - n_discovered: Number of discovered profiles
        - n_ground_truth: Number of ground truth cell types
        - n_matched: Number of exactly matched profiles
        - matched_pairs: List of (discovered_name, gt_name) tuples
        - unmatched_discovered: List of unmatched discovered profile names
        - missing_ground_truth: List of GT cell types not discovered
    """
    match_result = match_profiles_to_ground_truth(
        discovered_profiles, ground_truth_profiles
    )

    return {
        "profile_recovery_rate": match_result.profile_recovery_rate,
        "marker_precision": match_result.marker_precision,
        "marker_recall": match_result.marker_recall,
        "false_discovery_rate": match_result.false_discovery_rate,
        "n_discovered": len(discovered_profiles),
        "n_ground_truth": len(ground_truth_profiles),
        "n_matched": len(match_result.matched_profiles),
        "matched_pairs": list(match_result.matched_profiles.items()),
        "unmatched_discovered": match_result.unmatched_discovered,
        "missing_ground_truth": match_result.missing_ground_truth,
    }
