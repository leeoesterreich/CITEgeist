"""Coverage gap analysis for Module 3 and Module 3.5 inputs.

Compares Module 1 (marker interest) and Module 2 (colocalization) outputs
against the user's current cell_profile_dict and functional_marker_table,
flagging strongly-detected markers and pairs that have no coverage.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Optional

import pandas as pd

if TYPE_CHECKING:
    from CITEgeist.model.discovery.marker_interest import MarkerInterestResult
    from CITEgeist.model.discovery.spatial_colocalization import ColocalizationResult, ProfileDiscoveryResult

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CoverageCheckResult:
    """Coverage gap analysis results.

    Attributes:
        uncovered_markers: Interesting M1 markers not in any input dict.
            Columns: marker, morans_i, gmm_snr, kurtosis, interest_score
        uncovered_pairs: Strong M2 pairs where neither marker is covered.
            Columns: marker_a, marker_b, colocalization_score, suggested_type
        warning_lines: Human-readable warning strings (one per gap item).
        n_warnings: len(uncovered_markers) + len(uncovered_pairs).
        threshold_used: colocalization_threshold value used.
    """

    uncovered_markers: pd.DataFrame
    uncovered_pairs: pd.DataFrame
    warning_lines: list[str]
    n_warnings: int
    threshold_used: float

    def to_csv(self, output_dir: str) -> None:
        """Write coverage_check_markers.csv and coverage_check_pairs.csv to output_dir."""
        out = Path(output_dir)
        self.uncovered_markers.to_csv(out / "coverage_check_markers.csv", index=False)
        self.uncovered_pairs.to_csv(out / "coverage_check_pairs.csv", index=False)


def _build_covered_set(cell_profile_dict: dict, functional_marker_table: dict) -> set[str]:
    """Return set of all marker names present in either input dict."""
    covered: set[str] = set()
    for type_config in cell_profile_dict.values():
        if isinstance(type_config, dict):
            for markers in type_config.values():
                if isinstance(markers, list):
                    for item in markers:
                        covered.add(item[0] if isinstance(item, tuple) else str(item))
                elif isinstance(markers, str):
                    covered.add(markers)
        elif isinstance(type_config, list):
            for item in type_config:
                covered.add(item[0] if isinstance(item, tuple) else str(item))
    # functional_marker_table keys are marker names
    covered.update(str(m) for m in functional_marker_table)
    return covered


def check_module_coverage(
    m1_result: Optional["MarkerInterestResult"],
    m2_result: Optional["ColocalizationResult"],
    cell_profile_dict: dict,
    functional_marker_table: dict,
    profile_discovery_result: Optional["ProfileDiscoveryResult"] = None,
    colocalization_threshold: float = 0.75,
) -> "CoverageCheckResult":
    """Check for coverage gaps between M1/M2 outputs and Module 3/3.5 input dicts.

    Args:
        m1_result: Module 1 MarkerInterestResult. Pass None to skip M1 check.
        m2_result: Module 2 ColocalizationResult. Pass None to skip M2 check.
        cell_profile_dict: Current Module 3 cell-type-to-marker mapping.
        functional_marker_table: Current Module 3.5 functional marker table.
        profile_discovery_result: Optional Module 2b ProfileDiscoveryResult.
            When provided, used to populate suggested_type for uncovered pairs.
        colocalization_threshold: Pairs with colocalization_score above this
            value are treated as "clearly called". Default 0.75 (top quartile).

    Returns:
        CoverageCheckResult with uncovered_markers, uncovered_pairs, warning_lines.
    """
    covered = _build_covered_set(cell_profile_dict, functional_marker_table)
    warning_lines: list[str] = []

    # M1 gap analysis: interesting markers not in covered set
    uncovered_marker_rows: list[dict] = []
    if m1_result is not None:
        for m in m1_result.markers:
            if (m.passed_kurtosis or m.passed_morans) and m.name not in covered:
                uncovered_marker_rows.append({
                    "marker": m.name,
                    "morans_i": m.morans_i,
                    "gmm_snr": m.gmm_snr,
                    "kurtosis": m.kurtosis,
                    "interest_score": m.interest_score,
                })
                warning_lines.append(
                    f"COVERAGE WARNING [M1]: interesting marker '{m.name}' "
                    f"(morans_i={m.morans_i:.3f}, gmm_snr={m.gmm_snr:.2f}) "
                    "is not covered by cell_profile_dict or functional_marker_table"
                )

    uncovered_markers_df = (
        pd.DataFrame(uncovered_marker_rows, columns=["marker", "morans_i", "gmm_snr", "kurtosis", "interest_score"])
        .sort_values("interest_score", ascending=False)
        .reset_index(drop=True)
    )

    # M2 gap analysis: strong pairs where neither marker is covered
    uncovered_pair_rows: list[dict] = []
    if m2_result is not None:
        for pair in m2_result.pairs:
            if pair.colocalization_score < colocalization_threshold:
                continue
            if pair.marker_a in covered or pair.marker_b in covered:
                continue
            suggested_type = ""
            if profile_discovery_result is not None:
                profile = profile_discovery_result.get_profile_for_marker(pair.marker_a)
                if profile is None:
                    profile = profile_discovery_result.get_profile_for_marker(pair.marker_b)
                if profile is not None:
                    idx = profile_discovery_result.profiles.index(profile)
                    suggested_type = f"profile_{idx}"
            uncovered_pair_rows.append({
                "marker_a": pair.marker_a,
                "marker_b": pair.marker_b,
                "colocalization_score": pair.colocalization_score,
                "suggested_type": suggested_type,
            })
            warning_lines.append(
                f"COVERAGE WARNING [M2]: strong pair ('{pair.marker_a}', '{pair.marker_b}', "
                f"score={pair.colocalization_score:.3f}) is not covered by either input dict"
            )

    uncovered_pairs_df = (
        pd.DataFrame(uncovered_pair_rows, columns=["marker_a", "marker_b", "colocalization_score", "suggested_type"])
        .sort_values("colocalization_score", ascending=False)
        .reset_index(drop=True)
    )

    n_warnings = len(uncovered_markers_df) + len(uncovered_pairs_df)
    logger.debug("Coverage check: %d warnings (threshold=%.2f)", n_warnings, colocalization_threshold)
    return CoverageCheckResult(
        uncovered_markers=uncovered_markers_df,
        uncovered_pairs=uncovered_pairs_df,
        warning_lines=warning_lines,
        n_warnings=n_warnings,
        threshold_used=colocalization_threshold,
    )
