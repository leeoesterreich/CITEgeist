"""Module 2 proposal builder: rank candidate cell-type profiles into Module 3 (lineage)
and Module 3.5 (functional-marker) candidate tables for review."""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class MarkerRoleConfig:
    """Marker role assignments used to rank profiles.

    Attributes:
        lineage_anchors: marker -> parent cell-type name (lineage discriminators).
        functional_markers: markers treated as functional (Module 3.5) signals.
        ambiguous_markers: markers that penalize a profile's rank score.
        excluded_markers: markers dropped from every profile before ranking.
    """

    lineage_anchors: dict[str, str]
    functional_markers: set[str]
    ambiguous_markers: set[str]
    excluded_markers: set[str]


def _profile_rank_score(profile: dict, role_config: MarkerRoleConfig) -> tuple[float, int]:
    markers = [marker for marker in profile["markers"] if marker not in role_config.excluded_markers]
    ambiguity_count = sum(1 for marker in markers if marker in role_config.ambiguous_markers)
    rank_score = max(float(profile["selection_score"]) - 0.1 * ambiguity_count, 0.0)
    return rank_score, ambiguity_count


def build_candidate_rank_lists(
    profiles: list[dict],
    role_config: MarkerRoleConfig,
    ontology_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank profiles into Module 3 and Module 3.5 candidate tables.

    Args:
        profiles: list of profile dicts with keys ``profile_id``, ``markers``, ``selection_score``.
        role_config: MarkerRoleConfig assigning lineage/functional/ambiguous/excluded roles.
        ontology_name: label stamped onto every emitted row.

    Returns:
        (module3_df, module3_5_df): two DataFrames sorted by descending rank_score;
        module3_df has one row per profile with a lineage anchor, module3_5_df has one
        row per (profile, functional marker) pair.
    """
    module3_rows: list[dict[str, object]] = []
    module3_5_rows: list[dict[str, object]] = []

    for profile in profiles:
        markers = [marker for marker in profile["markers"] if marker not in role_config.excluded_markers]
        lineage_markers = [marker for marker in markers if marker in role_config.lineage_anchors]
        functional_markers = [marker for marker in markers if marker in role_config.functional_markers]

        if not lineage_markers:
            continue

        parent_type = role_config.lineage_anchors[lineage_markers[0]]
        rank_score, ambiguity_count = _profile_rank_score(profile, role_config)

        module3_rows.append(
            {
                "candidate_type": "module3",
                "profile_id": profile["profile_id"],
                "parent_type": parent_type,
                "markers": ",".join(markers),
                "rank_score": rank_score,
                "ambiguity_count": ambiguity_count,
                "ontology_name": ontology_name,
            }
        )

        for marker in functional_markers:
            module3_5_rows.append(
                {
                    "candidate_type": "module3_5",
                    "profile_id": profile["profile_id"],
                    "parent_type": parent_type,
                    "functional_marker": marker,
                    "pair_id": f"{parent_type}__{marker}",
                    "markers": ",".join(markers),
                    "rank_score": rank_score,
                    "ambiguity_count": ambiguity_count,
                    "ontology_name": ontology_name,
                }
            )

    module3_df = pd.DataFrame(module3_rows)
    module3_5_df = pd.DataFrame(module3_5_rows)

    if not module3_df.empty:
        module3_df = module3_df.sort_values("rank_score", ascending=False, kind="stable").reset_index(drop=True)
    if not module3_5_df.empty:
        module3_5_df = module3_5_df.sort_values("rank_score", ascending=False, kind="stable").reset_index(drop=True)

    return module3_df, module3_5_df
