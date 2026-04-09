from __future__ import annotations

from typing import Any, Mapping


def should_enrich_module3_5_output(summary: Mapping[str, Any]) -> bool:
    return bool(summary.get("benchmark_passed")) and bool(summary.get("validated_pairs"))


should_enrich_single_cell_output = should_enrich_module3_5_output


def normalized_validated_pairs(summary: Mapping[str, Any]) -> set[str]:
    validated_pairs = summary.get("validated_pairs", [])
    if isinstance(validated_pairs, set):
        return validated_pairs
    return {str(pair_id) for pair_id in validated_pairs}
