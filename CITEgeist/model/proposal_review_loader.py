from __future__ import annotations

from collections import OrderedDict

import pandas as pd


def _promoted_rows(reviewed: pd.DataFrame) -> pd.DataFrame:
    required = {"review_state"}
    missing = required - set(reviewed.columns)
    if missing:
        raise KeyError(f"Missing required review columns: {', '.join(sorted(missing))}")
    return reviewed.loc[reviewed["review_state"] == "promote"].copy()


def build_module3_profile_dict_from_review(reviewed: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    promoted = _promoted_rows(reviewed)
    promoted = promoted.sort_index(kind="stable")

    result: OrderedDict[str, dict[str, list[str]]] = OrderedDict()
    for _, row in promoted.iterrows():
        parent_type = str(row["parent_type"])
        markers = [marker for marker in str(row["markers"]).split(",") if marker]
        bucket = result.setdefault(parent_type, {"Major": [], "Minor": []})
        for marker in markers:
            if marker not in bucket["Major"]:
                bucket["Major"].append(marker)
    return dict(result)


def build_module3_5_table_from_review(reviewed: pd.DataFrame) -> dict[str, dict[str, list[str]]]:
    promoted = _promoted_rows(reviewed)
    promoted = promoted.sort_index(kind="stable")

    result: OrderedDict[str, dict[str, list[str]]] = OrderedDict()
    for _, row in promoted.iterrows():
        marker = str(row["functional_marker"])
        parent_type = str(row["parent_type"])
        bucket = result.setdefault(marker, {"function": "discovered", "active_types": []})  # type: ignore[dict-item]
        if parent_type not in bucket["active_types"]:
            bucket["active_types"].append(parent_type)
    return dict(result)
