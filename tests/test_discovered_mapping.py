import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Benchmarking" / "xenium_benchmarking" / "dictionary_discovery"))
from discovered_vs_hand import aggregate_to_types, leaves_of, map_leaves, nanaware_per_type, to_json_safe  # noqa: E402

pytestmark = pytest.mark.unit

PROFILE = {"A": {"Major": ["a1", "a2"]}, "B": {"Major": ["b1", "b2"]}, "C": {"Major": ["a1", "c1"]}}


def test_jaccard_pick():
    # {a1,a2} vs A: 1.0 ; vs C: 1/3 -> A
    assert map_leaves({"l": ["a1", "a2"]}, PROFILE) == {"l": "A"}


def test_tie_unmapped():
    # {a1} vs A: 1/2, vs C: 1/2 -> tie
    assert map_leaves({"l": ["a1"]}, PROFILE) == {"l": None}


def test_zero_overlap_unmapped():
    assert map_leaves({"l": ["z"]}, PROFILE) == {"l": None}


def test_suffix_insensitive():
    assert map_leaves({"l": ["b1-1", "b2-1"]}, PROFILE) == {"l": "B"}


def test_child_uses_parent_union_child():
    from CITEgeist.model.discovery.dictionary_discovery import DiscoveredDictionary

    d = DiscoveredDictionary(parents={"P1": ["a1", "a2"], "P2": ["b1"]}, children={"P1": {"x": ["c1"]}})
    lv = leaves_of(d)
    assert lv == {"P1/x": ["a1", "a2", "c1"], "P2": ["b1"]}


def test_aggregate_renormalizes_and_counts():
    lf = pd.DataFrame({"l1": [0.5, 0.0, 0.2], "l2": [0.3, 0.0, 0.2], "l3": [0.2, 1.0, 0.6]})
    mp = {"l1": "A", "l2": "A", "l3": None}
    df, um, nz = aggregate_to_types(lf, mp, ["A", "B"])
    assert np.allclose(df.iloc[0], [1.0, 0.0])
    assert np.allclose(df.iloc[1], [0.5, 0.5])  # zero mapped mass -> uniform
    assert nz == 1
    assert np.allclose(um, [0.2, 1.0, 0.6])
    assert np.allclose(df.sum(axis=1), 1.0)


def test_nanaware_and_json_safe():
    regs = [{"A_pearson_r": 0.5}, {"A_pearson_r": float("nan")}, {"A_pearson_r": 0.7}]
    m, n = nanaware_per_type(regs, ["A", "B"])
    assert abs(m["A"] - 0.6) < 1e-12 and n["A"] == 2 and m["B"] is None and n["B"] == 0
    assert to_json_safe({"x": np.float64("nan"), "y": [np.int64(2)]}) == {"x": None, "y": [2]}
