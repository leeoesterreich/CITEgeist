"""Tests for gt_rebuild/stability.py (synthetic)."""

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import stability  # noqa: E402

pytestmark = pytest.mark.unit


def test_agreement_identical_is_1_and_one_flip_is_n_minus_1_over_n():
    a = pd.Series(["x", "y", "x", "y"], index=list("abcd"))
    assert stability.label_agreement(a, a.copy()) == 1.0
    b = a.copy()
    b["a"] = "y"
    assert stability.label_agreement(a, b) == pytest.approx(3 / 4)


def test_id_mismatch_raises():
    a = pd.Series(["x", "y"], index=["a", "b"])
    b = pd.Series(["x", "y"], index=["a", "c"])
    with pytest.raises(ValueError):
        stability.label_agreement(a, b)


def test_spot_change_hand_computed():
    pa = pd.DataFrame({"T": [0.5, 0.5]}, index=["s1", "s2"])
    pb = pd.DataFrame({"T": [0.6, 0.2]}, index=["s1", "s2"])
    out = stability.spot_prop_change(pa, pb)
    assert out.loc["T", "mean_abs_change"] == pytest.approx(0.2)
    assert out.loc["T", "max_abs_change"] == pytest.approx(0.3)


def test_summary_reports_each_variant():
    base = pd.Series(["x", "y"], index=["a", "b"])
    out = stability.stability_summary({"leiden_0.6": base, "seed_1": base.replace("y", "x")}, base)
    assert out == {"leiden_0.6": 1.0, "seed_1": 0.5}
