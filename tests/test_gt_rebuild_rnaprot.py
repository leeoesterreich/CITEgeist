"""Tests for gt_rebuild/rnaprot.py (synthetic data only)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import rnaprot  # noqa: E402

pytestmark = pytest.mark.unit

LIN_COLS = ["CD20", "CD138", "CD3E", "CD68", "CD163", "CD31", "PanCK", "alphaSMA", "CD4", "CD8A", "E-Cadherin"]


def _z(rows, n=None):
    """z-frame with every channel at 0 except those given per cell: {cell: {chan: z}}."""
    df = pd.DataFrame(0.0, index=list(rows), columns=LIN_COLS)
    for cid, d in rows.items():
        for c, v in d.items():
            df.loc[cid, c] = v
    return df


def _labels(d6, d7=None):
    s6 = pd.Series(d6)
    if d7 is None:
        d7 = {k: v for k, v in d6.items()}
    return s6, pd.Series(d7)


def test_lineage_call_none_when_nothing_positive():
    z = _z({"a": {"CD3E": 2.9, "PanCK": 1.0}})
    assert rnaprot.lineage_call(z, Z=3.0)["a"] == "none"


def test_lineage_call_picks_highest_positive_lineage():
    z = _z({"a": {"CD3E": 5.0, "PanCK": 7.0}, "b": {"CD68": 4.0, "CD163": 6.0}, "c": {"CD138": 9.0, "CD20": 3.5}})
    c = rnaprot.lineage_call(z, Z=3.0)
    assert c["a"] == "Epithelial"
    assert c["b"] == "Macrophages"
    assert c["c"] == "B cells"


def test_z_boundary_equal_is_not_positive():
    z = _z({"a": {"CD3E": 3.0}, "b": {"CD3E": 3.0000001}})
    c = rnaprot.lineage_call(z, Z=3.0)
    assert c["a"] == "none"
    assert c["b"] == "T cells"


def test_ecadherin_is_unused():
    z = _z({"a": {"E-Cadherin": 50.0}})
    assert rnaprot.lineage_call(z)["a"] == "none"
    # and it never drops an RNA-labelled T cell
    s6, s7 = _labels({"a": "T cells"}, {"a": "CD4+ T cells"})
    out = rnaprot.apply_no_contradiction(s6, s7, z)
    assert out.loc["a", "keep6"] and out.loc["a", "keep7"]


def test_none_call_keeps_cell():
    z = _z({"a": {}})
    s6, s7 = _labels({"a": "Epithelial"})
    out = rnaprot.apply_no_contradiction(s6, s7, z)
    assert out.loc["a", ["keep6", "keep7"]].all()
    assert out.loc["a", "reason"] == ""


def test_other_lineage_call_drops_cell_with_reason():
    z = _z({"a": {"CD3E": 6.0}})
    s6, s7 = _labels({"a": "Epithelial"})
    out = rnaprot.apply_no_contradiction(s6, s7, z)
    assert not out.loc["a", "keep6"] and not out.loc["a", "keep7"]
    assert "T cells" in out.loc["a", "reason"]


def test_same_lineage_call_keeps_cell():
    z = _z({"a": {"PanCK": 6.0}})
    s6, s7 = _labels({"a": "Epithelial"})
    out = rnaprot.apply_no_contradiction(s6, s7, z)
    assert out.loc["a", ["keep6", "keep7"]].all()


def test_plasma_guard_b_cell_with_cd138_never_dropped():
    z = _z({"a": {"CD138": 5.0, "CD3E": 8.0}, "b": {"CD3E": 8.0}})
    s6, s7 = _labels({"a": "B cells", "b": "B cells"})
    out = rnaprot.apply_no_contradiction(s6, s7, z)
    assert out.loc["a", ["keep6", "keep7"]].all()
    assert not out.loc["b", "keep6"]


def test_plasma_guard_only_for_b_label():
    z = _z({"a": {"CD138": 5.0, "CD3E": 8.0}})
    s6, s7 = _labels({"a": "Epithelial"})
    out = rnaprot.apply_no_contradiction(s6, s7, z)
    assert not out.loc["a", "keep6"]


def test_plasma_guard_requires_cd138_strictly_above_z():
    z = _z({"a": {"CD138": 3.0, "CD3E": 8.0}})
    s6, s7 = _labels({"a": "B cells"})
    assert not rnaprot.apply_no_contradiction(s6, s7, z).loc["a", "keep6"]


def test_seven_type_subtype_contradiction():
    z = _z(
        {
            "cd4_but_cd8": {"CD3E": 6.0, "CD8A": 6.0, "CD4": 0.0},
            "cd8_but_cd4": {"CD3E": 6.0, "CD4": 6.0, "CD8A": 0.0},
            "cd4_both": {"CD3E": 6.0, "CD4": 6.0, "CD8A": 6.0},
            "cd4_neither": {"CD3E": 6.0},
            "cd4_agrees": {"CD3E": 6.0, "CD4": 6.0},
        }
    )
    d7 = {
        "cd4_but_cd8": "CD4+ T cells",
        "cd8_but_cd4": "CD8+ T cells",
        "cd4_both": "CD4+ T cells",
        "cd4_neither": "CD4+ T cells",
        "cd4_agrees": "CD4+ T cells",
    }
    d6 = {k: "T cells" for k in d7}
    out = rnaprot.apply_no_contradiction(pd.Series(d6), pd.Series(d7), z)
    assert out["keep6"].all()  # 6-type never looks at subtype
    assert not out.loc["cd4_but_cd8", "keep7"]
    assert not out.loc["cd8_but_cd4", "keep7"]
    assert out.loc[["cd4_both", "cd4_neither", "cd4_agrees"], "keep7"].all()
    assert "subtype" in out.loc["cd4_but_cd8", "reason"]


def test_lineage_drop_applies_to_7type_too():
    z = _z({"a": {"PanCK": 6.0}})
    out = rnaprot.apply_no_contradiction(pd.Series({"a": "T cells"}), pd.Series({"a": "CD8+ T cells"}), z)
    assert not out.loc["a", "keep6"] and not out.loc["a", "keep7"]


def test_label_indices_must_match():
    z = _z({"a": {}, "b": {}})
    with pytest.raises(ValueError):
        rnaprot.apply_no_contradiction(pd.Series({"a": "T cells"}), pd.Series({"b": "CD4+ T cells"}), z)


def test_negative_fit_uses_low_nonzero_component():
    rng = np.random.RandomState(1)
    low = rng.normal(1.0, 0.2, 5000)
    high = rng.normal(4.0, 0.3, 3000)
    x = np.concatenate([np.zeros(2000), low, high])
    m, sd = rnaprot.negative_fit(x, seed=0)
    assert abs(m - 1.0) < 0.05 and abs(sd - 0.2) < 0.05


def test_negative_fit_include_zeros_pulls_mean_down():
    rng = np.random.RandomState(1)
    x = np.concatenate([np.zeros(2000), rng.normal(1.0, 0.2, 5000), rng.normal(4.0, 0.3, 3000)])
    m0, _ = rnaprot.negative_fit(x, seed=0)
    m1, _ = rnaprot.negative_fit(x, seed=0, include_zeros=True)
    assert m1 < m0


def test_protein_z_is_per_channel_and_high_cells_are_positive():
    rng = np.random.RandomState(0)
    n = 4000
    bg = rng.poisson(2.0, n).astype(float)
    sig = bg.copy()
    sig[:200] += 500  # bright cells
    prot = pd.DataFrame({"X": sig, "Y": bg}, index=[f"c{i}" for i in range(n)])
    z = rnaprot.protein_z(prot)
    assert list(z.columns) == ["X", "Y"] and (z.index == prot.index).all()
    assert (z["X"].iloc[:200] > 3).all()
    assert (z["X"].iloc[200:] > 3).mean() < 0.05
    assert set(z.attrs["neg_fit"]) == {"X", "Y"}


def test_protein_z_keep_subset_but_fits_all_channels_in_order():
    rng = np.random.RandomState(0)
    prot = pd.DataFrame({c: rng.poisson(3, 1000).astype(float) for c in "ABC"})
    full = rnaprot.protein_z(prot)
    sub = rnaprot.protein_z(prot, keep=["C"])
    assert list(sub.columns) == ["C"]
    np.testing.assert_allclose(sub["C"], full["C"])
