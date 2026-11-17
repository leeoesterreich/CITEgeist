# tests/benchmarking/test_protein_perturb.py
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
from protein_sensitivity import perturb as P  # noqa: E402

pytestmark = pytest.mark.unit

PROFILE = {
    "B cells": {"Major": ["CD20"], "Minor": ["CD45RA"]},
    "Macrophages": {"Major": ["CD68", "CD163"], "Minor": ["CD16"]},
    "Epithelial": {"Major": ["PanCK"]},
}
# 7-type-like shared marker: CD3E is Major for both T types
SHARED = {
    "CD4+ T cells": {"Major": ["CD3E", "CD4"]},
    "CD8+ T cells": {"Major": ["CD3E", "CD8A"]},
}
SOFT = {"Myeloid": {"Major": ["CD68"], "Soft": [("CD14", 0.5)]}, "B cells": {"Major": ["CD20"]}}


def _X(seed=0):
    return np.random.default_rng(seed).poisson(20, size=(50, 4)).astype(np.int64)


def _markers(d):
    return {m for v in d.values() for m in v.get("Major", [])} | {s[0] for v in d.values() for s in v.get("Soft", [])}


def test_thin_keep_all_is_identity():
    X = _X()
    assert np.array_equal(P.thin_counts(X, 1.0, np.random.default_rng(1)), X)


def test_thin_reduces_mean_and_is_bounded():
    X = _X()
    Y = P.thin_counts(X, 0.25, np.random.default_rng(1))
    assert (Y <= X).all() and (Y >= 0).all()
    assert abs(Y.sum() / X.sum() - 0.25) < 0.05


def test_thin_reproducible_with_seed():
    X = _X()
    assert np.array_equal(
        P.thin_counts(X, 0.5, np.random.default_rng(7)), P.thin_counts(X, 0.5, np.random.default_rng(7))
    )


def test_thin_accepts_integral_float():
    assert P.thin_counts(_X().astype(float), 0.5, np.random.default_rng(0)).dtype == np.int64


def test_thin_rejects_nonintegral():
    with pytest.raises(ValueError, match="integer"):
        P.thin_counts(np.array([[1.5, 2.0]]), 0.5, np.random.default_rng(0))


def test_dropout_zeros_whole_channels_from_candidates_only():
    X = _X()
    names = ["CD20", "CD68", "CD163", "PanCK"]
    Y, dropped = P.dropout_channels(X, names, ["CD68", "CD163", "PanCK"], 0.5, np.random.default_rng(3))
    assert len(dropped) == 2 and "CD20" not in dropped
    for m in dropped:
        assert (Y[:, names.index(m)] == 0).all()
    assert np.array_equal(Y[:, 0], X[:, 0])


def test_dropout_survives_real_preprocess_antibody(tmp_path):
    from CITEgeist.model.citegeist_model import CitegeistModel

    names = ["CD20", "CD68", "CD163", "PanCK"]
    Y, dropped = P.dropout_channels(_X(), names, names, 0.5, np.random.default_rng(0))
    obs = pd.DataFrame(index=[f"s{i}" for i in range(Y.shape[0])])
    cite = ad.AnnData(X=Y.astype(float), obs=obs, var=pd.DataFrame(index=names))
    gex = ad.AnnData(
        X=np.random.default_rng(1).poisson(5, (Y.shape[0], 30)).astype(float),
        obs=obs.copy(),
        var=pd.DataFrame(index=[f"G{i}" for i in range(30)]),
    )
    m = CitegeistModel(
        sample_name="t",
        output_folder=str(tmp_path),
        simulation=True,
        gene_expression_adata=gex,
        antibody_capture_adata=cite,
    )
    m.preprocess_antibody()
    out = m.antibody_capture_adata
    assert np.isfinite(np.asarray(out.X)).all()
    raw = np.asarray(out.layers["raw_counts"])
    assert len(dropped) == 2
    for d in dropped:
        assert (raw[:, names.index(d)] == 0).all()


def test_background_noise_adds_nonnegative_counts():
    X = _X()
    bg = P.estimate_background(X)
    Y = P.add_background_noise(X, 2.0, bg, np.random.default_rng(0))
    assert (Y >= X).all() and bg.shape == (X.shape[1],) and (bg >= 0).all()


def test_estimate_background_bottom_half_mean():
    X = np.array([[0, 10], [2, 10], [4, 30], [100, 50]])
    # medians: col0 = 3.0 -> {0,2} mean 1.0 ; col1 = 20.0 -> {10,10} mean 10.0
    assert np.allclose(P.estimate_background(X), [1.0, 10.0])


def test_qp_active_markers_major_then_soft_ordered_unique():
    assert P.qp_active_markers(PROFILE) == ["CD20", "CD68", "CD163", "PanCK"]
    assert P.qp_active_markers(SOFT) == ["CD68", "CD20", "CD14"]


def test_active_profile_drops_minor_keeps_soft():
    ap = P.active_profile({**PROFILE, **SOFT})
    assert "Minor" not in str(ap) and ap["Myeloid"]["Soft"] == [("CD14", 0.5)]


@pytest.mark.parametrize("prof", [PROFILE, SHARED, SOFT])
def test_ablation_panels_balanced_nested_then_zero(prof):
    panels = P.ablation_panels(prof)
    names = [n for n, _ in panels]
    assert names == ["full", "one_per_type", "zero"]
    for name, d in panels[:-1]:  # every nonzero panel keeps every type with >=1 marker
        assert set(d) == set(prof)
        assert all(v.get("Major") or v.get("Soft") for v in d.values())
    assert panels[-1][1] == {}
    for (_, a), (_, b) in zip(panels, panels[1:]):
        assert _markers(b) <= _markers(a)


def test_one_per_type_prefers_type_specific_marker():
    one = dict(P.ablation_panels(SHARED))["one_per_type"]
    assert one["CD4+ T cells"]["Major"] == ["CD4"] and one["CD8+ T cells"]["Major"] == ["CD8A"]


def test_type_removal_series_cumulative_alphabetical():
    series = P.type_removal_series(PROFILE)
    assert [n for n, _ in series] == ["one_per_type", "minus_B cells", "minus_B cells+Epithelial", "zero"]
    assert [len(d) for _, d in series] == [3, 2, 1, 0]


def test_leave_one_type_out():
    loto = dict(P.leave_one_type_out(PROFILE))
    assert set(loto) == {"loto_B cells", "loto_Macrophages", "loto_Epithelial"}
    assert "Macrophages" not in loto["loto_Macrophages"]
