"""brightness_method plumbing in CitegeistModel.run_cell_proportion_model (truth table, raw-GEX source, diagnostics)."""
import logging

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import CITEgeist.model.deconvolution.brightness_estimators as BE
from CITEgeist.model.citegeist_model import CitegeistModel, _UNSET
from CITEgeist.model.deconvolution.marker_share import (
    estimate_marker_brightness_with_diagnostics,
    marker_share_proportions,
)

pytestmark = pytest.mark.unit

PROFILE = {
    "A": {"Major": ["m1", "m2"], "Minor": []},
    "B": {"Major": ["m3", "m4"], "Minor": []},
    "C": {"Major": ["m5"], "Minor": []},
}
MARKERS = ["m1", "m2", "m3", "m4", "m5"]
PANEL = {"A": ["g0", "g1", "g2"], "B": ["g3", "g4", "g5"], "C": ["g6", "g7", "g8"]}
N = 40
SPOTS = [f"s{i}" for i in range(N)]


def _build(tmp_path, raw_layer=True, gex_names=None):
    rng = np.random.default_rng(0)
    ab_X = rng.poisson(5, size=(N, len(MARKERS))).astype(float)
    ab = AnnData(X=ab_X, obs=pd.DataFrame(index=SPOTS), var=pd.DataFrame(index=MARKERS))
    gX = rng.poisson(6, size=(N, 12)).astype(float)
    gX[:, 8] = 1.0  # g8: mean over nonzero spots == 1.0, so filter_gex (> 1.0) drops it
    gex = AnnData(
        X=gX,
        obs=pd.DataFrame(index=gex_names if gex_names is not None else SPOTS),
        var=pd.DataFrame(index=[f"g{i}" for i in range(12)]),
    )
    if raw_layer:
        gex.layers["raw_counts"] = gex.X.copy()
    m = CitegeistModel(
        sample_name="t",
        output_folder=str(tmp_path),
        simulation=True,
        gene_expression_adata=gex,
        antibody_capture_adata=ab,
    )
    m.cell_profile_dict = PROFILE
    return m, ab_X


def _plain_expected(ab_X):
    df, _ = marker_share_proportions(ab_X, MARKERS, PROFILE, brightness=None)
    df.index = SPOTS
    return df.sort_index()


def _gmm_expected(ab_X):
    b, _ = estimate_marker_brightness_with_diagnostics(ab_X, MARKERS, MARKERS)
    df, _ = marker_share_proportions(ab_X, MARKERS, PROFILE, brightness=b)
    df.index = SPOTS
    return df.sort_index(), b


# (id, brightness_method, brightness_correction, expected)
ROWS = [
    (1, _UNSET, False, "plain"),
    (2, _UNSET, True, "gmm"),
    (3, None, False, "plain"),
    (4, None, True, "plain"),  # modern None (plain) takes precedence over legacy brightness_correction=True
    (5, "gmm", False, "gmm"),
    (6, "gmm", True, "gmm"),
    (7, "rna_slope", False, "rna"),
    (8, "rna_slope", True, ValueError),
    (9, "bogus", False, ValueError),
    (9.1, "bogus", True, ValueError),
]


@pytest.mark.parametrize("row,bm,bc,expected", ROWS, ids=[f"row{r[0]}" for r in ROWS])
def test_truth_table(tmp_path, row, bm, bc, expected):
    m, ab_X = _build(tmp_path)
    kw = dict(method="marker_share", brightness_correction=bc, rna_panel=PANEL)
    if bm is not _UNSET:
        kw["brightness_method"] = bm
    if expected is ValueError:
        with pytest.raises(ValueError):
            m.run_cell_proportion_model(**kw)
        return
    g, _ = m.run_cell_proportion_model(**kw)
    if expected == "plain":
        pd.testing.assert_frame_equal(g, _plain_expected(ab_X), check_exact=True)
        assert m.results["marker_brightness"] is None
    elif expected == "gmm":
        df, b = _gmm_expected(ab_X)
        pd.testing.assert_frame_equal(g, df, check_exact=True)
        assert m.results["marker_brightness"] == b
    else:
        assert set(m.results["marker_brightness"]) == set(MARKERS)
        assert np.allclose(g.sum(axis=1), 1.0)


def test_legacy_default_signature_values():
    import inspect

    p = inspect.signature(CitegeistModel.run_cell_proportion_model).parameters
    assert p["brightness_method"].default is _UNSET
    assert p["brightness_correction"].default is False
    assert p["allow_plain_fallback"].default is False and p["rna_panel"].default is None


def test_brightness_default_is_plain(tmp_path):
    """3.0.0 contract: no brightness args -> marker share with no brightness correction (the drop)."""
    m, ab_X = _build(tmp_path)
    g, _ = m.run_cell_proportion_model(method="marker_share")
    pd.testing.assert_frame_equal(g, _plain_expected(ab_X), check_exact=True)
    assert m.results["marker_brightness"] is None


def test_brightness_legacy_alias_warns_for_gmm(tmp_path):
    """Deprecated brightness_correction=True -> "gmm", with DeprecationWarning."""
    m, ab_X = _build(tmp_path)
    with pytest.warns(DeprecationWarning):
        g, _ = m.run_cell_proportion_model(method="marker_share", brightness_correction=True)
    df, _ = _gmm_expected(ab_X)
    pd.testing.assert_frame_equal(g, df, check_exact=True)


def test_brightness_modern_precedence_over_legacy(tmp_path):
    """Modern brightness_method=None (plain) wins over legacy brightness_correction=True (not gmm)."""
    m, ab_X = _build(tmp_path)
    with pytest.warns(DeprecationWarning):
        g, _ = m.run_cell_proportion_model(method="marker_share", brightness_method=None, brightness_correction=True)
    pd.testing.assert_frame_equal(g, _plain_expected(ab_X), check_exact=True)
    assert m.results["marker_brightness"] is None


def test_brightness_conflicting_nonnull_method_raises(tmp_path):
    """Legacy brightness_correction=True + an explicit non-"gmm" method is a genuine contradiction."""
    m, _ = _build(tmp_path)
    with pytest.warns(DeprecationWarning), pytest.raises(ValueError, match="contradicts brightness_correction"):
        m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", brightness_correction=True)


def test_rna_method_uses_raw_gex_not_filtered(tmp_path, monkeypatch):
    m, _ = _build(tmp_path, raw_layer=False)
    m.filter_gex(min_counts=1, keep_unfiltered=True)
    assert "g8" not in m.gene_expression_adata.var_names  # filtered out of the processed matrix
    seen = {}
    real = BE.estimate_brightness

    def spy(method, X, var_names, profile, gex=None, *a, **k):
        seen["gex"] = gex
        return real(method, X, var_names, profile, gex, *a, **k)

    monkeypatch.setattr(BE, "estimate_brightness", spy)
    m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)
    g = seen["gex"]
    assert "g8" in g.var_names and g.n_obs == N
    assert np.all(np.asarray(g.X)[:, list(g.var_names).index("g8")] == 1.0)  # raw counts, not normalized


def test_rna_method_obs_mismatch_raises(tmp_path):
    m, _ = _build(tmp_path, gex_names=[f"x{i}" for i in range(N)])
    with pytest.raises(ValueError, match="absent from raw GEX"):
        m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)


def test_rna_method_duplicate_obs_raises(tmp_path):
    m, _ = _build(tmp_path, gex_names=["s0"] + SPOTS[1:-1] + ["s0"])
    with pytest.raises(ValueError, match="duplicate"):
        m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)


def test_missing_raw_gex_raises_unless_allowed(tmp_path, caplog):
    m, ab_X = _build(tmp_path, raw_layer=False)
    with pytest.raises(ValueError, match="raw GEX"):
        m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)
    with caplog.at_level(logging.WARNING):
        g, _ = m.run_cell_proportion_model(
            method="marker_share", brightness_method="rna_slope", rna_panel=PANEL, allow_plain_fallback=True
        )
    assert "plain" in caplog.text
    pd.testing.assert_frame_equal(g, _plain_expected(ab_X))


def test_diagnostics_saved_in_results(tmp_path):
    m, _ = _build(tmp_path)
    m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)
    diag = m.results["marker_brightness_diagnostics"]
    assert set(diag) == set(MARKERS)
    assert all("used_brightness" in d for d in diag.values())


def test_default_filter_gex_holds_no_prefilter(tmp_path):
    m, _ = _build(tmp_path, raw_layer=False)
    m.filter_gex(min_counts=1)
    assert m._gex_prefilter is None


def test_filtered_without_keep_raises_with_hint_unless_allowed(tmp_path):
    m, ab_X = _build(tmp_path, raw_layer=False)
    m.filter_gex(min_counts=1)
    with pytest.raises(ValueError, match="keep_unfiltered=True"):
        m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)
    g, _ = m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL, allow_plain_fallback=True)
    pd.testing.assert_frame_equal(g, _plain_expected(ab_X))


def test_second_filter_gex_keeps_first_prefilter(tmp_path):
    m, _ = _build(tmp_path, raw_layer=False)
    m.filter_gex(min_counts=1, keep_unfiltered=True)
    first = m._gex_prefilter
    m.filter_gex(min_counts=1, keep_unfiltered=True)
    assert m._gex_prefilter is first and "g8" in first.var_names


def test_non_integer_x_without_layer_warns(tmp_path):
    m, _ = _build(tmp_path, raw_layer=False)
    m.filter_gex(min_counts=1, keep_unfiltered=True)
    m._gex_prefilter.X = m._gex_prefilter.X * 0.5 + 0.1
    with pytest.warns(UserWarning, match="not integer-valued"):
        m.run_cell_proportion_model(method="marker_share", brightness_method="rna_slope", rna_panel=PANEL)


def test_init_sets_prefilter_none(tmp_path):
    m, _ = _build(tmp_path)
    assert m._gex_prefilter is None
