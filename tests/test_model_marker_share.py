"""Tests for the marker-share default in CitegeistModel.run_cell_proportion_model."""
import inspect
import logging
import os
import sys
import warnings

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from CITEgeist.model.citegeist_model import CitegeistModel

pytestmark = pytest.mark.unit

PROFILE = {
    "A": {"Major": ["m1", "m2"], "Minor": []},
    "B": {"Major": ["m3", "m4"], "Minor": []},
    "C": {"Major": ["m5"], "Minor": []},
}
MARKERS = {m for v in PROFILE.values() for m in v["Major"]}


@pytest.fixture
def model(tmp_path):
    rng = np.random.default_rng(0)
    n = 30
    names = ["m1", "m2", "m3", "m4", "m5", "m6"]
    X = rng.poisson(5, size=(n, len(names))).astype(float)
    X[0, :5] = 0  # an all-zero-marker spot
    ab = AnnData(X=X, obs=pd.DataFrame(index=[f"s{i}" for i in range(n)]), var=pd.DataFrame(index=names))
    g = int(np.ceil(np.sqrt(n)))
    ab.obsm["spatial"] = np.array([[i % g, i // g] for i in range(n)], dtype=float) * 10
    gex = AnnData(X=rng.poisson(2, size=(n, 8)).astype(float), obs=ab.obs.copy(),
                  var=pd.DataFrame(index=[f"g{i}" for i in range(8)]))
    gex.obsm["spatial"] = ab.obsm["spatial"].copy()
    m = CitegeistModel(sample_name="t", output_folder=str(tmp_path), simulation=True,
                       gene_expression_adata=gex, antibody_capture_adata=ab)
    m.cell_profile_dict = PROFILE
    return m


def test_default_method_is_auto():
    sig = inspect.signature(CitegeistModel.run_cell_proportion_model)
    assert sig.parameters["method"].default == "auto"
    assert sig.parameters["brightness_correction"].default is False


def test_marker_share_runs_without_cuopt(model, monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "cuopt", None)
    g, f = model.run_cell_proportion_model(method="marker_share")
    prop = model.results["cell_prop"]
    assert np.allclose(prop.sum(axis=1), 1.0)
    assert list(prop.index) == sorted(model.antibody_capture_adata.obs_names)
    assert model.results["marker_brightness"] is None
    assert model.results["n_all_zero_marker_spots"] == 1
    assert "cell_prop_base" in model.results
    assert f.equals(g)
    for suffix in ("global", "finetuned"):
        assert os.path.exists(tmp_path / f"t_cell_prop_{suffix}_results.csv")


def test_brightness_flag_stores_dict(model):
    model.run_cell_proportion_model(method="marker_share", brightness_correction=True)
    assert set(model.results["marker_brightness"]) == MARKERS


def test_unknown_method_rejected(model):
    with pytest.raises(ValueError):
        model.run_cell_proportion_model(method="bogus")


def test_marker_share_needs_no_coords_or_radius(model, caplog):
    for ad in (model.antibody_capture_adata, model.gene_expression_adata):
        del ad.obsm["spatial"]
    with caplog.at_level(logging.INFO):
        model.run_cell_proportion_model(method="marker_share")
    assert np.allclose(model.results["cell_prop"].sum(axis=1), 1.0)
    assert "Auto-detected radius" not in caplog.text
    assert "No spatial coordinates" not in caplog.text


def test_qp_only_arg_warning(model, caplog):
    with caplog.at_level(logging.WARNING):
        model.run_cell_proportion_model(method="marker_share", lambda_laplacian=0.5)
    msgs = [r.getMessage() for r in caplog.records if "QP-only" in r.getMessage()]
    assert len(msgs) == 1
    assert "lambda_laplacian" in msgs[0]


def test_no_qp_only_warning_for_defaults(model, caplog):
    with caplog.at_level(logging.WARNING):
        model.run_cell_proportion_model(method="marker_share")
    assert "QP-only" not in caplog.text


def test_detection_gating_false_silent_for_marker_share(model):
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        model.run_cell_proportion_model(method="marker_share", use_detection_gating=False)
    assert not [x for x in w if "degrades proportion accuracy" in str(x.message)]


def test_model_marker_share_uses_raw_counts_after_preprocess(model):
    from CITEgeist.model.deconvolution.marker_share import marker_share_proportions

    ab = model.antibody_capture_adata
    raw = np.asarray(ab.X).copy()
    var_names = list(ab.var_names)
    expected = marker_share_proportions(raw, var_names, PROFILE)[0]
    expected.index = ab.obs_names
    expected = expected.sort_index()

    model.preprocess_antibody()
    ab = model.antibody_capture_adata
    clr = np.asarray(ab.X.toarray() if hasattr(ab.X, "toarray") else ab.X)
    wrong = marker_share_proportions(clr, list(ab.var_names), PROFILE)[0]
    assert not np.allclose(expected.values, wrong.values, atol=1e-6)  # test can fail

    model.run_cell_proportion_model(method="marker_share")
    got = model.results["cell_prop"]
    np.testing.assert_allclose(got.values, expected.loc[got.index].values, atol=1e-12)


def test_brightness_diagnostics_stored_with_flag(model):
    model.run_cell_proportion_model(method="marker_share", brightness_correction=True)
    diag = model.results["marker_brightness_diagnostics"]
    assert set(diag) == MARKERS
    for v in diag.values():
        base = {"used_brightness", "fallback", "reason"}
        assert set(v) == (base | {"fallback_value"} if v["fallback"] else base)


def test_brightness_diagnostics_none_without_flag(model):
    model.run_cell_proportion_model(method="marker_share")
    assert model.results["marker_brightness_diagnostics"] is None
