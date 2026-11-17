"""Stress-test perturbations must reach the marker-share path (otherwise robustness would look flat for no reason).

Uses the REAL CitegeistModel (CPU, marker_share needs no GPU) on small synthetic data; the perturbation is
applied with run_sensitivity.apply_condition exactly as the sweep does. Imports private Benchmarking code.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402

PROFILE = {
    "A": {"Major": ["m0", "m1"], "Minor": ["m4"]},
    "B": {"Major": ["m2"], "Minor": ["m5"]},
    "C": {"Major": ["m3"], "Minor": []},
}


@pytest.fixture
def data(monkeypatch):
    monkeypatch.setattr(B, "load_nuclei_counts", lambda region: pd.Series(dtype=int))
    rng = np.random.default_rng(0)
    n = 60
    spots = [f"spot_{i}" for i in range(n)]
    gex = ad.AnnData(X=rng.poisson(5, (n, 40)).astype(float), obs=pd.DataFrame(index=spots))
    gex.var_names = [f"g{i}" for i in range(40)]
    cite = ad.AnnData(X=rng.poisson(rng.uniform(2, 30, 6), (n, 6)).astype(float), obs=pd.DataFrame(index=spots))
    cite.var_names = [f"m{i}" for i in range(6)]
    return gex, cite


def _props(gex, cite, cond, profile=PROFILE):
    rng = np.random.default_rng(cond.get("seed", 0))
    c2, prof, meta = R.apply_condition(cite.copy(), profile, cond, rng)
    assert meta["status"] == "ok"
    extra = cond.get("setting")
    return B.run_xenium_qp(gex.copy(), c2, prof, 0, setting=extra, phase1="marker_share")


@pytest.mark.unit
@pytest.mark.parametrize(
    "cond",
    [
        {"axis": "thinning", "level": 0.1, "seed": 0},
        {"axis": "dropout", "level": 0.5, "seed": 0},
        {"axis": "noise", "level": 4.0, "seed": 0},
        {"axis": "ablation", "panel": "one_per_type", "profile": {"A": {"Major": ["m0"]}, "B": {"Major": ["m2"]}, "C": {"Major": ["m3"]}}},
        {"axis": "type_removal", "panel": "minus_C", "profile": {"A": {"Major": ["m0"]}, "B": {"Major": ["m2"]}}},
    ],
)
def test_perturbation_changes_marker_share_proportions(data, cond):
    gex, cite = data
    base = _props(gex, cite, {"axis": "baseline"})
    pert = _props(gex, cite, cond)
    cols = [c for c in pert.columns if c in base.columns]
    assert not np.allclose(base[cols].values, pert[cols].values), f"{cond['axis']} did not change marker-share output"


@pytest.mark.unit
def test_thinning_level_one_is_identity(data):
    gex, cite = data
    base = _props(gex, cite, {"axis": "baseline"})
    same = _props(gex, cite, {"axis": "thinning", "level": 1.0, "seed": 0})
    np.testing.assert_allclose(base.values, same.values)


@pytest.mark.unit
def test_markers_setting_reaches_marker_share_and_qp_knobs_are_ignored_loudly(data, caplog):
    gex, cite = data
    base = _props(gex, cite, {"axis": "baseline"})
    mk = _props(gex, cite, {"axis": "setting", "param": "markers", "level": "major_plus_minor", "setting": {"axis": "setting", "param": "markers", "level": "major_plus_minor"}})
    assert not np.allclose(base.values, mk.values)
    B._QP_ONLY_LOGGED.clear()
    knob = {"axis": "setting", "param": "lambda_reg", "level": 2}
    with caplog.at_level("WARNING"):
        k1 = _props(gex, cite, {"axis": "setting", "param": "lambda_reg", "level": 2, "setting": knob})
        _props(gex, cite, {"axis": "setting", "param": "lambda_reg", "level": 2, "setting": knob})
    np.testing.assert_allclose(base.values, k1.values)
    assert sum("lambda_reg" in r.message for r in caplog.records) == 1  # logged once, not per call
