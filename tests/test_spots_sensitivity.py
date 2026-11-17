"""Runner glue for SPOTS-under-perturbation: shared perturbation/seeding, zero-marker handling, matched-unit aggregation."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import perturb as P  # noqa: E402
from protein_sensitivity import spots_sensitivity as S  # noqa: E402

pytestmark = pytest.mark.unit

PROFILE = {"A": {"Major": ["A1", "A2"]}, "B": {"Major": ["B1"]}, "C": {"Major": ["C1"]}}
VARS = ["A1", "A2", "B1", "C1", "X1"]


@pytest.fixture
def cite(monkeypatch):
    monkeypatch.setitem(B.PROFILE, 6, PROFILE)
    rng = np.random.default_rng(42)
    X = rng.integers(1, 60, size=(12, 5)).astype(float)
    return ad.AnnData(X=X, obs=pd.DataFrame(index=[f"s{i}" for i in range(12)]), var=pd.DataFrame(index=VARS))


def metrics_capture(store):
    def fn(pred, u):
        store.append(pred.copy())
        return {"pearson_r": 0.5, "rmse": 0.1, "mae": 0.1, "jsd": 0.1, "per_type": {}}
    return fn


def test_same_perturbed_input_as_citegeist_seed(cite):
    store = []
    cond = [{"axis": "thinning", "level": 0.5, "seed": 1, "sentinel": ""}]
    S.run_unit(6, "0", cond, cite=cite.copy(), metrics_fn=metrics_capture(store))
    Xp = P.thin_counts(cite.X, 0.5, np.random.default_rng(1)).astype(float)  # the CITEgeist rng recipe
    expect, _ = S.spots_predict(Xp, VARS, PROFILE, list(PROFILE))
    assert np.allclose(store[0].to_numpy(), expect.to_numpy())
    other, _ = S.spots_predict(P.thin_counts(cite.X, 0.5, np.random.default_rng(2)).astype(float), VARS, PROFILE, list(PROFILE))
    assert not np.allclose(store[0].to_numpy(), other.to_numpy())  # seed matters


def test_dropout_zeroes_same_channels(cite):
    store = []
    cond = [{"axis": "dropout", "level": 0.5, "seed": 0, "sentinel": ""}]
    S.run_unit(6, "0", cond, cite=cite.copy(), metrics_fn=metrics_capture(store))
    Y, dropped = P.dropout_channels(cite.X, VARS, P.qp_active_markers(PROFILE), 0.5, np.random.default_rng(0))
    expect, _ = S.spots_predict(Y.astype(float), VARS, PROFILE, list(PROFILE))
    assert dropped and np.allclose(store[0].to_numpy(), expect.to_numpy())


def test_zero_markers_uniform_and_ok(cite):
    store = []
    cond = [{"axis": "ablation", "panel": "zero", "profile": {}, "sentinel": ""}]
    rows = S.run_unit(6, "0", cond, cite=cite.copy(), metrics_fn=metrics_capture(store))
    assert rows[0]["status"] == "ok"
    assert np.allclose(store[0].to_numpy(), 1 / 3)


def test_zero_markers_metric_failure_stays_visible(cite):
    def boom(pred, u):
        raise ValueError("constant input")
    cond = [{"axis": "ablation", "panel": "zero", "profile": {}, "sentinel": ""}]
    r = S.run_unit(6, "0", cond, cite=cite.copy(), metrics_fn=boom)[0]
    assert r["status"] == "failed_no_markers" and np.isnan(r["pearson_r"]) and "ValueError" in r["failure_message"]


def test_leave_one_type_out_zero_column_and_channel_removed(cite):
    store = []
    prof = P.leave_one_type_out(PROFILE)[1][1]  # drop B
    cond = [{"axis": "ablation", "panel": "loto_B", "profile": prof, "sentinel": ""}]
    S.run_unit(6, "0", cond, cite=cite.copy(), metrics_fn=metrics_capture(store))
    assert (store[0]["B"] == 0).all()
    assert np.allclose(store[0][["A", "C"]].sum(axis=1), 1.0)


def _frame(unit, r, status="ok", axis="dropout", level=0.5):
    return {"dataset": "xenium", "unit": unit, "n_types": 6, "axis": axis, "level": level, "seed": 0, "panel": "",
            "status": status, "pearson_r": r}


def test_matched_unit_aggregation_and_summary_assert():
    sp = pd.DataFrame([_frame(str(u), r) for u, r in enumerate([0.1, 0.2, 0.9, 0.4, 0.5])])
    cg = pd.DataFrame([_frame("0", 0.9), _frame("1", 0.8), _frame("2", np.nan, "failed_x"), _frame("3", 0.6), _frame("4", 0.7)])
    cg_mean = np.mean([0.9, 0.8, 0.6, 0.7])
    summ = {"xenium_6": {"dropout": {"0.5": {"pearson_r_mean": cg_mean, "status": "partial_failed"}}}}
    row = S.build_comparison(sp, cg, summ, "src.json", 6)[0]
    assert row["n_units"] == 5 and row["spots_r_mean"] == pytest.approx(0.42)
    assert row["spots_r_mean_on_citegeist_ok_units"] == pytest.approx(np.mean([0.1, 0.2, 0.4, 0.5]))
    assert row["citegeist_failed_units"] == ["2"] and row["n_citegeist_ok_units"] == 4
    assert row["delta_all_units"] == pytest.approx(cg_mean - 0.42)
    assert row["delta_matched_units"] == pytest.approx(cg_mean - 0.3)
    bad = {"xenium_6": {"dropout": {"0.5": {"pearson_r_mean": 0.123}}}}
    with pytest.raises(AssertionError):
        S.build_comparison(sp, cg, bad, "src.json", 6)
