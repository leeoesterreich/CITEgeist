import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402

pytestmark = pytest.mark.unit
PROFILE = {"B cells": {"Major": ["CD20"]}, "Epithelial": {"Major": ["PanCK"]}}


def _cite():
    X = np.random.default_rng(0).poisson(10, size=(20, 3))
    return ad.AnnData(
        X=X.astype(float),
        var=pd.DataFrame(index=["CD20", "PanCK", "CD3E"]),
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(20)]),
    )


def test_missing_type_filled_zero():
    pred = pd.DataFrame({"B cells": [0.5, 1.0]}, index=["a", "b"])
    out = B.align_prediction(pred, ["B cells", "Epithelial"])
    assert list(out.columns) == ["B cells", "Epithelial"]
    assert (out["Epithelial"] == 0).all()
    assert np.allclose(out.sum(1), 1.0)


def test_apply_condition_baseline_is_identity():
    a = _cite()
    b, prof, meta = R.apply_condition(a.copy(), PROFILE, {"axis": "baseline"}, np.random.default_rng(0))
    assert np.array_equal(a.X, b.X) and prof == PROFILE and meta["status"] == "ok"


def test_apply_condition_ablation_drops_columns_and_dict():
    b, prof, meta = R.apply_condition(
        _cite(),
        PROFILE,
        {"axis": "ablation", "panel": "x", "profile": {"B cells": {"Major": ["CD20"]}}},
        np.random.default_rng(0),
    )
    assert "PanCK" not in b.var_names and "Epithelial" not in prof


def test_apply_condition_type_removal_behaves_like_ablation():
    b, prof, meta = R.apply_condition(
        _cite(),
        PROFILE,
        {"axis": "type_removal", "panel": "minus_Epithelial", "profile": {"B cells": {"Major": ["CD20"]}}},
        np.random.default_rng(0),
    )
    assert "PanCK" not in b.var_names and "Epithelial" not in prof
    assert meta["status"] == "ok" and meta["dropped"] == ["PanCK"]
    _, prof0, meta0 = R.apply_condition(
        _cite(), PROFILE, {"axis": "type_removal", "panel": "zero", "profile": {}}, np.random.default_rng(0)
    )
    assert meta0["status"] == "failed_no_markers" and prof0 == {}


def _fake_kw(qp_fn=None):
    gex = ad.AnnData(X=np.ones((20, 5)), obs=_cite().obs.copy(), var=pd.DataFrame(index=list("abcde")))
    props = pd.DataFrame({"B cells": 0.5, "Epithelial": 0.5}, index=_cite().obs_names)
    return dict(
        loader=lambda unit: (gex, _cite()),
        qp_fn=qp_fn or (lambda g, c, prof, unit: props),
        sace_fn=lambda g, c, prof, unit: {},
        prop_metrics_fn=lambda pred, unit: {"pearson_r": 0.9, "rmse": 0.1, "mae": 0.05, "jsd": 0.2, "per_type": {}},
        gex_metrics_fn=lambda layers, unit: {"gex_flat_r": 0.7, "gex_nrmse": 0.06},
        profile=PROFILE,
    )


def test_evaluate_unit_rows_have_exact_keys():
    rows = R.evaluate_unit("xenium", "0", 6, [{"axis": "baseline", "sentinel": "baseline", "sace": True}], **_fake_kw())
    assert len(rows) == 1 and list(rows[0]) == R.ROW_KEYS
    assert rows[0]["status"] == "ok" and rows[0]["gex_flat_r"] == 0.7


def test_zero_panel_records_failure_through_evaluate_unit():
    rows = R.evaluate_unit(
        "xenium",
        "0",
        6,
        [{"axis": "ablation", "panel": "zero", "profile": {}, "sentinel": "zero", "sace": True}],
        **_fake_kw(),
    )
    assert rows[0]["status"] == "failed_no_markers"
    assert all(np.isnan(rows[0][k]) for k in ["pearson_r", "rmse", "mae", "jsd", "gex_flat_r", "gex_nrmse"])


def test_solver_error_becomes_failed_row():
    def boom(*a):
        raise ValueError("No matching antibody markers found")

    rows = R.evaluate_unit(
        "xenium", "0", 6, [{"axis": "baseline", "sentinel": "", "sace": False}], **_fake_kw(qp_fn=boom)
    )
    assert rows[0]["status"].startswith("failed_") and np.isnan(rows[0]["pearson_r"])


def test_zero_panel_records_failure():
    b, prof, meta = R.apply_condition(
        _cite(), PROFILE, {"axis": "ablation", "panel": "zero", "profile": {}}, np.random.default_rng(0)
    )
    assert meta["status"] == "failed_no_markers" and prof == {}


def test_dropout_only_touches_qp_active():
    b, _, meta = R.apply_condition(_cite(), PROFILE, {"axis": "dropout", "level": 1.0}, np.random.default_rng(0))
    assert set(meta["dropped"]) == {"CD20", "PanCK"}
    assert (b[:, "CD3E"].X == _cite()[:, "CD3E"].X).all()


def test_sim_alignment_does_not_renormalize():
    pred = pd.DataFrame({"B cells": [0.2, 0.3]}, index=["a", "b"])
    out = B.align_prediction(pred, ["B cells", "Epithelial"], renormalize=False)
    assert list(out.columns) == ["B cells", "Epithelial"]
    assert np.allclose(out.sum(1), [0.2, 0.3]) and (out["Epithelial"] == 0).all()


def test_unexpected_exception_recorded_and_next_condition_runs():
    calls = []

    def flaky(g, c, prof, unit):
        calls.append(1)
        if len(calls) == 1:
            raise RuntimeError("cuOPT exploded. details")
        return _fake_kw()["qp_fn"](g, c, prof, unit)

    seen = []
    rows = R.evaluate_unit(
        "xenium",
        "0",
        6,
        [{"axis": "baseline", "sentinel": ""}, {"axis": "baseline", "sentinel": "baseline"}],
        on_row=seen.append,
        **_fake_kw(qp_fn=flaky),
    )
    assert rows[0]["status"] == "failed_RuntimeError_cuOPT_exploded" and np.isnan(rows[0]["pearson_r"])
    assert rows[1]["status"] == "ok" and len(seen) == 2


def test_save_pred_baseline_one_file_and_failed_none():
    saved = {}
    conds = [
        {"axis": "baseline", "sentinel": "", "sace": False},
        {"axis": "thinning", "level": 0.5, "seed": 1, "sentinel": "", "sace": False},
        {"axis": "ablation", "panel": "zero", "profile": {}, "sentinel": "", "sace": False},
    ]
    R.evaluate_unit(
        "xenium", "0", 6, conds, save_pred=lambda cid, pred: saved.setdefault(cid, pred), **_fake_kw()
    )
    assert sorted(saved) == ["baseline", "thinning_0.5_s1"]  # failed zero panel saves nothing
    assert saved["baseline"].shape == (20, 2)
    assert R.condition_id({"axis": "ablation", "panel": "one_per_type"}) == "ablation_one_per_type"


def test_load_conditions_accepts_long_inline_json_and_path(tmp_path):
    conds = [{"axis": "baseline", "pad": "x" * 6000}]
    import json as _j

    assert R.load_conditions(_j.dumps(conds)) == conds  # >4096 chars inline: no stat()
    f = tmp_path / "c.json"
    f.write_text(_j.dumps(conds))
    assert R.load_conditions(str(f)) == conds


def test_save_pred_io_error_does_not_wipe_metrics():
    def bad(cid, pred):
        raise OSError("disk full")

    rows = R.evaluate_unit(
        "xenium", "0", 6, [{"axis": "baseline", "sentinel": "", "sace": False}], save_pred=bad, **_fake_kw()
    )
    assert rows[0]["status"] == "ok" and rows[0]["pearson_r"] == 0.9


def _layers_fake_kw(**over):
    kw = _fake_kw(**over)
    kw["sace_fn"] = lambda g, c, prof, u: {
        "B-cells": pd.DataFrame(np.ones((3, 4)), index=list("abc"), columns=list("wxyz")),
        "CAFs": pd.DataFrame(np.zeros((3, 4)), index=list("abc"), columns=list("wxyz")),
    }
    return kw


def test_save_layers_writes_names_shape_and_props(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "_build_wiring", lambda d, n, v=None: _layers_fake_kw())
    out, lay = tmp_path / "out", tmp_path / "cur" / "mixed" / "Wu_rep_0"
    R.main(
        ["--dataset", "xenium", "--unit", "0", "--conditions",
         '[{"axis":"baseline","sentinel":"baseline","sace":true}]', "--out", str(out), "--save-layers", str(lay)]
    )
    assert sorted(p.name for p in (lay / "layers").iterdir()) == ["B-cells_layer.csv", "CAFs_layer.csv"]
    df = pd.read_csv(lay / "layers" / "B-cells_layer.csv", index_col=0)
    assert df.shape == (3, 4) and list(df.columns) == list("wxyz")
    assert pd.read_csv(lay / "props.csv", index_col=0).shape == (20, 2)


def test_no_layers_written_without_option(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "_build_wiring", lambda d, n, v=None: _layers_fake_kw())
    out = tmp_path / "out"
    R.main(["--dataset", "xenium", "--unit", "0", "--conditions",
            '[{"axis":"baseline","sentinel":"baseline","sace":true}]', "--out", str(out)])
    assert not list(tmp_path.rglob("*_layer.csv")) and not list(tmp_path.rglob("props.csv"))


def test_save_layers_skips_non_baseline_conditions():
    saved = []
    R.evaluate_unit(
        "xenium", "0", 6, [{"axis": "thinning", "level": 0.5, "seed": 1, "sentinel": "", "sace": False}],
        save_layers=lambda *a: saved.append(a), **_layers_fake_kw(),
    )
    assert saved == []


def test_save_layers_error_does_not_wipe_metrics_and_save_pred_still_runs():
    def bad(cid, layers, pred):
        raise OSError("disk full")

    preds = []
    rows = R.evaluate_unit(
        "xenium", "0", 6, [{"axis": "baseline", "sentinel": "baseline", "sace": True}],
        save_layers=bad, save_pred=lambda cid, pred: preds.append(cid), **_layers_fake_kw(),
    )
    assert rows[0]["status"] == "ok" and rows[0]["pearson_r"] == 0.9
    assert preds == ["baseline"]


def test_failure_message_stored_untruncated():
    long = "Cell proportion validation failed for sample 'Xenium_region_1': type 'Epithelial' mean 0.2% < 1%. more"

    def boom(*a):
        raise ValueError(long)

    rows = R.evaluate_unit(
        "xenium", "0", 6, [{"axis": "baseline", "sentinel": "", "sace": False}], **_fake_kw(qp_fn=boom)
    )
    assert len(rows[0]["status"]) <= len("failed_") + 60  # status slug unchanged (truncated)
    assert rows[0]["failure_message"] == "ValueError: " + long  # full message kept
    ok = R.evaluate_unit("xenium", "0", 6, [{"axis": "baseline", "sentinel": "", "sace": False}], **_fake_kw())
    assert ok[0]["failure_message"] == ""


def test_run_sim_does_not_mutate_callers_adata(monkeypatch):
    """run_sim must work on copies: evaluate_unit reuses the same per-condition cite AnnData for the
    QP call and then the SACE call, and CitegeistModel preprocessing overwrites .X (CLR) and
    layers["raw_counts"] in place. Without copies the 2nd call sees double-CLR'd antibody data."""
    from CITEgeist.model.citegeist_model import CitegeistModel

    rng = np.random.default_rng(0)
    n, g, m = 30, 40, 4
    obs = pd.DataFrame(index=[f"s{i}" for i in range(n)])
    gex = ad.AnnData(
        X=rng.poisson(5, size=(n, g)).astype(np.float32),
        obs=obs.copy(),
        var=pd.DataFrame(index=[f"g{i}" for i in range(g)]),
    )
    gex.layers["raw_counts"] = gex.X.copy()
    cite = ad.AnnData(
        X=rng.poisson(20, size=(n, m)).astype(np.float64),
        obs=obs.copy(),
        var=pd.DataFrame(index=[f"m{i}" for i in range(m)]),
    )
    cite.layers["raw_counts"] = cite.X.copy()
    gex_X0, cite_X0, cite_raw0 = gex.X.copy(), cite.X.copy(), cite.layers["raw_counts"].copy()

    # CPU-only: stub the cuOPT QP solve (the preprocessing is the code under test).
    fake = pd.DataFrame(1.0 / 2, index=gex.obs_names, columns=["A", "B"])
    monkeypatch.setattr(CitegeistModel, "run_cell_proportion_model", lambda self, **kw: (fake.copy(), None))

    B.run_sim(gex, cite, {"A": {"Major": ["m0"]}, "B": {"Major": ["m1"]}}, sace=False)

    assert np.array_equal(cite.X, cite_X0)
    assert np.array_equal(cite.layers["raw_counts"], cite_raw0)
    assert np.array_equal(gex.X, gex_X0)
