"""Tests for the v40 settings-sensitivity axis (W3) and leave-one-region-out selection (W3b)."""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
from protein_sensitivity import aggregate as A  # noqa: E402
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import loro as L  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402
from protein_sensitivity import settings_conditions as SC  # noqa: E402

pytestmark = pytest.mark.unit

# Independent literal copy of the spec's one-at-a-time grid (NOT read from the protocol).
SPEC_LEVELS = {
    "detection_gate_ub": [0.005, 0.01, 0.02, 0.05],
    "lambda_laplacian": [0.0, 0.01, 0.03, 0.1],
    "lambda_reg": [0.5, 1, 2],
    "use_detection_gating": ["on", "off"],
    "use_gex_detection": ["on", "off"],
    "markers": ["major", "major_plus_minor"],
    "nuclei_source": ["label_derived", "all_segmented", "stardist_dapi", "none"],
}
SPEC_BASE = {
    "detection_gate_ub": 0.01,
    "lambda_laplacian": 0.0,
    "lambda_reg": 1,
    "use_detection_gating": "on",
    "use_gex_detection": "on",
    "markers": "major",
    "nuclei_source": "label_derived",
}
PROFILE = {
    "B cells": {"Major": ["CD20"], "Minor": ["CD45RA"]},
    "T cells": {"Major": ["CD3E"], "Minor": ["CD45RO", "CD3E"]},
    "Epithelial": {"Major": ["PanCK"], "Minor": []},
}
NUC = pd.Series({"spot_1": 3, "spot_2": 5})


# ---------------- condition set ----------------
def _independent_nonbaseline():
    return {(p, str(l)) for p, ls in SPEC_LEVELS.items() for l in ls if l != SPEC_BASE[p]}


def test_condition_set_matches_independent_derivation_and_literal_count():
    conds = SC.build_list(6)
    assert len(conds) == 15  # 14 non-baseline + 1 baseline (literal)
    assert sum(c["axis"] == "baseline" for c in conds) == 1
    got = {(c["param"], str(c["level"])) for c in conds if c["axis"] == "setting"}
    assert got == _independent_nonbaseline()
    assert len(got) == 14
    assert json.loads(SC.build(6)) == conds


def test_protocol_levels_and_baseline_equal_spec_literals():
    p = SC.load_protocol()
    assert {k: v["levels"] for k, v in p["settings"].items()} == SPEC_LEVELS
    assert p["baseline"]["levels"] == SPEC_BASE
    assert p["robustness_rule"]["mean_abs_threshold"] == 0.02
    assert p["robustness_rule"]["region_abs_threshold"] == 0.05
    assert p["expected_conditions"]["total_per_region_per_ntype"] == 15


def test_protocol_covers_five_regions_with_stardist_sha():
    p = SC.load_protocol()
    assert p["regions"] == [0, 1, 2, 3, 4]
    assert p["stardist_dapi"]["regions_covered"] == [0, 1, 2, 3, 4]
    assert len([k for k in p["stardist_dapi"]["sha256"] if k.startswith("stardist_counts_")]) == 5


# ---------------- profile / nuclei helpers ----------------
def test_major_plus_minor_promotes_minor_into_major():
    out = SC.promote_minor(PROFILE)
    assert out["B cells"]["Major"] == ["CD20", "CD45RA"]
    assert out["T cells"]["Major"] == ["CD3E", "CD45RO"]  # no duplicate CD3E
    assert out["Epithelial"]["Major"] == ["PanCK"]
    assert PROFILE["B cells"]["Major"] == ["CD20"]  # input untouched


def _mapping():
    idx = ["c1", "c2", "c3", "c4", "c5", "c6"]
    return (
        pd.DataFrame(
            {
                "spot_idx": [1, 1, 1, 2, -1, 1],
                "spot_id": ["spot_1", "spot_1", "spot_1", "spot_2", None, "spot_9"],
                "region_id": [0, 0, 0, 0, 0, 1],
            },
            index=idx,
        ),
        pd.DataFrame({"cell_type": ["T", "Unknown", "B", "Unknown", "T", "T"]}, index=idx),
    )


def test_all_segmented_counts_unknown_label_derived_does_not():
    m, ct = _mapping()
    allc = SC.all_segmented_counts(m, 0)
    lab = SC.label_derived_counts(m, ct, 0)
    assert allc.to_dict() == {"spot_1": 3, "spot_2": 1}
    assert lab.to_dict() == {"spot_1": 2}
    assert "spot_9" not in allc.index  # other region


def test_stardist_ids_map_to_spot_ids(tmp_path):
    (tmp_path / "stardist_counts_region_0.csv").write_text(",nuclei_count_raw\nspot_1,4\nspot_2,7\n")
    s = SC.stardist_counts(0, tmp_path)
    m, _ = _mapping()
    assert s.to_dict() == {"spot_1": 4, "spot_2": 7}
    assert set(s.index) == set(SC.all_segmented_counts(m, 0).index)
    with pytest.raises(FileNotFoundError):
        SC.stardist_counts(3, tmp_path)


@pytest.mark.skipif(not SC.STARDIST_DIR.exists(), reason="untracked StarDist counts absent (CI)")
def test_real_stardist_equals_oracle_index_and_all_segmented_equals_oracle():
    m = pd.read_csv(SC.XENIUM_DATA / "cell_to_spot_mapping.csv", index_col=0)
    for r in range(5):
        sd = SC.stardist_counts(r)
        oracle = pd.read_csv(SC.STARDIST_DIR / f"oracle_counts_region_{r}.csv", index_col=0)["n_cells"]
        allc = SC.all_segmented_counts(m, r)
        assert set(sd.index) == set(allc.index)
        assert (allc.reindex(oracle.index) == oracle).all() and len(allc) == len(oracle)


# ---------------- solver wiring (mock the model) ----------------
class FakeModel:
    last = None

    def __init__(self, **kw):
        FakeModel.last = self
        self.loaded = None
        self.run_kwargs = None
        self.results = {}

    def filter_gex(self, **kw):
        pass

    def preprocess_gex(self, **kw):
        pass

    def preprocess_antibody(self, **kw):
        pass

    def load_cell_profile_dict(self, p):
        self.loaded = p

    def run_cell_proportion_model(self, **kw):
        self.run_kwargs = kw
        self.results["cell_prop"] = pd.DataFrame({t: [0.5, 0.5] for t in PROFILE}, index=["spot_1", "spot_2"])


@pytest.fixture
def fake_model(monkeypatch):
    import CITEgeist.model.citegeist_model as cm

    monkeypatch.setattr(cm, "CitegeistModel", FakeModel)
    monkeypatch.setattr(B, "load_nuclei_counts", lambda region: NUC)
    FakeModel.last = None


def _adatas():
    g = ad.AnnData(X=np.ones((2, 2)), obs=pd.DataFrame(index=["spot_1", "spot_2"]))
    return g, g.copy()


def test_baseline_condition_kwargs_equal_current_published_call(fake_model):
    g, c = _adatas()
    B.run_xenium_qp(g, c, PROFILE, 0)
    published = FakeModel.last.run_kwargs
    assert set(published) == {"method", "lambda_laplacian", "nuclei_counts"} and published["method"] == "qp" and published["lambda_laplacian"] == 0.0
    assert published["nuclei_counts"] is NUC
    FakeModel.last = None
    B.run_xenium_qp(g, c, PROFILE, 0, setting={"axis": "baseline"})
    assert FakeModel.last.run_kwargs.keys() == published.keys()
    assert FakeModel.last.run_kwargs["nuclei_counts"] is NUC and FakeModel.last.loaded == PROFILE


@pytest.mark.parametrize("param,level", sorted(_independent_nonbaseline()))
def test_each_override_reaches_solver_call(fake_model, monkeypatch, param, level):
    lvl = next(l for l in SPEC_LEVELS[param] if str(l) == level)
    monkeypatch.setattr(SC, "all_segmented_counts", lambda m, r: pd.Series({"spot_1": 9}))
    monkeypatch.setattr(SC, "stardist_counts", lambda r, directory=SC.STARDIST_DIR: pd.Series({"spot_1": 11}))
    monkeypatch.setattr(pd, "read_csv", lambda *a, **k: pd.DataFrame())
    g, c = _adatas()
    B.run_xenium_qp(g, c, PROFILE, 0, setting={"axis": "setting", "param": param, "level": lvl})
    kw, model = FakeModel.last.run_kwargs, FakeModel.last
    if param in ("detection_gate_ub", "lambda_laplacian", "lambda_reg"):
        assert kw[param] == lvl
    elif param in ("use_detection_gating", "use_gex_detection"):
        assert kw[param] is False
    elif param == "markers":
        assert model.loaded["B cells"]["Major"] == ["CD20", "CD45RA"]
    else:
        expect = {"all_segmented": 9, "stardist_dapi": 11, "none": None}[lvl]
        got = kw["nuclei_counts"]
        assert (got is None) if expect is None else got.iloc[0] == expect
    # nothing else drifts from the published call: only the overridden kwarg is added/changed
    expected_keys = {"method", "lambda_laplacian", "nuclei_counts"} | ({param} if param not in ("markers", "nuclei_source") else set())
    assert set(kw) == expected_keys
    if param not in ("lambda_laplacian", "nuclei_source"):
        assert kw["lambda_laplacian"] == 0.0


# ---------------- runner support ----------------
def test_condition_id_for_setting_includes_param_and_handles_strings():
    assert R.condition_id({"axis": "setting", "param": "detection_gate_ub", "level": 0.01}) == "setting_detection_gate_ub_0.01"
    assert R.condition_id({"axis": "setting", "param": "nuclei_source", "level": "stardist_dapi"}) == "setting_nuclei_source_stardist_dapi"
    assert R.condition_id({"axis": "setting", "param": "lambda_reg", "level": 2}) == "setting_lambda_reg_2"
    assert R.condition_id({"axis": "thinning", "level": 0.5, "seed": 1}) == "thinning_0.5_s1"  # unchanged
    assert R.condition_id({"axis": "baseline"}) == "baseline"


def test_all_setting_condition_ids_unique():
    ids = [R.condition_id(c) for c in SC.build_list(6)]
    assert len(set(ids)) == len(ids) == 15


def test_apply_condition_setting_changes_nothing():
    X = np.random.default_rng(0).poisson(10, size=(5, 2)).astype(float)
    a = ad.AnnData(X=X, var=pd.DataFrame(index=["CD20", "PanCK"]))
    b, prof, meta = R.apply_condition(
        a.copy(), PROFILE, {"axis": "setting", "param": "lambda_reg", "level": 2}, np.random.default_rng(0)
    )
    assert np.array_equal(a.X, b.X) and prof == PROFILE and meta["status"] == "ok"


def test_row_keys_gain_param_and_existing_keys_unchanged():
    assert "param" in R.ROW_KEYS
    assert [k for k in R.ROW_KEYS if k != "param"] == [
        "dataset", "unit", "n_types", "axis", "level", "seed", "panel", "sentinel", "status",
        "pearson_r", "rmse", "mae", "jsd", "gex_flat_r", "gex_nrmse", "per_type_json", "failure_message",
    ]  # fmt: skip


def _evaluate(conds, calls):
    def qp_fn(g, c, prof, unit, *extra):
        calls.append(extra)
        return pd.DataFrame({"a": [1.0]})

    return R.evaluate_unit(
        "xenium", "0", 6, conds,
        loader=lambda u: (ad.AnnData(X=np.ones((2, 1))), ad.AnnData(X=np.ones((2, 1)), var=pd.DataFrame(index=["m"]))),
        qp_fn=qp_fn, sace_fn=None,
        prop_metrics_fn=lambda p, u: {"pearson_r": 0.5, "rmse": 1.0, "mae": 1.0, "jsd": 1.0},
        gex_metrics_fn=None, profile=PROFILE,
    )  # fmt: skip


def test_evaluate_unit_passes_condition_to_qp_fn_only_for_setting_axis():
    calls = []
    set_cond = {"axis": "setting", "param": "lambda_reg", "level": 2}
    rows = _evaluate([{"axis": "baseline"}, set_cond], calls)
    assert calls[0] == () and calls[1] == (set_cond,)
    assert rows[1]["param"] == "lambda_reg" and rows[1]["level"] == 2 and rows[0]["param"] == ""
    assert rows[1]["status"] == "ok"


def test_run_meta_records_settings_protocol_md5(tmp_path, monkeypatch):
    monkeypatch.setattr(R, "_build_wiring", lambda d, n, v=None: dict(
        loader=lambda u: (None, ad.AnnData(X=np.ones((2, 1)), var=pd.DataFrame(index=["m"]))),
        qp_fn=lambda *a: pd.DataFrame({"a": [1.0]}), sace_fn=None,
        prop_metrics_fn=lambda p, u: {"pearson_r": 0.5, "rmse": 1, "mae": 1, "jsd": 1},
        gex_metrics_fn=None, profile=PROFILE))  # fmt: skip
    R.main(["--dataset", "xenium", "--unit", "0", "--conditions", '[{"axis":"baseline"}]', "--out", str(tmp_path)])
    meta = json.loads((tmp_path / "run_meta.json").read_text())
    assert meta["settings_protocol_md5"] == SC.protocol_md5()


# ---------------- aggregation ----------------
def _rows(r_by):
    """r_by: {(axis,param,level): [r per region 0..4]} -> metrics frame."""
    rows = []
    for (axis, param, level), rs in r_by.items():
        for reg, r in enumerate(rs):
            rows.append(dict(dataset="xenium", unit=str(reg), n_types=6, axis=axis, param=param, level=level,
                             status="ok" if not pd.isna(r) else "failed_x", pearson_r=r))  # fmt: skip
    return pd.DataFrame(rows)


BASE_R = [0.70, 0.72, 0.74, 0.76, 0.78]


def test_summarize_settings_groups_by_param_level_and_applies_rule():
    df = _rows({
        ("baseline", "", float("nan")): BASE_R,
        ("setting", "lambda_reg", 2): [r - 0.01 for r in BASE_R],             # robust
        ("setting", "lambda_reg", 0.5): [r - 0.03 for r in BASE_R],           # mean breaks 0.02
        ("setting", "use_detection_gating", "off"): [r - 0.01 for r in BASE_R[:4]] + [BASE_R[4] - 0.06],  # one region >0.05
        ("setting", "nuclei_source", "none"): [r - 0.01 for r in BASE_R[:4]] + [float("nan")],           # failed cell
    })  # fmt: skip
    s = A.summarize_settings(df, 6, SC.load_protocol())
    by = {(x["param"], str(x["level"])): x for x in s["levels"]}
    assert by[("lambda_reg", "2")]["robust"] is True
    assert by[("lambda_reg", "0.5")]["robust"] is False
    assert by[("use_detection_gating", "off")]["robust"] is False
    assert by[("use_detection_gating", "off")]["worst_region_delta"] == pytest.approx(-0.06)
    assert by[("lambda_reg", "2")]["mean_delta_r"] == pytest.approx(-0.01)
    failed = by[("nuclei_source", "none")]
    assert failed["robust"] is False and failed["n_failed"] == 1 and "failed" in failed["status"]
    assert len(s["levels"]) == 14  # every prespecified level is reported, run or not
    absent = by[("lambda_laplacian", "0.1")]
    assert absent["robust"] is False and absent["n_missing"] == 5 and absent["mean_delta_r"] is None


def test_baseline_gate_check_blocks_deltas_when_baseline_drifts():
    df = _rows({("baseline", "", float("nan")): BASE_R})
    ok = A.check_baseline_gate(df, 6, {str(i): r for i, r in enumerate(BASE_R)}, 1e-4)
    assert ok["pass"] is True
    bad = A.check_baseline_gate(df, 6, {str(i): r + 0.001 for i, r in enumerate(BASE_R)}, 1e-4)
    assert bad["pass"] is False and len(bad["mismatches"]) == 5


# ---------------- LORO ----------------
def _grid(r):
    """r: {(param, level): [r per region 0..4]}"""
    return pd.DataFrame(
        [(6, reg, p, str(l), v) for (p, l), vs in r.items() for reg, v in enumerate(vs)],
        columns=["n_types", "region", "param", "level", "pearson_r"],
    )


def test_loro_select_picks_best_on_other_regions_and_excludes_heldout():
    # level "2" is best on regions 1-4 but terrible on region 0.
    grid = _grid({
        ("lambda_reg", 1): [0.70, 0.70, 0.70, 0.70, 0.70],
        ("lambda_reg", 2): [0.10, 0.80, 0.80, 0.80, 0.80],
        ("lambda_reg", 0.5): [0.90, 0.60, 0.60, 0.60, 0.60],
    })  # fmt: skip
    sel = L.select(grid, {"lambda_reg": 1}).set_index("held_out")
    # holding out region 0: others favour level 2 (0.80) -> chosen "2"; held-out r there is 0.10
    assert sel.loc[0, "chosen_level"] == "2"
    assert sel.loc[0, "heldout_r_chosen"] == pytest.approx(0.10)
    assert sel.loc[0, "heldout_r_baseline"] == pytest.approx(0.70)
    assert sel.loc[0, "delta"] == pytest.approx(-0.60)
    # holding out region 1: others (0,2,3,4) means: lvl1 .70, lvl2 (.1+.8*3)/4=.625, lvl0.5 (.9+.6*3)/4=.675 -> baseline
    assert sel.loc[1, "chosen_level"] == "1" and sel.loc[1, "delta"] == pytest.approx(0.0)
    # the held-out region's own r must not influence training mean
    assert sel.loc[0, "train_mean_r_chosen"] == pytest.approx(0.80)


def test_loro_tie_goes_to_baseline_and_failed_candidate_ineligible():
    grid = _grid({
        ("lambda_reg", 1): [0.7] * 5,
        ("lambda_reg", 2): [0.7] * 5,                      # exact tie
        ("lambda_reg", 0.5): [0.9, 0.9, 0.9, 0.9, np.nan],  # NaN in one region
    })  # fmt: skip
    sel = L.select(grid, {"lambda_reg": 1}).set_index("held_out")
    assert (sel.loc[[0, 1, 2, 3], "chosen_level"] == "1").all()  # 0.5 ineligible (NaN in region 4); 1 vs 2 tie -> baseline
    assert sel.loc[4, "chosen_level"] == "0.5"  # region 4 held out: 0.5 finite on 0-3 and best


def test_real_xenium_wiring_passes_setting_to_run_xenium_qp(monkeypatch):
    seen = []
    monkeypatch.setattr(B, "run_xenium_qp", lambda g, c, prof, region, setting=None: seen.append((region, setting)))
    qp_fn = R._build_wiring("xenium", 6)["qp_fn"]
    cond = {"axis": "setting", "param": "lambda_reg", "level": 2}
    qp_fn("g", "c", PROFILE, "3", cond)
    qp_fn("g", "c", PROFILE, "3")
    assert seen == [(3, cond), (3, None)]
