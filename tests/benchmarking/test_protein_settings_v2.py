"""Tests for settings round 2: new overrides, benchmark-only protein-GMM variants, joint block, inert status."""
import inspect
import itertools
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
from protein_sensitivity import baseline as B  # noqa: E402
from protein_sensitivity import loro as L  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402
from protein_sensitivity import settings_conditions as SC  # noqa: E402
from protein_sensitivity import settings_v2_report as V2  # noqa: E402

pytestmark = pytest.mark.unit

# Independent literal copies of the v2 design (NOT read from the protocol).
OAT = {
    "fusion_mode": ["intersection", "adaptive"],
    "gex_detection_k": [5, 20],
    "gex_detection_min_corr": [0.05, 0.30],
    "gmm_n_components": [3, 4],
    "gmm_base_threshold": [0.3, 0.7],
    "gmm_adaptive_threshold": ["off"],
}
ALPHA, REG, LAP = [0.0, 0.3, 0.5, 0.7, 1.0], [0.1, 0.5, 1, 2, 5], [0.0, 0.01, 0.03, 0.1, 0.3]
GATE = [0.005, 0.01, 0.02, 0.05]  # all four April levels
PROFILE = {"B cells": {"Major": ["CD20"], "Minor": []}, "T cells": {"Major": ["CD3E"], "Minor": []}}
NUC = pd.Series({"spot_1": 3, "spot_2": 5})
PROD = {"method": "qp", "lambda_laplacian": 0.0, "nuclei_counts": NUC}


# ---------------- condition set / protocol ----------------
def test_condition_count_literal_and_set_independent():
    conds = SC.build_list_v2()
    assert len(conds) == 511  # 1 baseline + 11 OAT detection + 499 joint (literal)
    assert sum(c["axis"] == "baseline" for c in conds) == 1
    oat = {(c["param"], str(c["level"])) for c in conds if c["axis"] == "setting"}
    assert oat == {(p, str(l)) for p, ls in OAT.items() for l in ls} and len(oat) == 11
    joint = {(c["alpha"], c["lambda_reg"], c["lambda_laplacian"], c["detection_gate_ub"]) for c in conds if c["axis"] == "joint"}
    expect = set(itertools.product(ALPHA, REG, LAP, GATE)) - {(0.5, 1, 0.0, 0.01)}
    assert joint == expect and len(joint) == 499
    ids = [V2._cid(c) for c in conds]
    assert len(set(ids)) == 511


def test_families_and_nuclei_prior_dropped():
    proto = SC.load_protocol_v2()
    assert "nuclei_prior_lambda" not in proto["settings"] and "nuclei_prior_lambda" in proto["dropped_settings"]
    fam = {p: s["family"] for p, s in proto["settings"].items()}
    assert {p for p, f in fam.items() if f == "benchmark_only_detection"} == {
        "gmm_n_components", "gmm_base_threshold", "gmm_adaptive_threshold",
    }  # fmt: skip
    assert proto["joint_block"]["levels"] == {"alpha": ALPHA, "lambda_reg": REG, "lambda_laplacian": LAP, "detection_gate_ub": GATE}
    assert proto["joint_block"]["decision_rule"]["mean_gain_max"] == 0.02
    assert proto["joint_block"]["decision_rule"]["region_gain_max"] == 0.05
    assert proto["expected_conditions"]["total_per_region_per_ntype"] == 511
    assert proto["r1_duplicates"]["expected_comparisons"] == 80


def test_production_values_sit_on_the_april_grid_baseline_is_one_cell():
    from CITEgeist.model.citegeist_model import CitegeistModel

    sig = inspect.signature(CitegeistModel.run_cell_proportion_model)
    assert sig.parameters["alpha"].default in ALPHA and sig.parameters["alpha"].default == 0.5
    assert sig.parameters["lambda_reg"].default in REG
    assert PROD["lambda_laplacian"] in LAP  # published call value
    assert "alpha_elastic" not in sig.parameters  # name exists only on the per_type_beta path
    assert sig.parameters["detection_gate_ub"].default == 0.01 and 0.01 in GATE  # production value is a level
    assert sig.parameters["fusion_mode"].default == "union"
    assert sig.parameters["gex_detection_k"].default == 10 and sig.parameters["gex_detection_min_corr"].default == 0.15


# ---------------- solver wiring ----------------
class FakeModel:
    last = None
    in_run = None  # hook called inside run_cell_proportion_model

    def __init__(self, **kw):
        FakeModel.last = self
        self.loaded, self.run_kwargs, self.results = None, None, {}

    def filter_gex(self, **kw):
        pass

    preprocess_gex = preprocess_antibody = filter_gex

    def load_cell_profile_dict(self, p):
        self.loaded = p

    def run_cell_proportion_model(self, **kw):
        self.run_kwargs = kw
        if FakeModel.in_run:
            FakeModel.in_run()
        self.results["cell_prop"] = pd.DataFrame({t: [0.5, 0.5] for t in PROFILE}, index=["spot_1", "spot_2"])


@pytest.fixture
def fake(monkeypatch):
    import CITEgeist.model.citegeist_model as cm

    monkeypatch.setattr(cm, "CitegeistModel", FakeModel)
    monkeypatch.setattr(B, "load_nuclei_counts", lambda region: NUC)
    FakeModel.last, FakeModel.in_run = None, None


def _ad():
    g = ad.AnnData(X=np.ones((2, 2)), obs=pd.DataFrame(index=["spot_1", "spot_2"]))
    return g, g.copy()


def _run(setting):
    g, c = _ad()
    return B.run_xenium_qp(g, c, PROFILE, 0, setting=setting)


def test_baseline_kwargs_equal_round1_published_call(fake):
    _run(None)
    assert FakeModel.last.run_kwargs.keys() == PROD.keys() and FakeModel.last.run_kwargs["lambda_laplacian"] == 0.0
    for s in ({"axis": "baseline"}, {"axis": "baseline", "diagnostics": True}):
        _run(s)
        assert FakeModel.last.run_kwargs.keys() == PROD.keys() and FakeModel.last.run_kwargs["nuclei_counts"] is NUC


@pytest.mark.parametrize("param,level", [(p, l) for p in ("fusion_mode", "gex_detection_k", "gex_detection_min_corr") for l in OAT[p]])
def test_production_oat_override_reaches_solver(fake, param, level):
    _run({"axis": "setting", "param": param, "level": level})
    assert FakeModel.last.run_kwargs == {**PROD, param: level}


@pytest.mark.parametrize("param,level", [(p, l) for p in ("gmm_n_components", "gmm_base_threshold", "gmm_adaptive_threshold") for l in OAT[p]])
def test_gmm_variants_do_not_add_solver_kwargs(fake, param, level):
    _run({"axis": "setting", "param": param, "level": level})
    assert FakeModel.last.run_kwargs.keys() == PROD.keys()


@pytest.mark.parametrize(
    "a,r,l,g", [(0.0, 0.1, 0.0, 0.01), (1.0, 5, 0.3, 0.02), (0.3, 1, 0.0, 0.005), (0.5, 2, 0.1, 0.01), (0.5, 1, 0.0, 0.02)]
)
def test_joint_combo_override_reaches_solver(fake, a, r, l, g):
    s = {"axis": "joint", "param": "joint", "level": SC.joint_level(a, r, l, g), "alpha": a, "lambda_reg": r,
         "lambda_laplacian": l, "detection_gate_ub": g}  # fmt: skip
    _run(s)
    kw = FakeModel.last.run_kwargs
    assert kw == {"method": "qp", "alpha": a, "lambda_reg": r, "lambda_laplacian": l, "detection_gate_ub": g, "nuclei_counts": NUC}


# ---------------- benchmark-only detection patches: effective arguments + scope ----------------
class SpyGMM:
    calls = []

    def __init__(self, *a, **kw):
        SpyGMM.calls.append(kw)


@pytest.fixture
def spies(monkeypatch):
    from CITEgeist.model.deconvolution import detection as D

    seen = {"detect": []}
    SpyGMM.calls = []

    def spy_detect(*a, **kw):
        seen["detect"].append(kw)
        SpyGMM(n_components=kw.get("_n", None))  # not used
        D.GaussianMixture(n_components=2, covariance_type="full")  # production-style keyword call
        return None

    monkeypatch.setattr(D, "detect_cell_types", spy_detect)
    monkeypatch.setattr(D, "GaussianMixture", SpyGMM)
    return D, seen


def _production_style_detect_call(D):
    # exactly how citegeist_model.py:989-995 calls it: threshold and adaptive_threshold by keyword
    D.detect_cell_types("X", {"a": [0]}, threshold=0.5, log_transform=True, adaptive_threshold=True)


def test_threshold_wrapper_replaces_keyword_not_partial(fake, spies):
    D, seen = spies
    FakeModel.in_run = lambda: _production_style_detect_call(D)
    _run({"axis": "setting", "param": "gmm_base_threshold", "level": 0.3})
    assert seen["detect"][-1]["threshold"] == 0.3 and seen["detect"][-1]["adaptive_threshold"] is True
    _run({"axis": "setting", "param": "gmm_adaptive_threshold", "level": "off"})
    assert seen["detect"][-1]["adaptive_threshold"] is False and seen["detect"][-1]["threshold"] == 0.5


def test_duplicate_keyword_would_be_caught():
    got = {}
    wrapped = SC.make_detect_wrapper(lambda *a, **kw: got.update(kw), {"threshold": 0.7})
    wrapped("X", {}, threshold=0.5, adaptive_threshold=True)  # must not raise and must replace
    assert got["threshold"] == 0.7
    import functools

    got.clear()
    functools.partial(lambda *a, **kw: got.update(kw), threshold=0.7)("X", {}, threshold=0.5)
    assert got["threshold"] == 0.5  # what a partial would do: override silently lost -> the test above would fail


@pytest.mark.parametrize("k", [3, 4])
def test_n_components_replaces_keyword(fake, spies, k):
    D, _ = spies
    FakeModel.in_run = lambda: _production_style_detect_call(D)
    _run({"axis": "setting", "param": "gmm_n_components", "level": k})
    assert SpyGMM.calls[-1]["n_components"] == k


def test_patches_do_not_leak_across_conditions_or_errors(fake, spies):
    D, seen = spies
    before = (D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold)
    FakeModel.in_run = lambda: _production_style_detect_call(D)
    _run({"axis": "setting", "param": "gmm_n_components", "level": 3})
    assert (D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold) == before
    _run({"axis": "baseline", "diagnostics": True})  # baseline AFTER a variant: original values
    assert SpyGMM.calls[-1]["n_components"] == 2 and seen["detect"][-1]["threshold"] == 0.5

    def boom():
        raise RuntimeError("solver failed")

    FakeModel.in_run = boom
    with pytest.raises(RuntimeError):
        _run({"axis": "setting", "param": "gmm_base_threshold", "level": 0.7})
    assert (D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold) == before


def test_stress_test_path_setting_none_installs_no_patches(fake, spies):
    D, _ = spies
    snap = (D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold)
    FakeModel.in_run = lambda: setattr(FakeModel, "inside", (D.detect_cell_types, D.GaussianMixture, D._compute_adaptive_threshold))
    _run(None)
    assert FakeModel.inside == snap


def test_recorder_records_real_adaptive_thresholds():
    from CITEgeist.model.deconvolution import detection as D

    rng = np.random.default_rng(0)
    X = np.vstack([rng.poisson(1, (150, 1)), rng.poisson(30, (50, 1))]).astype(float)
    rec = SC.ThresholdRecorder()
    with SC.detection_context({"threshold": 0.7}, rec):
        D.detect_cell_types(X, {"t": [0]}, threshold=0.5, log_transform=True, adaptive_threshold=True)
    assert len(rec.calls) == 1 and rec.calls[0][1] == 0.7  # base threshold actually replaced


# ---------------- diagnostics ----------------
def test_gex_candidate_counts_and_fallback():
    H = np.array([[0.9, 0.8, 0.7, 0.6, 0.0], [0.1, 0.1, 0.1, 0.1, 0.9]])
    assert SC.gex_candidate_counts(H, k=10, min_corr=0.15) == [4, 1]  # type 1 has 1 candidate -> fallback
    assert SC.gex_candidate_counts(H, k=2, min_corr=0.15) == [2, 1]
    assert SC.gex_candidate_counts(H, k=10, min_corr=0.75) == [2, 1]


def test_collect_diagnostics_hash_changes_with_mask():
    m = np.ones((4, 2), dtype=bool)
    m2 = m.copy()
    m2[0, 0] = False
    a = SC.collect_diagnostics({"detection_fused": m, "gex_detection_corr": np.eye(2)}, {}, None)
    b = SC.collect_diagnostics({"detection_fused": m2, "gex_detection_corr": np.eye(2)}, {}, None)
    assert a["mask_hash_fused"] != b["mask_hash_fused"]
    assert a["mask_hash_fused"] == SC.collect_diagnostics({"detection_fused": m.copy()}, {}, None)["mask_hash_fused"]
    assert a["gex_allTrue_fallback_types"] == [0, 1]


# ---------------- runner / real wiring ----------------
def test_real_wiring_passes_joint_and_setting_to_run_xenium_qp(monkeypatch):
    seen = []
    monkeypatch.setattr(B, "run_xenium_qp", lambda g, c, prof, region, setting=None: seen.append((region, setting)))
    qp = R._build_wiring("xenium", 6)["qp_fn"]
    j = {"axis": "joint", "param": "joint", "level": "a0_r1_l0_g0.02", "alpha": 0.0, "lambda_reg": 1, "lambda_laplacian": 0.0, "detection_gate_ub": 0.02}
    qp("g", "c", PROFILE, "2", j)
    qp("g", "c", PROFILE, "2")
    assert seen == [(2, j), (2, None)]


def test_evaluate_unit_extra_arg_rules_and_diag_callback():
    calls, diags = [], []

    def qp_fn(g, c, prof, unit, *extra):
        calls.append(extra)
        out = pd.DataFrame({"a": [1.0]})
        out.attrs["diagnostics"] = {"mask_hash_fused": "h"}
        return out

    j = {"axis": "joint", "param": "joint", "level": "a0_r1_l0_g0.02", "alpha": 0.0, "lambda_reg": 1, "lambda_laplacian": 0.0, "detection_gate_ub": 0.02}
    rows = R.evaluate_unit(
        "xenium", "0", 6, [{"axis": "baseline"}, {"axis": "baseline", "diagnostics": True}, j],
        loader=lambda u: (ad.AnnData(X=np.ones((2, 1))), ad.AnnData(X=np.ones((2, 1)), var=pd.DataFrame(index=["m"]))),
        qp_fn=qp_fn, sace_fn=None, prop_metrics_fn=lambda p, u: {"pearson_r": 0.5, "rmse": 1, "mae": 1, "jsd": 1},
        gex_metrics_fn=None, profile=PROFILE, save_diag=lambda cid, d: diags.append((cid, d)),
    )  # fmt: skip
    assert calls == [(), ({"axis": "baseline", "diagnostics": True},), (j,)]
    assert R.condition_id(j) == "joint_a0_r1_l0_g0.02" and rows[2]["axis"] == "joint" and rows[2]["level"] == "a0_r1_l0_g0.02"
    assert [c for c, _ in diags] == ["baseline", "baseline", "joint_a0_r1_l0_g0.02"]
    b, a = R.apply_condition(ad.AnnData(X=np.ones((2, 1)), var=pd.DataFrame(index=["m"])), PROFILE, j, np.random.default_rng(0))[:2]
    assert a == PROFILE


# ---------------- inert / headline ----------------
def _metrics(levels_r, nt=6):
    rows = []
    for (axis, param, level), rs in levels_r.items():
        for reg, r in enumerate(rs):
            rows.append(dict(dataset="xenium", unit=str(reg), n_types=nt, axis=axis, param=param, level=level,
                             status="ok", pearson_r=r))  # fmt: skip
    return pd.DataFrame(rows)


def test_inert_level_has_no_robust_verdict_and_leaves_denominators():
    proto = SC.load_protocol_v2()
    base = [0.7] * 5
    df = _metrics({
        ("baseline", "", float("nan")): base,
        ("setting", "gmm_base_threshold", 0.3): base,                 # identical mask in every region -> inert
        ("setting", "gmm_base_threshold", 0.7): [0.7, 0.7, 0.7, 0.7, 0.71],  # one region differs -> estimable
        ("setting", "fusion_mode", "intersection"): [0.6] * 5,        # not robust
        ("setting", "fusion_mode", "adaptive"): [0.69] * 5,           # robust
    })  # fmt: skip
    diag = {}
    for reg in range(5):
        diag[(6, reg, "baseline")] = {"mask_hash_fused": "B"}
        diag[(6, reg, "setting_gmm_base_threshold_0.3")] = {"mask_hash_fused": "B"}
        diag[(6, reg, "setting_gmm_base_threshold_0.7")] = {"mask_hash_fused": "B" if reg < 4 else "X"}
        diag[(6, reg, "setting_fusion_mode_intersection")] = {"mask_hash_fused": "I"}
        diag[(6, reg, "setting_fusion_mode_adaptive")] = {"mask_hash_fused": "A"}
    lv = {(x["param"], str(x["level"])): x for x in V2.summarize_oat(df, diag, 6, proto)}
    inert = lv[("gmm_base_threshold", "0.3")]
    assert inert["status"] == "inert" and inert["robust"] is None
    assert lv[("gmm_base_threshold", "0.7")]["status"] == "ok" and lv[("gmm_base_threshold", "0.7")]["robust"] is True
    assert lv[("fusion_mode", "intersection")]["robust"] is False and lv[("fusion_mode", "adaptive")]["robust"] is True
    # level never run / no diagnostics -> not estimable, never robust
    assert lv[("gex_detection_k", "5")]["status"].startswith("not_estimable") and lv[("gex_detection_k", "5")]["robust"] is None
    h = V2.headline_counts(list(lv.values()))
    assert h["production_estimable_levels"] == 2 and h["production_robust_levels"] == 1  # only fusion_mode adaptive/intersection counted
    assert h["benchmark_only_levels_reported_separately"] == 5


def test_missing_diagnostics_makes_mask_level_not_estimable():
    proto = SC.load_protocol_v2()
    df = _metrics({("baseline", "", float("nan")): [0.7] * 5, ("setting", "fusion_mode", "adaptive"): [0.7] * 5})
    lv = {(x["param"], str(x["level"])): x for x in V2.summarize_oat(df, {}, 6, proto)}
    assert lv[("fusion_mode", "adaptive")]["status"].startswith("not_estimable")


# ---------------- joint LORO + decision rule ----------------
def _jgrid(r):
    return pd.DataFrame(
        [(6, reg, c, v) for c, vs in r.items() for reg, v in enumerate(vs)],
        columns=["n_types", "region", "combo", "pearson_r"],
    )


def test_joint_loro_excludes_heldout_and_tie_goes_to_baseline():
    g = _jgrid({
        "base": [0.70] * 5,
        "tuned": [0.60, 0.80, 0.80, 0.80, 0.80],   # great on 1-4, bad on 0
        "tie": [0.70] * 5,
    })  # fmt: skip
    sel = L.select_joint(g, "base").set_index("held_out")
    assert sel.loc[0, "chosen_combo"] == "tuned" and sel.loc[0, "heldout_r_chosen"] == pytest.approx(0.60)
    assert sel.loc[0, "gain"] == pytest.approx(-0.10)
    assert sel.loc[0, "train_mean_r_chosen"] == pytest.approx(0.80)
    assert sel.loc[1, "chosen_combo"] == "base" or sel.loc[1, "chosen_combo"] == "tuned"
    g2 = _jgrid({"base": [0.7] * 5, "tie": [0.7] * 5})
    assert (L.select_joint(g2, "base")["chosen_combo"] == "base").all()
    g3 = _jgrid({"base": [0.7] * 5, "bad": [0.9, 0.9, 0.9, 0.9, np.nan]})  # NaN in region 4 -> ineligible for h != 4
    s3 = L.select_joint(g3, "base").set_index("held_out")
    assert (s3.loc[[0, 1, 2, 3], "chosen_combo"] == "base").all() and s3.loc[4, "chosen_combo"] == "bad"


def test_joint_decision_rule_boundaries():
    def sel(gains):
        return pd.DataFrame({"n_types": 6, "held_out": range(5), "gain": gains})

    assert L.joint_decision(sel([0.02] * 5))["answered"].iloc[0]  # mean == 0.02 passes
    assert not L.joint_decision(sel([0.021] * 5))["answered"].iloc[0]
    assert not L.joint_decision(sel([0.0, 0.0, 0.0, 0.0, 0.051]))["answered"].iloc[0]  # one region > 0.05
    assert L.joint_decision(sel([0.05, 0.0, 0.0, 0.0, 0.0]))["answered"].iloc[0]
    assert L.joint_decision(sel([-0.2, 0.0, 0.0, 0.0, 0.0]))["answered"].iloc[0]  # loss is not inflation
    assert not L.joint_decision(sel([0.0, 0.0, 0.0, 0.0, np.nan]))["answered"].iloc[0]  # missing -> not answered


def test_joint_grid_baseline_combo_slug():
    proto = SC.load_protocol_v2()
    jb = proto["joint_block"]
    assert SC.joint_level(*(jb["baseline_cell"][k] for k in SC.JOINT_KEYS)) == "a0.5_r1_l0_g0.01"


# ====================================================================================================
# Round 2 rev 3: fail-closed gates, run selection, primary/descriptive outputs
# ====================================================================================================
import json  # noqa: E402

from protein_sensitivity import aggregate as A  # noqa: E402

REGIONS = [str(i) for i in range(5)]
REF = {r: 0.70 + 0.01 * i for i, r in enumerate(REGIONS)}


def _bdf(nt=6, drop=None, status=None, r=None):
    rows = []
    for i, u in enumerate(REGIONS):
        if u == drop:
            continue
        rows.append(dict(dataset="xenium", unit=u, n_types=nt, axis="baseline", status=(status or {}).get(u, "ok"),
                         pearson_r=(r or {}).get(u, REF[u])))  # fmt: skip
    return pd.DataFrame(rows, columns=["dataset", "unit", "n_types", "axis", "status", "pearson_r"])


def test_baseline_gate_passes_only_on_complete_ok_finite_input():
    assert A.check_baseline_gate(_bdf(), 6, REF, 1e-4)["pass"] is True


@pytest.mark.parametrize(
    "case",
    ["missing_last_row", "failed_row", "nan_candidate", "nan_reference", "missing_reference_key", "extra_reference_key", "empty_df", "empty_reference"],
)
def test_baseline_gate_fails_closed(case):
    df, ref = _bdf(), dict(REF)
    if case == "missing_last_row":
        df = _bdf(drop="4")  # delete the LAST expected row
    elif case == "failed_row":
        df = _bdf(status={"2": "failed_x"})
    elif case == "nan_candidate":
        df = _bdf(r={"3": float("nan")})
    elif case == "nan_reference":
        ref["1"] = float("nan")
    elif case == "missing_reference_key":
        ref.pop("4")
    elif case == "extra_reference_key":
        ref["9"] = 0.5
    elif case == "empty_df":
        df = _bdf().iloc[0:0]
    elif case == "empty_reference":
        ref = {}
    assert A.check_baseline_gate(df, 6, ref, 1e-4)["pass"] is False


def test_baseline_gate_detects_drift():
    chk = A.check_baseline_gate(_bdf(r={"0": REF["0"] + 0.001}), 6, REF, 1e-4)
    assert chk["pass"] is False and "0" in chk["mismatches"]


def _r1_frames(drop=None, nan_r2=False, fail_r1=False, shift=0.0):
    proto = SC.load_protocol_v2()
    base = proto["joint_block"]["baseline_cell"]
    r2, r1 = [], []
    for nt in (6, 7):
        for u in REGIONS:
            for p, lvls in proto["r1_duplicates"]["levels"].items():
                for lv in lvls:
                    slug = SC.joint_level(*({**base, p: lv}[k] for k in SC.JOINT_KEYS))
                    r2.append(dict(axis="joint", level=slug, n_types=nt, unit=u, status="ok", pearson_r=0.7))
                    r1.append(dict(axis="setting", param=p, level=str(lv), n_types=nt, unit=u, status="ok", pearson_r=0.7))
    r2, r1 = pd.DataFrame(r2), pd.DataFrame(r1)
    if drop == "r2_last":
        r2 = r2.iloc[:-1]
    if drop == "r1_last":
        r1 = r1.iloc[:-1]
    if nan_r2:
        r2.loc[r2.index[3], "pearson_r"] = float("nan")
    if fail_r1:
        r1.loc[r1.index[5], "status"] = "failed_x"
    if shift:
        r2.loc[r2.index[0], "pearson_r"] += shift
    return r2, r1, proto


def test_r1_equivalence_all_80_present_and_equal():
    r2, r1, proto = _r1_frames()
    eq = V2.r1_equivalence(r2, r1, proto, 1e-4)
    assert eq["pass"] is True and eq["n_expected"] == 80


@pytest.mark.parametrize("kw", [{"drop": "r2_last"}, {"drop": "r1_last"}, {"nan_r2": True}, {"fail_r1": True}, {"shift": 0.001}])
def test_r1_equivalence_fails_closed(kw):
    r2, r1, proto = _r1_frames(**kw)
    assert V2.r1_equivalence(r2, r1, proto, 1e-4)["pass"] is False


def test_r1_equivalence_empty_input_fails():
    r2, r1, proto = _r1_frames()
    assert V2.r1_equivalence(r2.iloc[0:0], r1.iloc[0:0], proto, 1e-4)["pass"] is False


def test_r1_expected_count_anchored_to_protocol_literal():
    r2, r1, proto = _r1_frames()
    proto["r1_duplicates"]["expected_comparisons"] = 79
    assert V2.r1_equivalence(r2, r1, proto, 1e-4)["pass"] is False


def _mkrun(root, nt, reg, md5="M", commit="aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", cond="C", rows="x", suffix="a100"):
    d = Path(root) / f"xenium_{nt}_{reg}_race_T_{suffix}"
    d.mkdir(parents=True)
    (d / "metrics.csv").write_text(rows)
    (d / "run_meta.json").write_text(json.dumps({"settings_protocol_v2_md5": md5, "git_hash": commit, "conditions_sha256": cond}))
    return d


def _all_runs(root, **kw):
    for nt in (6, 7):
        for r in REGIONS:
            _mkrun(root, nt, r, **kw)


def test_select_runs_accepts_good_and_rejects_bad_provenance(tmp_path):
    _all_runs(tmp_path)
    ch, pr = V2.select_runs(tmp_path, REGIONS, "M", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "C")
    assert pr == [] and len(ch) == 10
    for args in (("BAD", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "C"), ("M", ["b" * 40], "C"), ("M", [], "C"), ("M", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "BAD")):
        assert V2.select_runs(tmp_path, REGIONS, args[0], args[1], args[2])[1]


def test_select_runs_missing_unit_and_empty_dir_fail(tmp_path):
    _all_runs(tmp_path)
    import shutil

    shutil.rmtree(next(tmp_path.glob("xenium_7_4_*")))  # delete the LAST expected unit
    assert V2.select_runs(tmp_path, REGIONS, "M", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "C")[1]
    empty = tmp_path / "empty"
    empty.mkdir()
    assert V2.select_runs(empty, REGIONS, "M", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "C")[1]


def test_select_runs_duplicates_identical_ok_different_rejected_never_mtime(tmp_path):
    _all_runs(tmp_path)
    _mkrun(tmp_path, 6, "0", suffix="a100_nvlink")  # byte-identical metrics
    assert V2.select_runs(tmp_path, REGIONS, "M", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "C")[1] == []
    _mkrun(tmp_path, 6, "1", rows="DIFFERENT", suffix="a100_nvlink")
    pr = V2.select_runs(tmp_path, REGIONS, "M", ["aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"], "C")[1]
    assert any("DIFFERENT metrics" in x for x in pr)


def test_gate_failure_writes_only_failure_artifact(tmp_path):
    out = tmp_path / "out"
    rc = V2.run(tmp_path / "no_such_frozen_protocol.json", tmp_path / "r", out)
    assert rc == 1 and sorted(p.name for p in out.iterdir()) == ["settings_v2_gate_failure.json"]
    # frozen protocol present but no authorized commits and no runs: still only the failure artifact
    proto = tmp_path / "frozen.json"
    unauthorized = json.loads(SC.PROTOCOL_V2.read_text())
    unauthorized["authorized_implementation_commits"] = []
    proto.write_text(json.dumps(unauthorized))
    (tmp_path / "r").mkdir()
    out2 = tmp_path / "out2"
    assert V2.run(proto, tmp_path / "r", out2) == 1
    assert sorted(p.name for p in out2.iterdir()) == ["settings_v2_gate_failure.json"]
    assert "authorized" in json.loads((out2 / "settings_v2_gate_failure.json").read_text())["problems"][0]


def test_conditions_hash_matches_what_submit_script_writes():
    proto = SC.load_protocol_v2()
    written = json.dumps(json.loads(SC.build_v2()))  # script: json.loads(CONDS_JSON) then json.dump
    import hashlib

    assert V2.conditions_sha256(proto) == hashlib.sha256(written.encode()).hexdigest()


def test_in_sample_best_descriptive_and_baseline_default():
    g = _jgrid({"base": [0.70] * 5, "better": [0.74] * 5, "worse": [0.6] * 5})
    best = V2.in_sample_best(g, "base")[0]
    assert best["best_combo"] == "better" and best["gain"] == pytest.approx(0.04) and best["baseline_mean_r"] == pytest.approx(0.70)
    g2 = _jgrid({"base": [0.70] * 5, "tie": [0.70] * 5})
    assert V2.in_sample_best(g2, "base")[0]["best_combo"] == "base"


def test_headline_wording_and_effective_thresholds_surface():
    lv = [
        {"param": "a", "level": 1, "family": "production", "status": "ok", "robust": True},
        {"param": "b", "level": 1, "family": "production", "status": "inert", "robust": None},
        {"param": "c", "level": 1, "family": "benchmark_only_detection", "status": "ok", "robust": True},
    ]
    h = V2.headline_counts(lv)
    assert "estimable production-setting levels" in h["phrase"]
    assert (h["production_estimable_levels"], h["production_robust_levels"]) == (1, 1)
    proto = SC.load_protocol_v2()
    df = _metrics({("baseline", "", float("nan")): [0.7] * 5, ("setting", "gmm_base_threshold", 0.3): [0.7] * 5})
    diag = {(6, r, "setting_gmm_base_threshold_0.3"): {"mask_hash_fused": "B", "effective_thresholds": [[0.6, 0.3, 0.3], [0.4, 0.3, 0.4]]} for r in range(5)}
    diag.update({(6, r, "baseline"): {"mask_hash_fused": "B"} for r in range(5)})
    lvl = {(x["param"], str(x["level"])): x for x in V2.summarize_oat(df, diag, 6, proto)}
    assert lvl[("gmm_base_threshold", "0.3")]["effective_thresholds_by_region"]["0"] == [0.3, 0.4]
    assert lvl[("gmm_base_threshold", "0.3")]["status"] == "inert"


# ====================================================================================================
# rev 5: exact commit authorization, pinned reference inputs, fresh output dir
# ====================================================================================================
import hashlib  # noqa: E402

FULL = "d71abb26fe6d5de86cc873f003c7287339166946"


@pytest.mark.parametrize("auth", [[FULL[:7]], [FULL[:12]], [""], [FULL.upper()], [FULL + "0"], [None], [123]])
def test_commit_authorization_rejects_prefixes_empty_and_malformed(tmp_path, auth):
    _all_runs(tmp_path, commit=FULL)
    pr = V2.select_runs(tmp_path, REGIONS, "M", auth, "C")[1]
    assert any("40 lowercase hex" in x for x in pr)


def test_commit_authorization_exact_equality_only(tmp_path):
    _all_runs(tmp_path, commit=FULL)
    assert V2.select_runs(tmp_path, REGIONS, "M", [FULL], "C")[1] == []
    other = "e" * 40
    assert V2.select_runs(tmp_path, REGIONS, "M", [other], "C")[1]  # valid hex but not the run's commit
    # a valid authorization list plus an invalid entry still fails
    assert V2.select_runs(tmp_path, REGIONS, "M", [FULL, ""], "C")[1]


def _ref_proto(tmp_path, gate_text="x", r1_text="y", good=True):
    files = {}
    for i in range(2):
        f = tmp_path / f"g{i}.csv"
        f.write_text(gate_text + str(i))
        files[f"g{i}.csv"] = hashlib.sha256(f.read_bytes()).hexdigest() if good else "0" * 64
    r1 = tmp_path / "r1.csv"
    r1.write_text(r1_text)
    return {"reference_files": {"a11_gate": {"files": files}, "r1_metrics": {"path": "r1.csv", "sha256": hashlib.sha256(r1.read_bytes()).hexdigest()}}}


def test_references_loaded_only_by_pinned_hash(tmp_path):
    (tmp_path / "g0.csv").write_text("unit,n_types\n0,6\n")
    (tmp_path / "g1.csv").write_text("unit,n_types\n1,6\n")
    (tmp_path / "r1.csv").write_text("unit,axis\n0,setting\n")
    proto = {"reference_files": {
        "a11_gate": {"files": {f: hashlib.sha256((tmp_path / f).read_bytes()).hexdigest() for f in ("g0.csv", "g1.csv")}},
        "r1_metrics": {"path": "r1.csv", "sha256": hashlib.sha256((tmp_path / "r1.csv").read_bytes()).hexdigest()},
    }}  # fmt: skip
    gate, r1, pr = V2.load_references(proto, tmp_path)
    assert pr == [] and len(gate) == 2 and len(r1) == 1
    (tmp_path / "g1.csv").write_text("unit,n_types\n1,6\n9,9\n")  # changed after pinning
    assert any("sha256 mismatch" in x for x in V2.load_references(proto, tmp_path)[2])
    (tmp_path / "g1.csv").unlink()
    assert any("file missing" in x for x in V2.load_references(proto, tmp_path)[2])
    assert V2.load_references({}, tmp_path)[2]  # no reference spec at all -> problems
    bad = {"reference_files": {"a11_gate": {"files": {"g0.csv": ""}}, "r1_metrics": {"path": "r1.csv", "sha256": "xyz"}}}
    assert any("missing or malformed" in x for x in V2.load_references(bad, tmp_path)[2])


@pytest.mark.skipif(not (SC.REPO / "Benchmarking/protein_sensitivity/results/settings_metrics.csv").exists(), reason="reference files absent")
def test_draft_protocol_pins_the_real_reference_files():
    proto = SC.load_protocol_v2()
    gate, r1, pr = V2.load_references(proto)
    assert pr == [] and len(proto["reference_files"]["a11_gate"]["files"]) == 10
    base = gate[(gate["axis"] == "baseline") & (gate["dataset"] == "xenium")]
    assert len(base) == 10 and not base.duplicated(["n_types", "unit"]).any()


def test_duplicate_a11_baseline_row_fails_gate():
    proto = SC.load_protocol_v2()
    gate = pd.concat([_bdf(6), _bdf(7), _bdf(6).iloc[[0]]], ignore_index=True)  # duplicate region 0, 6-type
    gate["axis"] = "baseline"
    df = pd.DataFrame({"dataset": ["xenium"], "unit": ["0"], "n_types": [6], "axis": ["baseline"], "status": ["ok"], "pearson_r": [0.7], "level": [np.nan]})
    probs, _ = V2.gate_problems(df, gate, pd.DataFrame(columns=["axis", "param", "n_types", "unit", "level", "status", "pearson_r"]), proto, 1e-4)
    assert any("exactly one baseline row per region" in x for x in probs)


def test_report_refuses_non_empty_output_dir_and_never_leaves_stale_files(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "settings_v2_summary.json").write_text("{stale}")
    rc = V2.run(tmp_path / "no_such.json", tmp_path / "r", out)
    assert rc == 2 and [p.name for p in out.iterdir()] == ["settings_v2_summary.json"]  # nothing added or removed
    out2 = tmp_path / "fresh"
    assert V2.run(tmp_path / "no_such.json", tmp_path / "r", out2) == 1
    assert [p.name for p in out2.iterdir()] == ["settings_v2_gate_failure.json"]
    assert V2.run(tmp_path / "no_such.json", tmp_path / "r", out2) == 2  # second run into the same dir refused
    empty = tmp_path / "empty"
    empty.mkdir()
    assert V2.run(tmp_path / "no_such.json", tmp_path / "r", empty) == 1  # empty existing dir is allowed
