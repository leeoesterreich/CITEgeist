"""Tests for brightness estimation panels and the RNA score helper."""
import logging

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from CITEgeist.model.deconvolution import brightness_estimators as be
from CITEgeist.model.deconvolution.marker_share import estimate_marker_brightness_with_diagnostics


def _toy(genes, n_spots=40, n_extra=200, seed=0):
    rng = np.random.default_rng(seed)
    names = list(genes) + [f"G{i}" for i in range(n_extra)]
    X = rng.poisson(2.0, size=(n_spots, len(names))).astype(np.float32)
    return ad.AnnData(X=X, obs=pd.DataFrame(index=[f"s{i}" for i in range(n_spots)]),
                      var=pd.DataFrame(index=names))


def test_panels_match_spec():
    P = be.ESTIMATION_PANELS
    assert set(P) == {"sim", "breast", "xenium6", "xenium7"}
    assert P["xenium7"]["CD4 T"] == ["CD28", "CTLA4", "FOXP3", "TIGIT"]
    assert P["breast"]["Macrophages"] == ["MS4A7", "CSF1R", "AIF1", "MS4A6A"]
    assert P["sim"]["B"] == ["BANK1", "CD22", "CD37", "TNFRSF13B"]
    assert P["sim"]["Macrophages"] == ["CD14", "MS4A7", "CSF1R", "AIF1"]
    assert P["xenium6"]["T"] == ["TRAC", "CD3G", "CD3D", "CD2"]
    assert P["xenium7"]["Epithelial"] == ["KRT7", "MUC1", "CDH1", "EHF"]
    assert P["breast"]["Cancer"] == ["KRT7", "ELF3", "CLDN4", "MUC1"]
    assert P["sim"]["PVL"] == ["RGS5", "ACTA2", "PDGFRB", "NOTCH3"]
    assert be.EVALUATION_PANEL["Cancer"] == ["EPCAM", "KRT8", "KRT18", "KRT19"]
    assert be.EVALUATION_PANEL["CD8 T"] == ["CD8A", "CD8B", "CD3E"]


def test_panels_disjoint_from_evaluation():
    ev = {g for genes in be.EVALUATION_PANEL.values() for g in genes}
    for ds, panel in be.ESTIMATION_PANELS.items():
        assert panel
        used = [g for genes in panel.values() for g in genes]
        assert not ev & set(used), ds
        assert all(len(g) >= 3 for g in panel.values())
    for ds in ("xenium7", "xenium6", "breast"):
        used = [g for genes in be.ESTIMATION_PANELS[ds].values() for g in genes]
        assert len(used) == len(set(used)), ds


def test_panel_to_profile_map_is_explicit_and_consistent():
    M = be.PANEL_TO_PROFILE_TYPES
    assert set(M) == set(be.ESTIMATION_PANELS)
    for ds, m in M.items():
        assert set(m) == set(be.ESTIMATION_PANELS[ds]), ds
    assert M["xenium7"]["CD4 T"] == ["CD4+ T cells"]
    assert M["xenium6"]["T"] == ["T cells"]
    assert M["breast"]["CD8 T"] == ["CD8_T_Cells"]
    assert M["sim"]["Epithelial"] == ["Cancer Epithelial", "Normal Epithelial"]
    assert be.profile_type_to_panel("sim")["Normal Epithelial"] == "Epithelial"
    assert be.profile_type_to_panel("xenium7")["B cells"] == "B"
    # every profile type is claimed by exactly one panel key
    for ds in M:
        types = [t for ts in M[ds].values() for t in ts]
        assert len(types) == len(set(types)), ds


def test_rna_scores_drops_absent_gene_and_reports():
    panel = {"A": ["A1", "A2", "A3", "A4"], "B": ["B1", "B2", "B3"]}
    gex = _toy(["A1", "A2", "A3", "B1", "B2", "B3"])  # A4 absent
    before = gex.X.copy()
    scores, rep = be.rna_scores(gex, panel)
    assert rep["A"]["dropped"] == ["A4"]
    assert rep["A"]["used"] == ["A1", "A2", "A3"]
    assert rep["A"]["estimable"] is True and rep["B"]["dropped"] == []
    assert list(scores.columns) == ["A", "B"]
    assert scores.index.equals(gex.obs_names)
    assert np.isfinite(scores.values).all()
    np.testing.assert_array_equal(gex.X, before)  # input untouched


def test_rna_scores_two_genes_not_estimable():
    panel = {"A": ["A1", "A2", "A3"], "B": ["B1", "B2", "B3"]}
    gex = _toy(["A1", "A2", "B1", "B2", "B3"])
    scores, rep = be.rna_scores(gex, panel)
    assert rep["A"]["estimable"] is False
    assert scores["A"].isna().all()
    assert rep["B"]["estimable"] is True


def test_rna_score_tracks_signal():
    gex = _toy(["A1", "A2", "A3", "A4"])
    gex.X[:20, :4] += 20
    scores, _ = be.rna_scores(gex, {"A": ["A1", "A2", "A3", "A4"]})
    assert scores["A"].iloc[:20].min() > scores["A"].iloc[20:].max()


def test_combined_proxy_sums_listing_types():
    s = pd.DataFrame({"a": [1.0, 2.0], "b": [10.0, 20.0], "c": [np.nan, np.nan]})
    np.testing.assert_allclose(be.combined_proxy(s, ["a", "b"]), [11.0, 22.0])
    np.testing.assert_allclose(be.combined_proxy(s, ["a", "c"]), [1.0, 2.0])  # non-estimable skipped
    assert np.isnan(be.combined_proxy(s, ["c"])).all()
    assert be.combined_proxy(s, ["c"]).shape == (2,)


# ---------------------------------------------------------------- Task 2: b* estimand
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

pytestmark = pytest.mark.unit

_REPO = Path(__file__).resolve().parents[1]


def _toy_bstar(n=60, seed=1):
    """Noiseless-ish toy: 2 types; A lists mA1,mA2; B lists mB; both list SH (shared)."""
    rng = np.random.default_rng(seed)
    prop = rng.dirichlet([1, 1], size=n)
    prop = pd.DataFrame(prop, columns=["A", "B"], index=[f"s{i}" for i in range(n)])
    vn = ["mA1", "mA2", "mB", "SH", "other"]
    Xc = np.column_stack([
        0.5 * prop["A"], 2.0 * prop["A"], 1.5 * prop["B"] + 0.01,
        0.8 * (prop["A"] + prop["B"]) + 0.3 * prop["A"], rng.random(n)])
    profile = {"A": {"Major": ["mA1", "mA2", "SH"]}, "B": {"Major": ["mB", "SH"]}}
    return Xc, vn, profile, prop


def test_bstar_equals_oracle_script():
    pytest.importorskip("scipy")
    sys.path.insert(0, str(_REPO / "Benchmarking" / "xenium_benchmarking" / "spots_baseline"))
    sys.path.insert(0, str(_REPO / "Benchmarking" / "xenium_benchmarking"))
    try:
        from brightness_oracle_test import oracle_brightness
    except Exception as e:  # oracle script needs the xenium benchmark tree/data imports
        pytest.skip(f"oracle script not importable here: {e}")
    Xc, vn, profile, prop = _toy_bstar()
    ref, _ = oracle_brightness(Xc, vn, profile, prop)
    got, diag = be.effective_brightness(Xc, vn, profile, prop)
    assert set(got) == set(ref)
    for m in ref:
        assert got[m] == pytest.approx(ref[m], rel=1e-12)


def test_bstar_shared_marker_uses_summed_proportion():
    Xc, vn, profile, prop = _toy_bstar()
    # summed A+B is constant (1) in the Dirichlet toy, so rescale rows to vary the sum
    prop2 = prop * np.random.default_rng(3).uniform(0.2, 1.0, size=(len(prop), 1))
    Xc2 = Xc.copy()
    Xc2[:, vn.index("SH")] = 0.8 * (prop2["A"] + prop2["B"]).to_numpy()
    got2, diag2 = be.effective_brightness(Xc2, vn, profile, prop2)
    assert got2["SH"] == pytest.approx(0.8, rel=1e-9)
    assert diag2["SH"]["median_fallback"] is False
    assert diag2["SH"]["types"] == ["A", "B"]


def test_bstar_slope_below_threshold_uses_median_of_others():
    Xc, vn, profile, prop = _toy_bstar()
    Xc[:, vn.index("mB")] = 0.5  # flat -> slope 0
    prop2 = prop * np.random.default_rng(4).uniform(0.2, 1.0, size=(len(prop), 1))
    Xc[:, vn.index("mA1")] = 0.5 * prop2["A"].to_numpy()
    Xc[:, vn.index("mA2")] = 2.0 * prop2["A"].to_numpy()
    Xc[:, vn.index("SH")] = 0.8 * (prop2["A"] + prop2["B"]).to_numpy()
    got, diag = be.effective_brightness(Xc, vn, profile, prop2)
    assert diag["mB"]["median_fallback"] is True
    assert diag["mB"]["raw_slope"] <= 0.001
    assert got["mB"] == pytest.approx(np.median([0.5, 2.0, 0.8]), rel=1e-9)
    assert diag["mA1"]["median_fallback"] is False


def test_bstar_threshold_is_strict():
    n = 20
    prop = pd.DataFrame({"A": np.linspace(0, 1, n), "B": np.linspace(1, 0, n) * 0.5})
    vn = ["m1", "m2"]
    Xc = np.column_stack([0.001 * prop["A"], 0.5 * prop["A"]])
    got, diag = be.effective_brightness(Xc, vn, {"A": {"Major": ["m1", "m2"]}}, prop)
    assert diag["m1"]["raw_slope"] == pytest.approx(0.001, abs=1e-12)
    assert diag["m1"]["median_fallback"] is True  # slope exactly 0.001 is NOT > min_slope
    assert got["m1"] == pytest.approx(0.5)


def test_bstar_excludes_nan_gt_rows():
    Xc, vn, profile, prop = _toy_bstar()
    prop2 = prop * np.random.default_rng(3).uniform(0.2, 1.0, size=(len(prop), 1))
    Xc[:, vn.index("mA1")] = 0.5 * prop2["A"].to_numpy()
    bad = prop2.copy()
    bad.iloc[:5] = np.nan
    Xbad = Xc.copy()
    Xbad[:5] = 99.0  # would wreck the slope if not excluded
    got, diag = be.effective_brightness(Xbad, vn, profile, bad)
    assert diag["mA1"]["n_spots_used"] == len(prop2) - 5
    assert got["mA1"] == pytest.approx(0.5, rel=1e-9)


def test_bstar_absent_type_is_neutral_not_crash():
    # a type absent from the whole sample (all-zero proportion) made np.polyfit scale by 0 -> NaN -> LinAlgError
    Xc, vn, profile, prop = _toy_bstar()
    prop2 = prop * np.random.default_rng(3).uniform(0.2, 1.0, size=(len(prop), 1))
    prop2["B"] = 0.0
    Xc[:, vn.index("mA1")] = 0.5 * prop2["A"].to_numpy()
    Xc[:, vn.index("mA2")] = 2.0 * prop2["A"].to_numpy()
    Xc[:, vn.index("SH")] = 0.8 * prop2["A"].to_numpy()  # SH is listed by A too -> identifiable via A
    got, diag = be.effective_brightness(Xc, vn, profile, prop2)
    assert diag["mB"]["degenerate_abundance"] is True and np.isnan(diag["mB"]["raw_slope"])
    assert diag["mB"]["median_fallback"] is True
    assert got["mB"] == pytest.approx(np.median([0.5, 2.0, 0.8]), rel=1e-9)
    assert "degenerate_abundance" not in diag["mA1"] and diag["mA1"]["median_fallback"] is False
    assert got["mA1"] == pytest.approx(0.5, rel=1e-9)


def test_bstar_constant_abundance_degenerate_and_all_degenerate_raises():
    n = 10
    prop = pd.DataFrame({"A": np.full(n, 0.5), "B": np.linspace(0, 1, n)})
    vn = ["m1", "m2"]
    Xc = np.column_stack([np.linspace(0, 1, n), 2.0 * prop["B"]])
    got, diag = be.effective_brightness(Xc, vn, {"A": {"Major": ["m1"]}, "B": {"Major": ["m2"]}}, prop)
    assert diag["m1"]["degenerate_abundance"] is True and got["m1"] == pytest.approx(2.0)
    with pytest.raises(RuntimeError):  # no identifiable slope anywhere -> no median to borrow
        be.effective_brightness(Xc, vn, {"A": {"Major": ["m1", "m2"]}}, prop.assign(A=0.0))


def test_geomean_normalize_gm_is_one():
    out = be.geomean_normalize({"a": 1.0, "b": 4.0, "c": 16.0})
    assert np.exp(np.mean(np.log(list(out.values())))) == pytest.approx(1.0)
    assert out["b"] == pytest.approx(1.0)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            be.geomean_normalize({"a": 1.0, "b": bad})


def test_neutral_fill_never_one():
    out = be.neutral_fill({"a": 0.2, "b": 0.2, "c": None}, {"c"})
    assert out["c"] == pytest.approx(0.2)
    assert out["a"] == 0.2
    out = be.neutral_fill({"a": 0.2, "b": float("nan"), "c": 0.6, "d": 0.4}, {"b"})
    assert out["b"] == pytest.approx(0.4)
    with pytest.raises(ValueError):
        be.neutral_fill({"a": None}, {"a"})


def test_real_profile_types_covered_by_panel_map():
    sys.path.insert(0, str(_REPO / "Benchmarking" / "xenium_benchmarking"))
    sys.path.insert(0, str(_REPO))
    try:
        import benchmark_constants as C
        from CITEgeist.examples.breast_dictionary import BREAST_PROFILE
    except Exception as e:
        pytest.skip(f"dictionaries not importable: {e}")
    M = be.PANEL_TO_PROFILE_TYPES

    def claimed(ds):
        return {t for ts in M[ds].values() for t in ts}

    assert set(C.MS_V41_7_CELL_PROFILE_DICT) <= claimed("xenium7")
    assert set(C.MS_V41_6_CELL_PROFILE_DICT) <= claimed("xenium6")
    assert set(BREAST_PROFILE) <= claimed("breast")


# ---------------------------------------------------------------- Task 3: RNA-gate candidates
def _gate_toy(n=400, seed=5):
    """Noiseless toy: markers a1 (b=2), a2 (b=0.5) on type A; c1 (b=1) on type C.

    The latent abundance vectors are permutations of the same values so the top-decile medians match.
    """
    rng = np.random.default_rng(seed)
    base = np.sort(rng.uniform(0.1, 1.0, n))
    a = rng.permutation(base)
    c = rng.permutation(base)
    scores = pd.DataFrame({"A": a, "C": c})
    vn = ["a1", "a2", "c1"]
    Xc = np.column_stack([2.0 * a, 0.5 * a, 1.0 * c])
    profile = {"A": {"Major": ["a1", "a2"]}, "C": {"Major": ["c1"]}}
    return Xc, vn, profile, scores, {"A": True, "C": True}


GATE_KEYS = {"used_brightness", "method", "proxy_types", "n_genes_used", "n_spots_used", "fallback",
             "reason", "selected_spot_purity_proxy", "selected_idx", "all_fallback"}


def test_gate_recovers_known_brightness_on_noiseless_toy():
    Xc, vn, profile, scores, est = _gate_toy()
    b, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=False)
    assert np.exp(np.mean(np.log(list(b.values())))) == pytest.approx(1.0, abs=1e-12)
    assert b["a1"] / b["a2"] == pytest.approx(4.0, rel=1e-6)
    assert b["a1"] / b["c1"] == pytest.approx(2.0, rel=1e-6)
    assert set(diag["a1"]) == GATE_KEYS
    assert diag["a1"]["method"] == "rna_gate" and diag["a1"]["fallback"] is False
    assert diag["a1"]["proxy_types"] == ["A"]
    assert diag["a1"]["used_brightness"] == pytest.approx(b["a1"])
    assert diag["a1"]["n_spots_used"] == len(diag["a1"]["selected_idx"]) >= 40
    assert diag["a1"]["selected_spot_purity_proxy"] > 1.0  # top decile of a standardised proxy


def _crosstalk_toy(n=3000, seed=7):
    rng = np.random.default_rng(seed)
    a, bq, c = rng.uniform(0, 1, (3, n))
    # RNA score of A carries B cross-talk; B has no markers of its own but is scored
    scores = pd.DataFrame({"A": a + 3.0 * bq, "B": bq, "C": c})
    vn = ["a1", "c1"]
    Xc = np.column_stack([2.0 * a, 1.0 * c])
    profile = {"A": {"Major": ["a1"]}, "C": {"Major": ["c1"]}}
    return Xc, vn, profile, scores, {"A": True, "B": True, "C": True}


def test_gate_adj_removes_competitor_composition():
    Xc, vn, profile, scores, est = _crosstalk_toy()
    bu, _ = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=False)
    ba, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=True)
    truth = 2.0
    assert abs(bu["a1"] / bu["c1"] - truth) / truth > 0.15  # unadjusted gate picks wrong spots
    assert abs(ba["a1"] / ba["c1"] - truth) / truth < 0.10
    assert diag["a1"]["method"] == "rna_gate_adj"


def test_gate_adj_dedupes_identical_competitor_columns():
    Xc, vn, profile, scores, est = _crosstalk_toy()
    scores = scores.copy()
    scores["B2"] = scores["B"]  # identical panel (e.g. Cancer / Normal Epithelial) -> would be singular
    est = {**est, "B2": True}
    b1, _ = be.est_rna_gate(Xc, vn, profile, scores.drop(columns="B2"), {k: v for k, v in est.items() if k != "B2"},
                            adjusted=True)
    b2, _ = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=True)
    assert b2["a1"] == pytest.approx(b1["a1"], rel=1e-9)


def test_gate_fallback_fewer_than_20_spots():
    Xc, vn, profile, scores, est = _gate_toy(n=300)
    # a2 is NaN on all but 10 of the spots in its own top decile -> <20 usable spots
    top = np.argsort(scores["A"].to_numpy())[-30:]
    Xc = Xc.copy()
    Xc[top[:20], vn.index("a2")] = np.nan
    b, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=False)
    assert diag["a2"]["fallback"] is True and diag["a2"]["reason"] == "n_spots<20"
    assert diag["a1"]["fallback"] is False
    assert np.isfinite(b["a2"]) and b["a2"] > 0 and b["a2"] != 1.0
    # neutral = median of the other markers' raw b (a1 raw 2*med, c1 raw med) before geomean-1
    raw = {m: diag[m]["used_brightness"] for m in b}
    assert raw["a2"] == pytest.approx(np.median([raw["a1"], raw["c1"]]), rel=1e-9)


def test_gate_150_spots_all_fallback_is_uniform_no_exception():
    Xc, vn, profile, scores, est = _gate_toy(n=150)
    for adj in (False, True):
        b, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=adj)
        assert all(v == 1.0 for v in b.values())
        assert all(d["fallback"] and d["reason"] == "n_spots<20" and d["all_fallback"] for d in diag.values())


def test_gate_non_estimable_type_uses_neutral():
    Xc, vn, profile, scores, est = _gate_toy()
    scores = scores.copy()
    scores["C"] = np.nan
    b, diag = be.est_rna_gate(Xc, vn, profile, scores, {"A": True, "C": False}, adjusted=True)
    assert diag["c1"]["fallback"] is True and diag["c1"]["reason"] == "no_estimable_type"
    assert diag["c1"]["proxy_types"] == []
    assert all(np.isfinite(v) and v > 0 for v in b.values())
    assert b["c1"] != 1.0


def test_gate_geomean_one_over_all_markers_when_some_fall_back():
    rng = np.random.default_rng(11)
    n = 400
    lat = {t: rng.uniform(0.1, 1, n) for t in "ABCD"}
    scores = pd.DataFrame(lat)
    scores["D"] = np.nan  # D not estimable -> its marker falls back
    vn = ["a", "b", "c", "d"]
    Xc = np.column_stack([1.0 * lat["A"], 3.0 * lat["B"], 0.7 * lat["C"], 5.0 * lat["D"]])
    profile = {t: {"Major": [t.lower()]} for t in "ABCD"}
    est = {"A": True, "B": True, "C": True, "D": False}
    for adj in (False, True):
        b, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=adj)
        assert diag["d"]["fallback"] is True
        assert sum(d["fallback"] for d in diag.values()) == 1
        assert np.exp(np.mean(np.log(list(b.values())))) == pytest.approx(1.0, abs=1e-12)


def test_gate_constant_column_flagged_finite():
    Xc, vn, profile, scores, est = _gate_toy()
    Xc = Xc.copy()
    Xc[:, vn.index("a2")] = 0.0
    b, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=False)
    assert diag["a2"]["fallback"] is True and diag["a2"]["reason"] == "nonpositive_b"
    assert np.isfinite(b["a2"]) and b["a2"] > 0 and b["a2"] != 1.0


def test_gate_zero_variance_proxy_flagged():
    Xc, vn, profile, scores, est = _gate_toy()
    scores = scores.copy()
    scores["C"] = 0.3
    b, diag = be.est_rna_gate(Xc, vn, profile, scores, est, adjusted=False)
    assert diag["c1"]["fallback"] is True and diag["c1"]["reason"] == "proxy_constant"
    assert all(np.isfinite(v) and v > 0 for v in b.values())


def _slope_toy(n=300, seed=3):
    rng = np.random.default_rng(seed)
    a, c = rng.uniform(0.1, 1, (2, n))
    scores = pd.DataFrame({"A": a, "C": c})
    vn = ["a1", "a2", "c1"]
    Xc = np.column_stack([2.0 * a + 0.3, 0.5 * a, 1.0 * c])  # intercept on a1: slope must ignore it
    profile = {"A": {"Major": ["a1", "a2"]}, "C": {"Major": ["c1"]}}
    return Xc, vn, profile, scores, {"A": True, "C": True}


def test_slope_recovers_known_brightness():
    Xc, vn, profile, scores, est = _slope_toy()
    for adj in (False, True):
        b, diag = be.est_rna_slope(Xc, vn, profile, scores, est, adjusted=adj)
        assert np.exp(np.mean(np.log(list(b.values())))) == pytest.approx(1.0, abs=1e-12)
        assert b["a1"] / b["a2"] == pytest.approx(4.0, rel=1e-6)
        assert b["a1"] / b["c1"] == pytest.approx(2.0, rel=1e-6)
        assert set(diag["a1"]) == GATE_KEYS
        assert diag["a1"]["method"] == ("rna_slope_adj" if adj else "rna_slope")
        assert diag["a1"]["fallback"] is False and diag["a1"]["all_fallback"] is False
        assert diag["a1"]["selected_spot_purity_proxy"] is None
        assert diag["a1"]["used_brightness"] == pytest.approx(b["a1"])
        assert diag["a1"]["n_spots_used"] == 300


def test_slope_adj_beats_unadjusted_under_leakage():
    rng = np.random.default_rng(21)
    n = 1500
    a, u, c = rng.uniform(0, 1, (3, n))
    bq = 0.5 * (1 - a) + 0.5 * u  # competitor composition correlated with A
    scores = pd.DataFrame({"A": a + 0.15 * bq, "B": bq, "C": c})  # A's RNA score leaks 0.15 of B
    Xc = np.column_stack([2.0 * a, 1.0 * c])
    profile = {"A": {"Major": ["a1"]}, "C": {"Major": ["c1"]}}
    est = {"A": True, "B": True, "C": True}
    err = {}
    for adj in (False, True):
        b, _ = be.est_rna_slope(Xc, ["a1", "c1"], profile, scores, est, adjusted=adj)
        err[adj] = abs(np.log(b["a1"] / b["c1"] / 2.0))
    assert err[True] <= err[False]
    assert err[True] < 1e-3 < err[False]


def test_slope_adj_dedupes_identical_and_own_columns():
    rng = np.random.default_rng(22)
    n = 400
    a, bq, c = rng.uniform(0, 1, (3, n))
    scores = pd.DataFrame({"A": a, "A2": a, "B": bq, "B2": bq, "C": c})
    Xc = np.column_stack([2.0 * a, 1.0 * c])
    profile = {"A": {"Major": ["a1"]}, "C": {"Major": ["c1"]}}
    est = {t: True for t in scores.columns}
    b, diag = be.est_rna_slope(Xc, ["a1", "c1"], profile, scores, est, adjusted=True)
    assert diag["a1"]["fallback"] is False
    assert b["a1"] / b["c1"] == pytest.approx(2.0, rel=1e-6)


def test_slope_negative_slope_neutral_fallback():
    Xc, vn, profile, scores, est = _slope_toy()
    Xc = Xc.copy()
    Xc[:, vn.index("a2")] = -0.5 * scores["A"].to_numpy()
    for adj in (False, True):
        b, diag = be.est_rna_slope(Xc, vn, profile, scores, est, adjusted=adj)
        assert diag["a2"]["fallback"] is True and diag["a2"]["reason"] == "non_positive_slope"
        assert diag["a1"]["fallback"] is False
        assert np.isfinite(b["a2"]) and b["a2"] > 0 and b["a2"] != 1.0
        raw = {m: diag[m]["used_brightness"] for m in b}
        assert raw["a2"] == pytest.approx(np.median([raw["a1"], raw["c1"]]), rel=1e-9)
        assert np.exp(np.mean(np.log(list(b.values())))) == pytest.approx(1.0, abs=1e-12)


def test_slope_zero_variance_proxy_flagged():
    Xc, vn, profile, scores, est = _slope_toy()
    scores = scores.copy()
    scores["C"] = 0.3
    for adj in (False, True):
        b, diag = be.est_rna_slope(Xc, vn, profile, scores, est, adjusted=adj)
        assert diag["c1"]["fallback"] is True and diag["c1"]["reason"] == "proxy_constant"
        assert all(np.isfinite(v) and v > 0 for v in b.values())


def test_slope_non_estimable_few_spots_and_all_fallback():
    Xc, vn, profile, scores, est = _slope_toy()
    s2 = scores.copy()
    s2["C"] = np.nan
    b, diag = be.est_rna_slope(Xc, vn, profile, s2, {"A": True, "C": False}, adjusted=True)
    assert diag["c1"]["reason"] == "no_estimable_type" and b["c1"] != 1.0
    Xs, vs, ps, ss, es = _slope_toy(n=10)
    for adj in (False, True):
        b, diag = be.est_rna_slope(Xs, vs, ps, ss, es, adjusted=adj)
        assert all(v == 1.0 for v in b.values())
        assert all(d["fallback"] and d["reason"] == "n_spots<20" and d["all_fallback"] for d in diag.values())


# ---------------- Task 5: dispatcher ----------------

_PROFILE = {"A": {"Major": ["mA"]}, "B": {"Major": ["mB"]}, "C": {"Major": ["mC"]}}
_RNA_METHODS = ("rna_gate", "rna_gate_adj", "rna_slope", "rna_slope_adj")
_DIAG_KEYS = {"used_brightness", "method", "proxy_types", "n_genes_used", "n_spots_used",
              "fallback", "reason", "selected_spot_purity_proxy", "selected_idx"}


def _dispatch_data(n=120, seed=3):
    rng = np.random.default_rng(seed)
    prop = rng.dirichlet(np.ones(3), size=n)
    bright = np.array([1.0, 2.0, 4.0])
    X = prop * bright + rng.normal(0, 0.01, (n, 3)) + 0.05
    X = np.abs(X)
    names = [f"s{i}" for i in range(n)]
    scores = pd.DataFrame(prop + rng.normal(0, 0.05, prop.shape), columns=list("ABC"), index=names)
    est = {t: True for t in "ABC"}
    return X, ["mA", "mB", "mC"], names, scores, est


def test_gmm_applied_b_identical_to_current():
    X, vn, names, scores, est = _dispatch_data()
    ref, _ = estimate_marker_brightness_with_diagnostics(X, vn, ["mA", "mB", "mC"], seed=1234, max_components=3)
    res = be.estimate_from_scores("gmm", X, vn, _PROFILE, scores, est)
    assert res.applied_b == ref


def test_gmm_reported_b_is_geomean_one_and_proportional():
    X, vn, names, scores, est = _dispatch_data()
    res = be.estimate_from_scores("gmm", X, vn, _PROFILE, scores, est)
    rep = np.array([res.reported_b[m] for m in vn])
    app = np.array([res.applied_b[m] for m in vn])
    assert np.exp(np.mean(np.log(rep))) == pytest.approx(1.0)
    ratio = rep / app
    assert np.allclose(ratio, ratio[0])


def test_plain_is_all_ones():
    X, vn, names, scores, est = _dispatch_data()
    res = be.estimate_from_scores("plain", X, vn, _PROFILE, scores, est)
    assert res.applied_b == {"mA": 1.0, "mB": 1.0, "mC": 1.0}
    assert res.reported_b == res.applied_b
    assert set(res.diag) == {"mA", "mB", "mC"}


def test_unknown_method_raises():
    X, vn, names, scores, est = _dispatch_data()
    with pytest.raises(ValueError, match="unknown"):
        be.estimate_from_scores("nope", X, vn, _PROFILE, scores, est)
    with pytest.raises(ValueError, match="unknown"):
        be.estimate_brightness("nope", X, vn, _PROFILE)


@pytest.mark.parametrize("method", _RNA_METHODS)
def test_diag_schema_has_required_keys(method):
    X, vn, names, scores, est = _dispatch_data()
    res = be.estimate_from_scores(method, X, vn, _PROFILE, scores, est)
    assert set(res.diag) == set(vn)
    for m in vn:
        assert _DIAG_KEYS <= set(res.diag[m]), (method, m)
    assert res.applied_b == res.reported_b
    assert np.exp(np.mean(np.log(list(res.reported_b.values())))) == pytest.approx(1.0)


def _wrapper_inputs(n=40):
    gex = _toy(["g1", "g2", "g3", "g4"], n_spots=n)
    names = list(gex.obs_names)
    rng = np.random.default_rng(5)
    X = np.abs(rng.normal(1, 0.3, (n, 3)))
    panel = {"A": ["g1", "g2", "g3"], "B": ["g2", "g3", "g4"], "C": ["g1", "g3", "g4"]}
    return X, ["mA", "mB", "mC"], names, gex, panel


def test_rna_method_misaligned_obs_raises():
    X, vn, names, gex, panel = _wrapper_inputs()
    bad = list(reversed(names))
    with pytest.raises(ValueError, match="obs|spot_names"):
        be.estimate_brightness("rna_gate", X, vn, _PROFILE, gex=gex, rna_panel=panel, spot_names=bad)
    with pytest.raises(ValueError, match="obs|spot_names"):
        be.estimate_brightness("rna_gate", X, vn, _PROFILE, gex=gex, rna_panel=panel, spot_names=names[:-1])


def test_rna_method_duplicate_obs_raises():
    X, vn, names, gex, panel = _wrapper_inputs()
    dup = names[:-1] + [names[0]]
    gex2 = gex.copy()
    gex2.obs_names = dup
    with pytest.raises(ValueError, match="unique|duplicate"):
        be.estimate_brightness("rna_gate", X, vn, _PROFILE, gex=gex2, rna_panel=panel, spot_names=dup)


def test_missing_gex_raises_unless_allowed(caplog):
    X, vn, names, gex, panel = _wrapper_inputs()
    with pytest.raises(ValueError, match="gex"):
        be.estimate_brightness("rna_slope", X, vn, _PROFILE, gex=None, rna_panel=panel, spot_names=names)
    with caplog.at_level(logging.WARNING):
        res = be.estimate_brightness("rna_slope", X, vn, _PROFILE, gex=None, rna_panel=panel,
                                     spot_names=names, allow_plain_fallback=True)
    assert res.applied_b == {"mA": 1.0, "mB": 1.0, "mC": 1.0}
    assert any("plain" in r.message for r in caplog.records)


def test_panel_set_required_without_explicit_panel():
    X, vn, names, gex, panel = _wrapper_inputs()
    with pytest.raises(ValueError, match="panel_set"):
        be.estimate_brightness("rna_gate", X, vn, _PROFILE, gex=gex, spot_names=names)


def test_wrapper_and_score_api_same_implementation(monkeypatch):
    X, vn, names, gex, panel = _wrapper_inputs()
    calls = []
    real = be.est_rna_gate

    def recorder(*a, **k):
        calls.append(1)
        return real(*a, **k)

    monkeypatch.setattr(be, "est_rna_gate", recorder)
    res_w = be.estimate_brightness("rna_gate", X, vn, _PROFILE, gex=gex, rna_panel=panel, spot_names=names)
    assert len(calls) == 1
    scores, report = be.rna_scores(gex, panel)
    est = {t: r["estimable"] for t, r in report.items()}
    res_s = be.estimate_from_scores("rna_gate", X, vn, _PROFILE, scores, est)
    assert len(calls) == 2
    assert res_w.applied_b == res_s.applied_b
    assert res_w.reported_b == res_s.reported_b


def test_panel_set_maps_panel_keys_to_profile_types():
    # sim panel set: panel key "Epithelial" scores both epithelial profile types
    gex = _toy(["KRT7", "ELF3", "CLDN4", "MUC1", "BANK1", "CD22", "CD37", "TNFRSF13B"], n_spots=40)
    names = list(gex.obs_names)
    rng = np.random.default_rng(2)
    X = np.abs(rng.normal(1, 0.3, (40, 2)))
    prof = {"Cancer Epithelial": {"Major": ["mE"]}, "B-cells": {"Major": ["mB"]}}
    res = be.estimate_brightness("rna_gate", X, ["mE", "mB"], prof, gex=gex, panel_set="sim", spot_names=names)
    assert res.diag["mE"]["proxy_types"] == ["Cancer Epithelial"]
    assert res.diag["mB"]["proxy_types"] == ["B-cells"]


@pytest.mark.parametrize("method", _RNA_METHODS)
@pytest.mark.parametrize("panel_set", [None, "sim"])
def test_wrapper_diag_n_genes_used_is_int(method, panel_set):
    if panel_set is None:
        X, vn, names, gex, panel = _wrapper_inputs()
        res = be.estimate_brightness(method, X, vn, _PROFILE, gex=gex, rna_panel=panel, spot_names=names)
    else:
        gex = _toy(["KRT7", "ELF3", "CLDN4", "MUC1", "BANK1", "CD22", "CD37", "TNFRSF13B"], n_spots=40)
        X = np.abs(np.random.default_rng(2).normal(1, 0.3, (40, 2)))
        prof = {"Cancer Epithelial": {"Major": ["mE"]}, "B-cells": {"Major": ["mB"]}}
        res = be.estimate_brightness(method, X, ["mE", "mB"], prof, gex=gex, panel_set="sim",
                                     spot_names=list(gex.obs_names))
    for d in res.diag.values():
        assert isinstance(d["n_genes_used"], int) and d["n_genes_used"] >= 3
