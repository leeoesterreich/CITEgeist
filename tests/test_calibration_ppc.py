"""Tests for calibration/ppc_gates.py (Task 16). Toy data only."""
import copy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "simulation_generation" / "calibration"))
import ppc_gates as pg  # noqa: E402

pytestmark = pytest.mark.unit

LABS = ["A", "B", "C"]


def _ok(m):
    return m <= 1.0 + pg.EPS


def test_each_threshold_boundary():
    f = pg._factor_margin
    assert _ok(f(1.25, 1.0, 1.25)) and not _ok(f(1.2501, 1.0, 1.25))
    assert _ok(f(1.0, 1.25, 1.25)) and not _ok(f(1.0, 1.2501, 1.25))  # symmetric
    # zero fraction 0.05 / 0.0501
    real = np.ones(10000)
    sim = np.ones(10000)
    sim[:500] = 0
    assert pg.ppc_entry_margins(real, sim)["zero"] <= 1 + pg.EPS
    sim[:501] = 0
    assert pg.ppc_entry_margins(real, sim)["zero"] > 1 + pg.EPS
    # VMR factor 1.5
    assert _ok(f(1.5, 1.0, 1.5)) and not _ok(f(1.5001, 1.0, 1.5))
    # quantile log-ratio 0.3 (pseudocount 1): (q_sim+1)/(q_real+1) = exp(0.3) boundary
    r = np.full(50, 9.0)
    ok = np.full(50, 10 * np.exp(0.3) - 1)
    bad = np.full(50, 10 * np.exp(0.3001) - 1)
    assert pg.ppc_entry_margins(r, ok)["q50"] <= 1 + 1e-9
    assert pg.ppc_entry_margins(r, bad)["q50"] > 1
    # leakage 0.05 absolute on off/on ratio
    assert pg.leakage_ok(0.10, 1.0, 0.15, 1.0)
    assert not pg.leakage_ok(0.10, 1.0, 0.1501, 1.0)
    assert not pg.leakage_ok(0.1, 0.0, 0.1, 1.0)
    # corr 0.1: 2 markers, matrix diff is a single pair
    rng = np.random.default_rng(0)
    x = rng.normal(size=(200, 2))
    assert pg.corr_margin(x, x) == 0.0
    assert pg.corr_ok(x, x)
    # exact boundary: monkeypatch the matrices so the single pair differs by 0.1 / 0.1001
    import unittest.mock as um
    for d, expect in ((0.1, True), (0.1001, False)):
        mats = iter([np.array([[1, 0.5], [0.5, 1]]), np.array([[1, 0.5 + d], [0.5 + d, 1]])])
        with um.patch.object(pg, "spearman_matrix", lambda _x: next(mats)):
            assert pg.corr_ok(x, x) is expect
    assert pg.corr_margin(x, np.column_stack([x[:, 0], rng.normal(size=200)])) > 1  # real decorrelation


def test_gate_requires_90_percent_and_median():
    mk = lambda npass, n, m_pass=0.5, m_fail=2.0: (
        [{"pass": True, "margin": m_pass}] * npass + [{"pass": False, "margin": m_fail}] * (n - npass))
    assert pg.ppc_gate(mk(90, 100))["passed"]
    assert not pg.ppc_gate(mk(89, 100))["passed"]
    # 90% pass but median margin > 1 (entries flagged pass with margin > 1 are inconsistent input): median rule bites
    e = [{"pass": True, "margin": 1.5}] * 90 + [{"pass": False, "margin": 3.0}] * 10
    g = pg.ppc_gate(e)
    assert g["frac_pass"] == 0.9 and not g["median_pass"] and not g["passed"]
    assert not pg.ppc_gate([])["passed"]
    assert not pg.ppc_gate([{"pass": True, "margin": float("nan")}])["passed"]


def _frozen(k=20.0, mu_scale=1.0):
    labels = LABS
    ch = []
    mus = {"M1": [8.0, 0.4, 0.4], "M2": [0.4, 6.0, 0.4], "M3": [0.4, 0.4, 5.0], "M4": [1.0, 1.0, 1.0]}
    for i, (m, mu) in enumerate(mus.items()):
        st = ["tA", "tB", "tC", gp_NONSPEC][i]
        ch.append({"name": f"{st}_P1", "sim_type": st, "rank": 1, "source_marker": m, "borrowed": False,
                   "jitter": None, "mu": [x * mu_scale for x in mu], "phi": 0.3})
    return {"labels": labels, "k": {lb: k for lb in labels}, "channels": ch,
            "channel_names": [c["name"] for c in ch],
            "sim_type_labels": {"tA": {"A": 1.0}, "tB": {"B": 1.0}, "tC": {"C": 1.0}}}


gp_NONSPEC = pg.gp.NONSPECIFIC


def _table(frozen, n_per, regions=(0, 1), seed=5, scale=1.0, markers=("M1", "M2", "M3", "M4")):
    """Real toy cells drawn from the model with mu scaled by `scale`."""
    pm = pg.ppc_markers(frozen)
    rng = np.random.default_rng(seed)
    cnt, reg, lab = [], [], []
    for r in regions:
        for lb in LABS:
            c = pg.simulate_label(frozen, list(markers), {m: dict(pm[m], mu=pm[m]["mu"] * scale) for m in markers},
                                  lb, n_per, rng)
            cnt.append(c); reg += [r] * n_per; lab += [lb] * n_per
    return pg.fx.CellTable(counts=np.vstack(cnt), cell_id=np.arange(len(reg)).astype(object),
                           region=np.array(reg), label=np.array(lab, dtype=object), markers=list(markers))


def test_ppc_passes_on_data_drawn_from_the_model():
    fz = _frozen()
    res = pg.run_ppc(fz, [0, 1], table=_table(fz, 3000))
    assert res["all_gates_passed"], {k: v for k, v in res["gates"].items() if not v["passed"]}
    assert res["seed"] == [412, 77]


def test_ppc_fails_when_mean_off_by_two():
    fz = _frozen()
    res = pg.run_ppc(fz, [0, 1], table=_table(fz, 3000, scale=2.0))
    assert not res["gates"]["mean"]["passed"]
    assert not res["all_gates_passed"]


def test_ppc_reproducible_with_fixed_seed():
    fz = _frozen()
    t = _table(fz, 300)
    a, b = pg.run_ppc(fz, [0, 1], table=t), pg.run_ppc(fz, [0, 1], table=t)
    assert a["gates"] == b["gates"]


def test_denominators_per_metric():
    fz = _frozen()
    res = pg.run_ppc(fz, [0, 1], table=_table(fz, 50))
    n_reg, n_lab, n_mk = 2, 3, 4
    for g in ("mean", "zero", "vmr", "q25", "q50", "q75", "q95"):
        assert res["gates"][g]["n_entries"] == n_reg * n_lab * n_mk  # each quantile counted separately
    assert res["gates"]["corr"]["n_entries"] == n_reg * n_lab  # (region, label)
    # leakage: M1-M3 each have 2 off-target labels; M4 (nonspecific) has none
    assert res["gates"]["leakage"]["n_entries"] == n_reg * 3 * 2


def test_jitter_removed_for_ppc_markers():
    fz = _frozen()
    fz["channels"][0]["jitter"] = 2.0
    fz["channels"][0]["mu"] = [x * 2.0 for x in fz["channels"][0]["mu"]]
    assert np.allclose(pg.ppc_markers(fz)["M1"]["mu"], [8.0, 0.4, 0.4])


def _cell_protein(fz, meta, seed):
    """Toy generated cell protein: Poisson, on-target channel high for its own cell type (independent of the removed NB generator)."""
    rng = np.random.default_rng(seed)
    types = meta[pg.gp.CELL_TYPE_COLUMN].to_numpy()
    cols = fz["channel_names"]
    mean = np.full((len(meta), len(cols)), 0.5)
    for j, c in enumerate(fz["channels"]):
        if c["sim_type"] != pg.gp.NONSPECIFIC:
            mean[types == c["sim_type"], j] = 6.0
    return pd.DataFrame(rng.poisson(mean), index=meta["Cell"].to_numpy(), columns=cols)


def test_sim_bstar_report_geomean_one():
    fz = _frozen()
    n = 600
    types = np.array(["tA", "tB", "tC"] * (n // 3))
    meta = pd.DataFrame({"Cell": [f"c{i}" for i in range(n)], pg.gp.CELL_TYPE_COLUMN: types})
    rng = np.random.default_rng(0)
    st = pd.DataFrame({"Cell": meta["Cell"], "spot": rng.integers(0, 60, n)})
    sim_markers = {"tA": ["tA_P1"], "tB": ["tB_P1"], "tC": ["tC_P1"], pg.gp.NONSPECIFIC: [f"{pg.gp.NONSPECIFIC}_P1"]}
    prot = _cell_protein(fz, meta, 1)
    rep = pg.sim_bstar_report(fz, {"s0": {"sc_meta": meta, "st_index": st, "cell_protein": prot}}, sim_markers)
    b = rep["sims"]["s0"]["bstar"]
    assert set(b) == {"tA_P1", "tB_P1", "tC_P1"}
    assert np.isclose(np.exp(np.mean(np.log(list(b.values())))), 1.0)
    # b* depends on the true proportions: changing the cell-type composition of the spots changes it
    st2 = st.copy()
    st2["spot"] = np.where(types == "tA", st["spot"] % 10, st["spot"])  # tA cells crowd into 10 spots
    b2 = pg.sim_bstar_report(fz, {"s0": {"sc_meta": meta, "st_index": st2, "cell_protein": prot}},
                             sim_markers)["sims"]["s0"]["bstar"]
    assert any(not np.isclose(b[k], b2[k], rtol=1e-3) for k in b)


def test_missing_marker_empty_stratum_and_uncovered_type_raise():
    fz = _frozen()
    t = _table(fz, 20)
    t2 = pg.fx.CellTable(**{**t.__dict__, "markers": ["M1", "M2", "M3", "X"]})
    with pytest.raises(ValueError, match="M4"):
        pg.run_ppc(fz, [0, 1], table=t2)
    keep = ~((t.region == 1) & (t.label == "C"))
    t3 = pg.fx.CellTable(**{**t.__dict__, "counts": t.counts[keep], "cell_id": t.cell_id[keep],
                            "region": t.region[keep], "label": t.label[keep]})
    with pytest.raises(ValueError, match="zero real cells"):
        pg.run_ppc(fz, [0, 1], table=t3)
    res = pg.run_ppc(fz, [0, 1], table=t3, allow_empty_strata={(1, "C"): "none mapped"})
    assert res["skipped_strata"] == [{"region": 1, "label": "C", "reason": "none mapped"}]
    meta = pd.DataFrame({"Cell": ["a", "b", "c"], pg.gp.CELL_TYPE_COLUMN: ["tA", "tB", "tZ"]})
    st = pd.DataFrame({"Cell": ["a", "b", "c"], "spot": [0, 0, 1]})
    with pytest.raises(ValueError, match="tZ"):
        pg.sim_bstar_report(dict(fz, sim_type_labels={**fz["sim_type_labels"], "tZ": {"A": 1.0}}),
                            {"s": {"sc_meta": meta, "st_index": st,
                                   "cell_protein": pd.DataFrame(1, index=meta["Cell"], columns=fz["channel_names"])}},
                            {"tA": ["tA_P1"], "tB": ["tB_P1"]})


def test_shared_source_marker_disagreement_raises():
    fz = _frozen()
    fz["channels"].append(dict(fz["channels"][0], name="dup", phi=0.31))
    with pytest.raises(ValueError, match="disagree"):
        pg.ppc_markers(fz)


def test_descriptive_resample_report_is_ungated_and_reproducible():
    fz = _frozen()
    fz["sim_type_labels"] = {"tA": ["A"], "tB": ["B"], "tC": ["C"]}
    t = _table(fz, 400)
    r1 = pg.descriptive_resample_report(fz, t, [0, 1])
    r2 = pg.descriptive_resample_report(fz, t, [0, 1])
    assert r1 == r2 and r1["gated"] is False and r1["descriptive"] is True
    assert set(r1["per_region"]) == {0, 1} and set(r1["pooled"]) == set(pg.GATES)
    # resampling real cells reproduces real marginals: with 200 held-out cells per stratum the mean gate passes for most entries
    assert r1["pooled"]["mean"]["frac_pass"] > 0.8
    # tiny strata are skipped and recorded, not silently dropped
    keep = ~((t.region == 1) & (t.label == "C")) | (np.arange(len(t.region)) % 400 == 0)
    t2 = pg.fx.CellTable(**{**t.__dict__, "counts": t.counts[keep], "cell_id": t.cell_id[keep],
                            "region": t.region[keep], "label": t.label[keep]})
    assert pg.descriptive_resample_report(fz, t2, [1])["skipped_strata"] == [{"region": 1, "label": "C", "n_cells": 1}]
    with pytest.raises(ValueError, match="M4"):
        pg.descriptive_resample_report(fz, pg.fx.CellTable(**{**t.__dict__, "markers": ["M1", "M2", "M3", "X"]}), [0])


def test_spearman_matrix_constant_leading_column_not_scalar():
    # scipy's legacy check looks only at columns 0/1: a constant column there returned one scalar NaN (crash in corr_margin)
    rng = np.random.default_rng(0)
    x = rng.poisson(3, (40, 5)).astype(float)
    x[:, 0] = 0
    c = pg.spearman_matrix(x)
    assert c.shape == (5, 5) and np.isnan(c[0]).all() and np.isnan(c[:, 0]).all()
    from scipy.stats import spearmanr
    assert np.allclose(c[1:, 1:], spearmanr(x[:, 1:]).statistic)
    y = rng.poisson(3, (40, 5)).astype(float)
    assert np.isfinite(pg.corr_margin(x, y))


def test_spearman_matrix_two_usable_columns_and_non_degenerate_unchanged():
    from scipy.stats import spearmanr
    rng = np.random.default_rng(1)
    x = rng.poisson(3, (30, 4)).astype(float)
    assert np.array_equal(pg.spearman_matrix(x), spearmanr(x).statistic)  # no constant column: identical to scipy
    x2 = x[:, :2]
    c = pg.spearman_matrix(x2)
    assert c.shape == (2, 2) and c[0, 1] == spearmanr(x2).statistic
    x[:, 2:] = 1.0  # only 2 usable columns of 4
    c = pg.spearman_matrix(x)
    assert c.shape == (4, 4) and np.isfinite(c[0, 1]) and np.isnan(c[2:]).all()
    x[:, :] = 1.0
    assert np.isnan(pg.spearman_matrix(x)).all() and pg.corr_margin(x, x) == float("inf")
