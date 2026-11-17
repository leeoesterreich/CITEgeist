"""Strict Xenium evaluator contract (evaluate_benchmark.py): common mask, evaluability, per-type r."""

import importlib.util
import os
import stat
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
SRC = XEN / "evaluation" / "src" / "evaluate_benchmark.py"
sys.path.insert(0, str(XEN))

from gt_rebuild.manifest import sha256_file  # noqa: E402

pytestmark = pytest.mark.unit


def _load():
    spec = importlib.util.spec_from_file_location("eb_under_test", SRC)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


eb = _load()

A4 = np.array([0.1, 0.2, 0.3, 0.4])
SPOTS4 = ["s1", "s2", "s3", "s4"]


def gt_frame(a=A4, spots=SPOTS4):
    a = np.asarray(a, dtype=float)
    df = pd.DataFrame({"A": a, "B": 1 - a}, index=spots)
    df["n_cells"] = 5
    df["spot_x"] = np.arange(len(spots), dtype=float)
    df["spot_y"] = 0.0
    return df


def pred_frame(a=A4, spots=SPOTS4, with_b=True):
    a = np.asarray(a, dtype=float)
    df = pd.DataFrame({"A": a}, index=spots)
    if with_b:
        df["B"] = 1 - a
    return df


def write_gt(root, frames, prefix="Xenium"):
    d = Path(root) / "ground_truth"
    d.mkdir(parents=True, exist_ok=True)
    for i, df in frames.items():
        df.to_csv(d / f"{prefix}_region_{i}_prop.csv")
    return str(root)


def write_pred(root, frames, prefix="Xenium"):
    Path(root).mkdir(parents=True, exist_ok=True)
    for i, df in frames.items():
        df.to_csv(Path(root) / f"{prefix}_region_{i}_cell_prop_finetuned_results.csv")
    return str(root)


def setup(tmp_path, gt, preds, n_regions=2):
    """gt: frame or {region: frame}; preds: {method: frame or {region: frame}}."""
    gts = gt if isinstance(gt, dict) else {i: gt for i in range(n_regions)}
    gdir = write_gt(tmp_path / "gt", gts)
    pdirs = {}
    for m, p in preds.items():
        ps = p if isinstance(p, dict) else {i: p for i in range(n_regions)}
        pdirs[m] = write_pred(tmp_path / f"pred_{m}", ps)
    return gdir, pdirs


def run(gdir, pdirs, types=("A", "B"), n_regions=2, max_drop_frac=0.01):
    types = list(types)
    ev = eb.evaluability_table(gdir, types, n_regions=n_regions)
    masks = eb.common_spot_masks(gdir, pdirs, ev, n_regions=n_regions, max_drop_frac=max_drop_frac)
    return ev, masks, {m: eb.evaluate_strict(gdir, d, types, masks, ev, n_regions=n_regions) for m, d in pdirs.items()}


def test_perfect_prediction_r1_rmse0(tmp_path):
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": pred_frame()})
    ev, masks, res = run(gdir, pdirs)
    s = res["X"]["summary"]
    assert s["A_mean_r"] == pytest.approx(1.0)
    assert s["B_mean_r"] == pytest.approx(1.0)
    assert s["per_type_mean_r"] == pytest.approx(1.0)
    assert s["n_types_scored"] == 2
    assert s["overall_mean_pearson_r"] == pytest.approx(1.0)
    assert s["overall_mean_rmse"] == pytest.approx(0.0, abs=1e-12)
    assert s["n_regions"] == 2 and s["total_spots"] == 8 and s["cell_types"] == ["A", "B"]
    assert s["n_zero_filled_missing_type"] == 0 and s["n_zero_const_pred"] == 0


def test_missing_type_is_zero_column_r0(tmp_path):
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": pred_frame(with_b=False)})
    _, _, res = run(gdir, pdirs)
    s = res["X"]["summary"]
    assert s["B_mean_r"] == 0.0
    assert s["A_mean_r"] == pytest.approx(1.0)
    assert s["per_type_mean_r"] == pytest.approx(0.5)
    assert s["n_zero_filled_missing_type"] == 2
    t = res["X"]["region_type_table"]
    assert set(t[t["type"] == "B"]["reason"]) == {"missing_type_r0"}
    # zero column contributes to error metrics: RMSE of B column is nonzero
    assert res["X"]["regions"][0]["B_rmse"] > 0


def test_constant_gt_excluded_for_every_method(tmp_path):
    g0 = gt_frame()
    g0["A"] = 0.25  # A constant in region 0, B varies
    g0["B"] = A4
    gt = {0: g0, 1: gt_frame()}
    pred0 = pred_frame()
    pred0["B"] = A4
    preds = {"X": {0: pred0, 1: pred_frame()}, "Y": {0: pred0, 1: pred_frame()}}
    gdir, pdirs = setup(tmp_path, gt, preds)
    ev, _, res = run(gdir, pdirs)
    assert not ev.loc[0, "A"] and ev.loc[0, "B"] and ev.loc[1, "A"]
    for m in ("X", "Y"):
        t = res[m]["region_type_table"]
        row = t[(t["region"] == 0) & (t["type"] == "A")].iloc[0]
        assert np.isnan(row["r"]) and row["reason"] == "gt_constant_excluded" and not row["evaluable"]
        assert res[m]["summary"]["A_mean_r"] == pytest.approx(1.0)  # region 1 only


def test_pred_constant_gt_varies_r0(tmp_path):
    p = pred_frame()
    p["A"] = 0.3
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": p})
    _, _, res = run(gdir, pdirs)
    s = res["X"]["summary"]
    assert s["A_mean_r"] == 0.0
    assert s["n_zero_const_pred"] == 2
    t = res["X"]["region_type_table"]
    assert set(t[t["type"] == "A"]["reason"]) == {"pred_constant_r0"}
    assert not np.isnan(s["per_type_mean_r"])


def test_common_mask_shared(tmp_path):
    px = pred_frame().iloc[:3]  # X lacks s4
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": px, "Y": pred_frame()})
    ev, masks, res = run(gdir, pdirs, max_drop_frac=0.5)
    for i in (0, 1):
        assert list(masks[i]) == ["s1", "s2", "s3"]
    for m in ("X", "Y"):
        assert res[m]["summary"]["total_spots"] == 6
    rep = res["Y"]["mask_report"]
    assert rep[0]["n_gt_spots"] == 4
    assert rep[0]["n_mask_spots"] == 3
    assert rep[0]["spots_dropped"] == 1


def test_evaluability_computed_before_predictions_load(tmp_path):
    gdir = write_gt(tmp_path / "gt", {0: gt_frame(), 1: gt_frame()})
    ev = eb.evaluability_table(gdir, ["A", "B"], n_regions=2)
    assert not (tmp_path / "does_not_exist").exists()
    assert ev.shape == (2, 2) and ev.values.all()
    assert list(ev.index) == [0, 1] and list(ev.columns) == ["A", "B"]


def test_method_dropping_varying_spots_does_not_make_type_unevaluable_for_everyone(tmp_path):
    a = np.array([0.25, 0.25, 0.1, 0.4])  # A varies only at s3, s4
    gt = gt_frame(a)
    px = pred_frame(a).iloc[:2]  # X omits every spot where A varies
    gdir, pdirs = setup(tmp_path, gt, {"X": px, "Y": pred_frame(a)})
    ev = eb.evaluability_table(gdir, ["A", "B"], n_regions=2)
    assert ev["A"].all()  # fixed rule: decided on the full GT
    with pytest.raises(eb.MaskDegenerate):
        eb.common_spot_masks(gdir, pdirs, ev, n_regions=2, max_drop_frac=0.6)
    assert ev["A"].all()


def _big(n=200):
    spots = [f"s{i}" for i in range(n)]
    return np.linspace(0.05, 0.6, n), spots


def test_mask_dropping_more_than_1pct_raises(tmp_path):
    a, spots = _big()
    px = pred_frame(a, spots).iloc[3:]  # drops 3/200 = 1.5%
    gdir, pdirs = setup(tmp_path, gt_frame(a, spots), {"X": px, "Y": pred_frame(a, spots)})
    ev = eb.evaluability_table(gdir, ["A", "B"], n_regions=2)
    with pytest.raises(eb.MaskDegenerate):
        eb.common_spot_masks(gdir, pdirs, ev, n_regions=2)


def test_mask_dropping_under_1pct_reports_not_raises(tmp_path):
    a, spots = _big()
    px = pred_frame(a, spots).iloc[1:]  # drops 1/200 = 0.5%
    gdir, pdirs = setup(tmp_path, gt_frame(a, spots), {"X": px, "Y": pred_frame(a, spots)})
    ev, masks, res = run(gdir, pdirs)
    assert len(masks[0]) == 199
    assert res["Y"]["mask_report"][0]["spots_dropped"] == 1


def test_baseline_csv_gz_loaded_natively_and_exact_file_hashed(tmp_path):
    d = tmp_path / "pred" / "Xenium_region_0" / "predictions"
    d.mkdir(parents=True)
    f = d / "baseline.csv.gz"
    pred_frame().to_csv(f, compression="gzip")
    os.chmod(f, 0o444)
    os.chmod(d, 0o555)
    try:
        before = sorted(str(p) for p in (tmp_path / "pred").rglob("*"))
        digest_before = sha256_file(f)
        df, src = eb.load_prediction_with_source(tmp_path / "pred", 0)
        assert src == f
        assert sha256_file(src) == digest_before
        assert list(df.index) == SPOTS4 and list(df.columns) == ["A", "B"]
        df2 = eb.load_prediction(tmp_path / "pred", 0)
        pd.testing.assert_frame_equal(df, df2)
        assert sorted(str(p) for p in (tmp_path / "pred").rglob("*")) == before
        assert not (f.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
    finally:
        os.chmod(d, 0o755)


def test_four_regions_raises(tmp_path):
    gdir = write_gt(tmp_path / "gt", {i: gt_frame() for i in range(4)})
    with pytest.raises(FileNotFoundError):
        eb.evaluability_table(gdir, ["A", "B"], n_regions=5)


def test_missing_region_file_raises_not_skipped(tmp_path):
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": {0: pred_frame()}})
    ev = eb.evaluability_table(gdir, ["A", "B"], n_regions=2)
    with pytest.raises(FileNotFoundError):
        eb.common_spot_masks(gdir, pdirs, ev, n_regions=2)


def test_spot_id_format_mismatch_raises(tmp_path):
    p = pred_frame()
    p.index = [f"AAAC{i}-1" for i in range(4)]
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": p})
    ev = eb.evaluability_table(gdir, ["A", "B"], n_regions=2)
    with pytest.raises(ValueError):
        eb.common_spot_masks(gdir, pdirs, ev, n_regions=2)


def test_flat_r_is_mean_over_regions_of_flattened_r(tmp_path):
    g = {0: gt_frame(), 1: gt_frame([0.4, 0.1, 0.3, 0.2])}
    p0 = pred_frame([0.15, 0.2, 0.35, 0.3])
    p1 = pred_frame([0.3, 0.2, 0.25, 0.1])
    gdir, pdirs = setup(tmp_path, g, {"X": {0: p0, 1: p1}})
    _, _, res = run(gdir, pdirs)

    def flat(gdf, pdf):
        gv = gdf[["A", "B"]].values.flatten()
        pv = pdf[["A", "B"]].values.flatten()
        return np.corrcoef(gv, pv)[0, 1]

    expected = np.mean([flat(g[0], p0), flat(g[1], p1)])
    s = res["X"]["summary"]
    assert s["overall_mean_pearson_r"] == pytest.approx(expected)
    assert s["overall_std_pearson_r"] == pytest.approx(np.std([flat(g[0], p0), flat(g[1], p1)]))


def test_mae_and_jsd_reported(tmp_path):
    p = pred_frame([0.2, 0.2, 0.3, 0.5])
    gdir, pdirs = setup(tmp_path, gt_frame(), {"X": p})
    _, _, res = run(gdir, pdirs)
    s = res["X"]["summary"]
    exp_mae = np.mean(np.abs(np.r_[A4, 1 - A4] - np.r_[p["A"], p["B"]]))
    assert s["overall_mean_mae"] == pytest.approx(exp_mae)
    assert s["overall_mean_jsd"] > 0
    assert "mean_jsd" in res["X"]["regions"][0] and "A_mae" in res["X"]["regions"][0]


# ---- legacy path must be unchanged -------------------------------------------------
LEGACY_GOLDEN = {'summary': {'n_regions': 2,
                 'total_spots': 16,
                 'cell_types': ['A', 'B', 'C'],
                 'A_mean_r': 0.5431216746947476,
                 'A_std_r': 0.01736571202003656,
                 'B_mean_r': -0.010648301968110133,
                 'B_std_r': 0.17823295719596038,
                 'C_mean_r': 0.15085380989039232,
                 'C_std_r': 0.04918540124147275,
                 'overall_mean_pearson_r': 0.16085847660523642,
                 'overall_std_pearson_r': 0.094828725313826,
                 'overall_mean_jsd': 0.24004432103536272,
                 'overall_mean_rmse': 0.23446959857204064,
                 'overall_mean_mae': 0.1951642884626461},
     'regions': [{'A_pearson_r': 0.5604873867147842,
                  'A_pearson_p': 0.14843734495233737,
                  'A_rmse': 0.1234047771664394,
                  'A_mae': 0.10685979260896665,
                  'B_pearson_r': 0.16758465522785027,
                  'B_pearson_p': 0.6916123918214802,
                  'B_rmse': 0.24619869244874112,
                  'B_mae': 0.1985105471821508,
                  'C_pearson_r': 0.10166840864891957,
                  'C_pearson_p': 0.810681275317406,
                  'C_rmse': 0.285694926124023,
                  'C_mae': 0.24726717273510918,
                  'mean_jsd': 0.22471026193308305,
                  'median_jsd': 0.22007780747844663,
                  'overall_pearson_r': 0.2556872019190624,
                  'overall_rmse': 0.22910268876930806,
                  'overall_mae': 0.18421250417540888,
                  'region_id': 0,
                  'n_spots': 8,
                  'n_cell_types': 3,
                  'cell_types': ['A', 'B', 'C'],
                  'non_eval_mass': 0.0},
                 {'A_pearson_r': 0.525755962674711,
                  'A_pearson_p': 0.18080454411781066,
                  'A_rmse': 0.1933907489270605,
                  'A_mae': 0.18442201651241252,
                  'B_pearson_r': -0.18888125916407053,
                  'B_pearson_p': 0.6541806774764378,
                  'B_rmse': 0.20329117287031664,
                  'B_mae': 0.17879961259042276,
                  'C_pearson_r': 0.20003921113186507,
                  'C_pearson_p': 0.6348122442712125,
                  'C_rmse': 0.30632885843591295,
                  'C_mae': 0.2551265891468147,
                  'mean_jsd': 0.2553783801376424,
                  'median_jsd': 0.2504572445597346,
                  'overall_pearson_r': 0.06602975129141042,
                  'overall_rmse': 0.2398365083747732,
                  'overall_mae': 0.2061160727498833,
                  'region_id': 1,
                  'n_spots': 8,
                  'n_cell_types': 3,
                  'cell_types': ['A', 'B', 'C'],
                  'non_eval_mass': 0.0}]}
# (above captured by running the PRE-change evaluate_all_regions on _legacy_fixture)


def _legacy_fixture(tmp_path):
    rng = np.random.default_rng(0)
    spots = [f"s{i}" for i in range(8)]
    gt, pr = {}, {}
    for r in range(2):
        g = pd.DataFrame(rng.dirichlet([1, 2, 3], 8), index=spots, columns=["A", "B", "C"])
        g["n_cells"] = 4
        g["spot_x"] = 1.0
        g["spot_y"] = 2.0
        gt[r] = g
        pr[r] = pd.DataFrame(rng.dirichlet([2, 2, 2], 8), index=spots, columns=["A", "B", "C"])
    return write_gt(tmp_path / "gt", gt), write_pred(tmp_path / "pred", pr)


def test_legacy_evaluate_all_regions_unchanged(tmp_path):
    gdir, pdir = _legacy_fixture(tmp_path)
    res = eb.evaluate_all_regions(gdir, pdir, n_regions=2, use_achievable_7=False)
    assert LEGACY_GOLDEN is not None
    assert set(res) == set(LEGACY_GOLDEN) == {"summary", "regions"}
    assert set(res["summary"]) == set(LEGACY_GOLDEN["summary"])
    for k, v in LEGACY_GOLDEN["summary"].items():
        if isinstance(v, float):
            assert res["summary"][k] == pytest.approx(v, rel=1e-12), k
        else:
            assert res["summary"][k] == v, k
    assert len(res["regions"]) == len(LEGACY_GOLDEN["regions"])
    for got, exp in zip(res["regions"], LEGACY_GOLDEN["regions"]):
        assert set(got) == set(exp)
        for k, v in exp.items():
            if isinstance(v, float):
                assert got[k] == pytest.approx(v, rel=1e-12), k
            else:
                assert got[k] == v, k
