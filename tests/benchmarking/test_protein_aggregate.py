import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "Benchmarking"))
from protein_sensitivity import aggregate as A  # noqa: E402

pytestmark = pytest.mark.unit


def test_seed_then_unit_aggregation():
    df = pd.DataFrame(
        {
            "dataset": "xenium",
            "n_types": 6,
            "axis": "thinning",
            "level": 0.5,
            "panel": "",
            "unit": ["0", "0", "1", "1"],
            "seed": [0, 1, 0, 1],
            "pearson_r": [0.6, 0.8, 0.9, 0.9],
            "status": "ok",
        }
    )
    s = A.summarize(df)
    cell = s["xenium_6"]["thinning"]["0.5"]
    assert np.isclose(cell["pearson_r_mean"], (0.7 + 0.9) / 2)
    assert np.isclose(cell["pearson_r_sd"], np.std([0.7, 0.9], ddof=1))
    assert cell["n_units"] == 2 and cell["status"] == "ok"


def test_failed_status_not_averaged():
    df = pd.DataFrame(
        {
            "dataset": "xenium",
            "n_types": 6,
            "axis": "ablation",
            "level": np.nan,
            "panel": "zero",
            "unit": ["0", "1"],
            "seed": [0, 0],
            "pearson_r": [np.nan, np.nan],
            "status": "failed_no_markers",
        }
    )
    cell = A.summarize(df)["xenium_6"]["ablation"]["zero"]
    assert cell["status"] == "failed_no_markers"
    assert "pearson_r_mean" not in cell


def _row(**kw):
    base = dict(
        dataset="xenium", n_types=7, axis="thinning", level=0.25, panel="", unit="0", seed=0,
        sentinel="", status="ok", pearson_r=0.5, rmse=0.1, mae=0.05, jsd=0.3,
        gex_flat_r=np.nan, gex_nrmse=np.nan, per_type_json="{}",
    )  # fmt: skip
    base.update(kw)
    return base


def test_sace_sentinel_cells_and_sim_keys():
    df = pd.DataFrame(
        [
            _row(sentinel="moderate", gex_flat_r=0.6, gex_nrmse=0.07, per_type_json='{"B": 0.5}'),
            _row(unit="1", sentinel="moderate", gex_flat_r=0.8, gex_nrmse=0.09, per_type_json='{"B": 0.7}'),
            _row(axis="baseline", level=np.nan, sentinel="baseline", gex_flat_r=0.7, gex_nrmse=0.07),
            _row(dataset="sim", n_types=6, unit="mixed:0", sentinel="moderate", gex_flat_r=0.5, gex_nrmse=0.1),
            _row(dataset="sim", n_types=6, unit="high_seg:0", gex_flat_r=np.nan),
        ]
    )
    s = A.summarize(df)
    c = s["xenium_7"]["sace_sentinels"]["moderate"]["thinning"]
    assert np.isclose(c["gex_flat_r_mean"], 0.7) and c["n_units"] == 2
    assert "gex_flat_r_mean" in s["xenium_7"]["sace_sentinels"]["baseline"]["baseline"]
    assert "sim_mixed" in s and "sim_high_seg" in s
    assert "metric_note" in s["sim_mixed"]
    assert "sace_sentinels" not in s["sim_high_seg"] or not s["sim_high_seg"]["sace_sentinels"]
    # per_type averaged
    assert np.isclose(s["xenium_7"]["thinning"]["0.25"]["per_type"]["B"], 0.6)


def test_noise_note_and_xenium6_no_sace():
    # xenium_6 sentinel rows exist but SACE never ran (NaN gex) -> no sace_sentinels block
    df = pd.DataFrame(
        [_row(n_types=6, axis="noise", level=1.0, sentinel="moderate"), _row(n_types=6, axis="baseline", level=np.nan, sentinel="baseline")]
    )
    s = A.summarize(df)
    assert "noise_note" in s["xenium_6"]
    assert "sace_sentinels" not in s["xenium_6"]


def test_partial_failed_cell():
    df = pd.DataFrame(
        [_row(unit="0", pearson_r=0.5), _row(unit="1", pearson_r=np.nan, status="failed_x")]
    )
    c = A.summarize(df)["xenium_7"]["thinning"]["0.25"]
    assert c["status"] == "partial_failed" and c["n_failed_units"] == 1 and c["n_units"] == 1
    assert np.isclose(c["pearson_r_mean"], 0.5)


def test_noise_note_is_dataset_specific(tmp_path):
    tsv = tmp_path / "bg.tsv"
    tsv.write_text(
        "dataset\tunit\tprofile\tmarker\tbackground\n"
        "xenium\t0\txenium6\tCD20\t0.0\n"
        "xenium\t2\txenium6\tPanCK\t0.0\n"
        "xenium\t1\txenium6\tCD3E\t5.0\n"
        "sim\tmixed:1\tsim\tPVL_Protein_2\t0.0\n"
        "sim\thigh_seg:2\tsim\tPlasmablasts_Protein_2\t0.0\n"
    )
    x = A.noise_note("xenium_6", tsv)
    assert "CD20 (units 0)" in x and "PanCK (units 2)" in x and "PVL" not in x
    m = A.noise_note("sim_mixed", tsv)
    assert "PVL_Protein_2 (units mixed:1)" in m and "Plasmablasts" not in m and "CD20" not in m


def _mk_unit(root, name, mtime, with_meta=True, value=0.5):
    import os

    d = root / name
    d.mkdir()
    pd.DataFrame([_row(pearson_r=value)]).to_csv(d / "metrics.csv", index=False)
    if with_meta:
        (d / "run_meta.json").write_text("{}")
        os.utime(d / "run_meta.json", (mtime, mtime))
    return d


def test_collect_dedups_races_and_warns_on_incomplete(tmp_path, caplog):
    _mk_unit(tmp_path, "xenium_7_0_race_A", 1000, value=0.1)
    _mk_unit(tmp_path, "xenium_7_0_race_B", 2000, value=0.9)  # newest wins
    _mk_unit(tmp_path, "xenium_7_1", 1000, value=0.3)  # non-race dir
    _mk_unit(tmp_path, "xenium_7_2_race_A", 1000, with_meta=False)  # crashed
    (tmp_path / "conditions.json").write_text("[]")  # stray file ignored
    with caplog.at_level("WARNING"):
        df = A.collect(tmp_path)
    assert sorted(df["pearson_r"]) == [0.3, 0.9]
    assert any("xenium_7_2_race_A" in r.message for r in caplog.records)


def test_row_level_partial_failed_and_counts():
    ok = dict(status="ok", pearson_r=0.8)
    bad = dict(status="failed_Cell_proportion_validation_failed_for_sample_'Xenium_region_", pearson_r=np.nan)
    df = pd.DataFrame(
        [
            _row(axis="dropout", level=0.1, unit="0", seed=0, **ok),
            _row(axis="dropout", level=0.1, unit="0", seed=1, **ok),
            _row(axis="dropout", level=0.1, unit="1", seed=0, **bad),
            _row(axis="dropout", level=0.1, unit="1", seed=1, **ok),
        ]
    )
    c = A.summarize(df)["xenium_7"]["dropout"]["0.1"]
    assert c["status"] == "partial_failed"
    assert (c["n_rows_attempted"], c["n_rows_ok"], c["n_rows_failed"]) == (4, 3, 1)
    assert c["n_units"] == 2 and c["n_units_contributing"] == 2 and c["n_failed_units"] == 0
    assert np.isclose(c["pearson_r_mean"], 0.8)
    assert c["failed_statuses"] == ["failed_Cell_proportion_validation_failed_for_sample_'Xenium_region...[truncated]"]


def test_all_ok_and_all_failed_unchanged_semantics():
    df = pd.DataFrame([_row(unit="0"), _row(unit="1")])
    c = A.summarize(df)["xenium_7"]["thinning"]["0.25"]
    assert c["status"] == "ok" and c["n_rows_failed"] == 0 and c["n_units_contributing"] == 2
    z = pd.DataFrame([_row(axis="ablation", level=np.nan, panel="zero", status="failed_no_markers", pearson_r=np.nan)])
    cz = A.summarize(z)["xenium_7"]["ablation"]["zero"]
    assert cz["status"] == "failed_no_markers" and cz["n_rows_failed"] == 1 and cz["n_rows_ok"] == 0


def test_key_set_gate_catches_swap_not_just_count():
    conds = {("xenium", 7): [{"axis": "baseline"}, {"axis": "thinning", "level": 0.5, "seed": 0}]}
    units = {("xenium", 7): ["0"]}
    good = pd.DataFrame(
        [_row(axis="baseline", level=np.nan, seed=-1), _row(axis="thinning", level=0.5, seed=0)]
    )
    assert A.check_key_set(good, conds, units) == []
    bad = pd.DataFrame(  # same count (2), but one condition duplicated and one missing
        [_row(axis="baseline", level=np.nan, seed=-1), _row(axis="baseline", level=np.nan, seed=-1)]
    )
    assert A.check_key_set(bad, conds, units)
    assert A.check_key_set(good.iloc[:0], conds, units)  # empty is a failure, not a pass


def test_recompute_dropped_matches_runner_path():
    import anndata as ad  # noqa: F401
    from protein_sensitivity import perturb as P
    from protein_sensitivity import run_sensitivity as R

    prof = {"Epi": {"Major": ["PanCK"]}, "Mac": {"Major": ["CD68", "CD163"]}, "B": {"Major": ["CD20"]}}
    names = ["CD20", "PanCK", "CD68", "CD163"]
    cond = {"axis": "dropout", "level": 0.25, "seed": 1}
    got = A.recompute_dropped(cond, names, prof)
    X = np.zeros((2, 4))
    _, want = P.dropout_channels(X, names, P.qp_active_markers(prof), 0.25, np.random.default_rng(1))
    assert got == want and len(got) == 1


def test_failure_reason_honest_about_truncation():
    prof = {"Epi": {"Major": ["PanCK"]}, "Mac": {"Major": ["CD68", "CD163"]}}
    s = "failed_Cell_proportion_validation_failed_for_sample_'Xenium_region_"
    r = A.failure_reason(s, ["PanCK"], prof)
    assert "truncated" in r and "Epi" in r
    r2 = A.failure_reason(s, ["CD163"], prof)
    assert "not derivable" in r2
    assert A.failure_reason(s, [], prof, "ValueError: full text") == "ValueError: full text"
    assert A.failure_reason("ok", [], prof) == ""
