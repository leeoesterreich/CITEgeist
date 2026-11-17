"""Tests for gt_rebuild/integrity_gates.py (synthetic data only)."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import gt_config  # noqa: E402
from gt_rebuild import integrity_gates as ig  # noqa: E402

pytestmark = pytest.mark.unit

STEMS = {
    "B cells": "B_cells",
    "CD4+ T cells": "CD4pos_T_cells",
    "CD8+ T cells": "CD8pos_T_cells",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
    "Macrophages": "Macrophages",
}
N_SPOTS = 5
GATES = [
    "type_has_cluster",
    "cell_id_join_and_duplicates",
    "mapped_counts_match_mapping",
    "rows_sum_to_one",
    "type_columns_present",
    "five_regions",
    "spot_order_consistent",
    "type_evaluable_in_3_regions",
    "gex_layers_match_canonical_spots",
]
COLLAPSE = {"CD4+ T cells": "T cells", "CD8+ T cells": "T cells"}


def _props(cells: pd.DataFrame, types):
    counts = cells.groupby(["spot_id", "label"]).size().unstack(fill_value=0)
    counts = counts.reindex(index=[f"spot_{i}" for i in range(N_SPOTS)], columns=types, fill_value=0)
    out = counts.div(counts.sum(axis=1), axis=0).fillna(0)
    out["n_cells"] = counts.sum(axis=1)
    out["spot_x"] = np.arange(N_SPOTS)
    out["spot_y"] = np.arange(N_SPOTS)
    out.index.name = "spot_id"
    return out


@pytest.fixture
def world(tmp_path):
    rng = np.random.default_rng(0)
    rows = []
    for r in range(5):
        for k in range(70):
            rows.append((f"r{r}c{k}", r, int(rng.integers(0, N_SPOTS)), gt_config.CT7[int(rng.integers(0, 7))]))
    cells = pd.DataFrame(rows, columns=["cell_id", "region_id", "spot_idx", "label"])
    cells["spot_id"] = [f"spot_{i}" for i in cells["spot_idx"]]
    gt = tmp_path / "gt"
    for n, types in ((6, gt_config.CT6), (7, gt_config.CT7)):
        d = gt / f"ground_truth_{n}type" / "ground_truth"
        d.mkdir(parents=True)
        for r in range(5):
            c = cells[cells["region_id"] == r].copy()
            if n == 6:
                c["label"] = c["label"].replace(COLLAPSE)
            _props(c, types).to_csv(d / f"Xenium_region_{r}_prop.csv")
    for r in range(5):
        d = gt / "ground_truth_7type" / "ground_truth_gex" / f"Xenium_region_{r}"
        d.mkdir(parents=True)
        for stem in STEMS.values():
            pd.DataFrame(
                rng.random((N_SPOTS, 3)),
                index=pd.Index([f"spot_{i}" for i in range(N_SPOTS)], name="spot_id"),
                columns=["G1", "G2", "G3"],
            ).to_csv(d / f"{stem}_GT.csv")
    # mapping: all labelled cells + 3 unmapped labelled cells + 2 dropped (zero-count) cells with no label
    mp = cells[["cell_id", "region_id", "spot_idx", "spot_id"]].set_index("cell_id")
    extra = pd.DataFrame(
        {
            "region_id": [0, np.nan, 1, 2, 3],
            "spot_idx": [-1, -1, 0, 1, 2],
            "spot_id": [None, None, "spot_0", "spot_1", "spot_2"],
        },
        index=["u0", "u1", "z0", "z1", "z2"],
    )
    mp = pd.concat([mp, extra])
    mp.index.name = None
    mp_path = tmp_path / "mapping.csv"
    mp.to_csv(mp_path)
    dropped = ["z0", "z1", "z2"]
    lab = cells[["cell_id", "label"]].rename(columns={"label": "singler_label"})
    lab["cluster"] = "0"
    lab = pd.concat(
        [lab, pd.DataFrame({"cell_id": ["u0", "u1"], "singler_label": ["B cells"] * 2, "cluster": ["0"] * 2})]
    )
    lab_path = tmp_path / "cell_labels_7type.csv"
    lab.to_csv(lab_path, index=False)
    names = pd.DataFrame({"cluster": [str(i) for i in range(6)], "label": gt_config.CT6, "n_cells": 10})
    return dict(gt=gt, labels=lab_path, mapping=mp_path, names=names, dropped=dropped, tmp=tmp_path, lab=lab, mp=mp)


def _run(w, **over):
    kw = dict(
        gt_root=w["gt"], labels_csv=w["labels"], mapping_csv=w["mapping"], names=w["names"], dropped_ids=w["dropped"]
    )
    kw.update(over)
    return ig.run_gates(**kw)


def _only_fails(res, name):
    assert [g.name for g in res] == GATES
    failed = [g.name for g in res if not g.ok]
    assert failed == [name], [(g.name, g.detail) for g in res if not g.ok]


def test_all_nine_gates_always_evaluated(world):
    res = _run(world)
    assert [g.name for g in res] == GATES
    assert all(g.ok for g in res), [(g.name, g.detail) for g in res if not g.ok]


def test_all_nine_gates_evaluated_even_on_empty_gt(world, tmp_path):
    res = _run(world, gt_root=tmp_path / "nothing")
    assert [g.name for g in res] == GATES
    assert not all(g.ok for g in res)


def test_zero_clusters_for_a_type_stops(world):
    names = world["names"][world["names"]["label"] != "Fibroblasts"]
    _only_fails(_run(world, names=names), "type_has_cluster")


def test_duplicate_cell_id_stops(world):
    lab = pd.concat([world["lab"], world["lab"].iloc[[0]]])
    lab.to_csv(world["labels"], index=False)
    _only_fails(_run(world), "cell_id_join_and_duplicates")


def test_mapped_count_off_by_one_stops(world):
    # one mapped labelled cell loses its label: mapped count is one below the mapping-derived count
    world["lab"][world["lab"]["cell_id"] != "r1c0"].to_csv(world["labels"], index=False)
    _only_fails(_run(world), "mapped_counts_match_mapping")


def test_row_sum_0_99_stops(world):
    p = world["gt"] / "ground_truth_7type" / "ground_truth" / "Xenium_region_2_prop.csv"
    df = pd.read_csv(p, index_col=0)
    assert df.loc["spot_1", "n_cells"] > 0
    df.loc["spot_1", gt_config.CT7] = df.loc["spot_1", gt_config.CT7] * 0.99
    df.to_csv(p)
    _only_fails(_run(world), "rows_sum_to_one")


def test_missing_type_column_stops(world):
    p = world["gt"] / "ground_truth_7type" / "ground_truth" / "Xenium_region_0_prop.csv"
    pd.read_csv(p, index_col=0).drop(columns=["CD8+ T cells"]).to_csv(p)
    res = _run(world)
    assert [g.name for g in res] == GATES
    assert not {g.name: g.ok for g in res}["type_columns_present"]


def test_four_regions_stops(world):
    for n in (6, 7):
        (world["gt"] / f"ground_truth_{n}type" / "ground_truth" / "Xenium_region_4_prop.csv").unlink()
    _only_fails(_run(world), "five_regions")


def test_spot_order_differs_between_6_and_7_stops(world):
    p = world["gt"] / "ground_truth_6type" / "ground_truth" / "Xenium_region_3_prop.csv"
    pd.read_csv(p, index_col=0).iloc[::-1].to_csv(p)
    _only_fails(_run(world), "spot_order_consistent")


def test_type_evaluable_in_two_regions_stops(world):
    for r in (0, 1, 2):
        p = world["gt"] / "ground_truth_7type" / "ground_truth" / f"Xenium_region_{r}_prop.csv"
        df = pd.read_csv(p, index_col=0)
        df["CD4+ T cells"] = df["CD4+ T cells"] + df["CD8+ T cells"]
        df["CD8+ T cells"] = 0.0
        df.to_csv(p)
    res = _run(world)
    assert [g.name for g in res] == GATES
    failed = [g.name for g in res if not g.ok]
    assert failed == ["type_evaluable_in_3_regions"], failed


def test_gex_layer_spot_count_differs_stops(world):
    p = world["gt"] / "ground_truth_7type" / "ground_truth_gex" / "Xenium_region_1" / "Epithelial_GT.csv"
    pd.read_csv(p, index_col=0).iloc[:-1].to_csv(p)
    _only_fails(_run(world), "gex_layers_match_canonical_spots")


def test_write_report_default_path_and_all_ok(world, monkeypatch, tmp_path):
    monkeypatch.setattr(gt_config, "OUT_ROOT", tmp_path / "out")
    res = _run(world)
    path = ig.write_report(res)
    assert path == tmp_path / "out" / "reports" / "integrity_report.json"
    rep = json.loads(path.read_text())
    assert rep["all_ok"] is True and [g["name"] for g in rep["gates"]] == GATES
    res[0] = ig.GateResult(res[0].name, False, "x")
    p2 = ig.write_report(res, tmp_path / "o2" / "r.json")
    assert json.loads(p2.read_text())["all_ok"] is False


def test_report_only(world):
    names = pd.DataFrame(
        {
            "cluster": ["0", "1", "2"],
            "label": ["B cells", "Epithelial", "Fibroblasts"],
            "n_cells": [10, 20, 30],
            "delta_next": [0.3, 0.1, 0.2],
        }
    )
    rep = ig.report_only(names, world["gt"], qc={"n_genes": 405, "n_cells_total": 99})
    assert rep["cluster_sizes"] == {"0": 10, "1": 20, "2": 30}
    assert rep["detected"] == {"n_genes": 405, "n_cells_total": 99}
    d = rep["singler_deltas"]
    assert [x["cluster"] for x in d] == ["1", "2", "0"]
    assert [round(x["cumulative_cell_fraction"], 4) for x in d] == [round(20 / 60, 4), round(50 / 60, 4), 1.0]
    assert set(rep["region_type_presence"]) == {str(r) for r in range(5)}
    assert set(rep["region_type_presence"]["0"]) == set(gt_config.CT7)


def test_integrity_gates_does_not_import_protein_report():
    code = (
        "import sys; sys.path.insert(0, %r); import gt_rebuild.integrity_gates; "
        "assert 'gt_rebuild.quality_report' not in sys.modules and 'quality_report' not in sys.modules"
    ) % str(REPO / "Benchmarking" / "xenium_benchmarking")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_raising_gate_is_reported_failed(world, monkeypatch):
    def boom(*a, **k):
        raise RuntimeError("kaboom")

    monkeypatch.setattr(ig, "_gate_rows_sum", boom)
    res = _run(world)
    assert [g.name for g in res] == GATES
    bad = [g for g in res if not g.ok]
    assert [g.name for g in bad] == ["rows_sum_to_one"] and "kaboom" in bad[0].detail


def test_missing_labels_file_returns_all_nine_failed_dependents(world, tmp_path):
    res = _run(world, labels_csv=tmp_path / "nope.csv")
    assert [g.name for g in res] == GATES
    d = {g.name: g for g in res}
    assert not d["cell_id_join_and_duplicates"].ok and "labels_csv" in d["cell_id_join_and_duplicates"].detail
    assert not d["mapped_counts_match_mapping"].ok
    assert d["rows_sum_to_one"].ok
    assert ig.json.loads(ig.write_report(res, tmp_path / "r.json").read_text())["all_ok"] is False


def test_labels_without_cell_id_returns_all_nine(world):
    pd.DataFrame({"x": [1]}).to_csv(world["labels"], index=False)
    res = _run(world)
    assert [g.name for g in res] == GATES and not {g.name: g.ok for g in res}["cell_id_join_and_duplicates"]


def test_constant_float_noise_column_not_evaluable(world):
    for r in range(5):
        p = world["gt"] / "ground_truth_7type" / "ground_truth" / f"Xenium_region_{r}_prop.csv"
        df = pd.read_csv(p, index_col=0)
        df["Endothelial"] = 0.2 + np.arange(len(df)) * 1e-15
        df.to_csv(p)
    res = {g.name: g.ok for g in _run(world)}
    assert not res["type_evaluable_in_3_regions"]


def test_six_type_t_cells_evaluable_path(world):
    for r in (0, 1, 2):
        p = world["gt"] / "ground_truth_6type" / "ground_truth" / f"Xenium_region_{r}_prop.csv"
        df = pd.read_csv(p, index_col=0)
        df["Macrophages"] = df["Macrophages"] + df["T cells"]
        df["T cells"] = 0.0
        df.to_csv(p)
    res = _run(world)
    d = {g.name: g for g in res}
    assert not d["type_evaluable_in_3_regions"].ok and "6type T cells" in d["type_evaluable_in_3_regions"].detail


def test_extra_region_file_stops(world):
    d = world["gt"] / "ground_truth_7type" / "ground_truth"
    shutil.copy(d / "Xenium_region_0_prop.csv", d / "Xenium_region_5_prop.csv")
    _only_fails(_run(world), "five_regions")


def test_unreadable_prop_keeps_exception_text(world):
    p = world["gt"] / "ground_truth_7type" / "ground_truth" / "Xenium_region_2_prop.csv"
    p.write_bytes(b"")
    d = {g.name: g for g in _run(world)}
    assert not d["five_regions"].ok and "EmptyDataError" in d["five_regions"].detail


def test_empty_names_frame_fails_type_has_cluster(world):
    names = world["names"].iloc[0:0]
    _only_fails(_run(world, names=names), "type_has_cluster")


def test_row_sum_tolerance_is_absolute(world):
    p = world["gt"] / "ground_truth_7type" / "ground_truth" / "Xenium_region_2_prop.csv"
    df = pd.read_csv(p, index_col=0)
    df.loc["spot_1", gt_config.CT7] = df.loc["spot_1", gt_config.CT7] * (1 - 1e-5)
    df.to_csv(p)
    _only_fails(_run(world), "rows_sum_to_one")


# ----------------------------------------------------------------------------- low-confidence cluster exclusion
def _with_exclusion(w, ids=("r0c0", "r0c1")):
    """Labels without ``ids`` (as the assemble stage writes them) plus an independent excluded_cells.csv."""
    lab = w["lab"][~w["lab"]["cell_id"].isin(ids)]
    lab_path = w["tmp"] / "labels_excl.csv"
    lab.to_csv(lab_path, index=False)
    exc_path = w["tmp"] / "excluded_cells.csv"
    pd.DataFrame({"cell_id": list(ids), "cluster": "9"}).to_csv(exc_path, index=False)
    return lab_path, exc_path


def test_independent_count_includes_exclusions(world):
    lab_path, exc_path = _with_exclusion(world)
    res = _run(world, labels_csv=lab_path, excluded_csv=exc_path)
    assert [g.name for g in res] == GATES and all(g.ok for g in res), [(g.name, g.detail) for g in res if not g.ok]


def test_mutation_labels_dropped_cells_but_exclusion_not_recorded_stops(world):
    lab_path, _ = _with_exclusion(world)
    _only_fails(_run(world, labels_csv=lab_path), "mapped_counts_match_mapping")  # no excluded_csv -> counts disagree


def test_mutation_excluded_cell_still_labelled_stops(world):
    _, exc_path = _with_exclusion(world)
    _only_fails(_run(world, excluded_csv=exc_path), "mapped_counts_match_mapping")  # labels still hold the 2 cells


def test_mutation_one_excluded_id_missing_from_exclusion_file_stops(world):
    lab_path, exc_path = _with_exclusion(world)
    pd.read_csv(exc_path).iloc[:1].to_csv(exc_path, index=False)
    _only_fails(_run(world, labels_csv=lab_path, excluded_csv=exc_path), "mapped_counts_match_mapping")


def test_unreadable_excluded_file_fails_only_the_count_gate(world):
    _only_fails(_run(world, excluded_csv=world["tmp"] / "no_such.csv"), "mapped_counts_match_mapping")


def test_all_excluded_spot_with_zero_cells_is_skipped_by_row_sum_gate(world):
    for r in range(5):
        for n in (6, 7):
            p = ig._prop_path(world["gt"], n, r)
            df = pd.read_csv(p, index_col=0)
            df.loc["spot_0", :] = 0.0  # a spot with n_cells == 0 and all-zero proportions
            df.loc["spot_0", ["spot_x", "spot_y"]] = [0, 0]
            df.to_csv(p)
    res = {g.name: g for g in _run(world)}
    assert res["rows_sum_to_one"].ok, res["rows_sum_to_one"].detail


def test_report_only_records_excluded_fraction_without_gating(world, tmp_path):
    names = pd.DataFrame(
        {
            "cluster": ["0", "1", "9"],
            "label": ["B cells", "Epithelial", "Fibroblasts"],
            "n_cells": [60, 30, 10],
            "delta_next": [0.3, 0.1, 0.2],
            "excluded": [False, False, True],
        }
    )
    exc = tmp_path / "excluded_cells.csv"
    pd.DataFrame({"cell_id": [f"x{i}" for i in range(10)], "cluster": "9"}).to_csv(exc, index=False)
    ps = world["gt"] / "ground_truth_7type" / "excluded_cells_per_spot.csv"
    pd.DataFrame({"spot_id": ["a", "b"], "all_excluded": [True, False]}).to_csv(ps, index=False)
    rep = ig.report_only(names, world["gt"], qc=None, excluded_csv=exc)
    assert rep["exclusion"] == {
        "clusters": ["9"],
        "n_excluded_cells": 10,
        "excluded_fraction": pytest.approx(0.1),
        "n_spots_all_excluded": 1,
    }
