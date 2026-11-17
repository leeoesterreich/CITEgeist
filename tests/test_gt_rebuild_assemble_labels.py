"""Tests for Benchmarking/xenium_benchmarking/gt_rebuild/assemble_labels.py (synthetic data only)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import gt_config  # noqa: E402
from gt_rebuild.assemble_labels import assign_cell_labels, collapse7to6, write_label_files  # noqa: E402

pytestmark = pytest.mark.unit

PY = sys.executable
AGG = REPO / "Benchmarking/xenium_benchmarking/ground_truth_singler/src/aggregate_to_spots.py"


def _inputs():
    cells = [f"c{i}" for i in range(8)]
    level1 = pd.Series(["0", "0", "1", "1", "2", "3", "4", "4"], index=cells, name="cluster")
    names = pd.DataFrame(
        {
            "cluster": ["0", "1", "2", "3", "4"],
            "label": ["Epithelial", "T cells", "B cells", "Macrophages", "Fibroblasts"],
        }
    )
    tsplit = pd.DataFrame(
        {"cell_id": ["c2", "c3"], "subcluster": ["0", "1"], "label7": ["CD4+ T cells", "CD8+ T cells"]}
    )
    return level1, names, tsplit


def test_collapse_is_deterministic_and_total():
    s = pd.Series(gt_config.CT7)
    out = collapse7to6(s)
    assert list(out) == ["B cells", "T cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]
    assert set(out) == set(gt_config.CT6)
    assert collapse7to6(s).equals(out)


def test_unknown_label_raises():
    with pytest.raises(ValueError):
        collapse7to6(pd.Series(["B cells", "T cells"]))  # "T cells" is not a 7-type label
    with pytest.raises(ValueError):
        collapse7to6(pd.Series(["Neutrophils"]))


def test_every_cell_inherits_cluster_label():
    df = assign_cell_labels(*_inputs())
    assert list(df.columns) == ["cell_id", "cluster", "label7", "label6"]
    d = df.set_index("cell_id")
    assert d.loc["c0", "label7"] == "Epithelial" and d.loc["c1", "label6"] == "Epithelial"
    assert d.loc["c4", "label7"] == "B cells" and d.loc["c7", "label7"] == "Fibroblasts"
    assert d.loc["c2", "label7"] == "CD4+ T cells" and d.loc["c3", "label6"] == "T cells"
    assert len(df) == 8


def test_no_unassigned_class():
    df = assign_cell_labels(*_inputs())
    assert df["label7"].notna().all() and df["label6"].notna().all()
    assert set(df["label7"]) <= set(gt_config.CT7)
    assert set(df["label6"]) <= set(gt_config.CT6)
    level1, names, tsplit = _inputs()
    with pytest.raises(ValueError):
        assign_cell_labels(level1, names, tsplit.iloc[:1])  # c3 is a T cell with no split label
    with pytest.raises(ValueError):
        assign_cell_labels(level1, names.iloc[:4], tsplit)  # cluster 4 has no name


def test_duplicate_cell_ids_raise():
    level1, names, tsplit = _inputs()
    dup = pd.concat([level1, level1.iloc[:1]])
    with pytest.raises(ValueError):
        assign_cell_labels(dup, names, tsplit)
    with pytest.raises(ValueError):
        assign_cell_labels(level1, names, pd.concat([tsplit, tsplit.iloc[:1]]))


def test_aggregate_6type_equals_collapse_of_7type(tmp_path):
    df = assign_cell_labels(*_inputs())
    write_label_files(df, tmp_path)
    l7 = pd.read_csv(tmp_path / "cell_labels_7type.csv")
    l6 = pd.read_csv(tmp_path / "cell_labels_6type.csv")
    assert list(l7.columns) == ["cell_id", "singler_label", "cluster"]
    assert list(l6.columns) == ["cell_id", "singler_label", "cluster"]

    # mapping in the real file's format: unnamed index + spot_idx,spot_id,spot_x,spot_y,region_id
    rows = [
        ("c0", 0, "s0", 1.0, 2.0, 0),
        ("c1", 0, "s0", 1.0, 2.0, 0),
        ("c2", 0, "s0", 1.0, 2.0, 0),
        ("c3", 0, "s0", 1.0, 2.0, 0),
        ("c4", 1, "s1", 3.0, 4.0, 0),
        ("c5", 1, "s1", 3.0, 4.0, 0),
        ("c6", 2, "s2", 5.0, 6.0, 1),
        ("c7", -1, None, None, None, None),
    ]
    m = pd.DataFrame(rows, columns=["", "spot_idx", "spot_id", "spot_x", "spot_y", "region_id"]).set_index("")
    mapping = tmp_path / "mapping.csv"
    m.to_csv(mapping)

    out = {}
    for n in (6, 7):
        od = tmp_path / f"agg{n}"
        subprocess.run(
            [
                PY,
                str(AGG),
                "--labels-csv",
                str(tmp_path / f"cell_labels_{n}type.csv"),
                "--mapping-csv",
                str(mapping),
                "--output-dir",
                str(od),
                "--n-types",
                str(n),
            ],
            check=True,
            capture_output=True,
        )
        out[n] = pd.read_csv(od / "ground_truth" / "Xenium_region_0_prop.csv", index_col=0)
    p6, p7 = out[6], out[7]
    assert list(p6.index) == list(p7.index) and len(p6) == 2
    assert ((p6["T cells"] - (p7["CD4+ T cells"] + p7["CD8+ T cells"])).abs() < 1e-12).all()
    assert (p7["CD4+ T cells"] > 0).any()
