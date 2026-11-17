"""Tests for evaluation/src/generate_singler_gex_gt.py (synthetic data only)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "Benchmarking" / "xenium_benchmarking" / "evaluation" / "src"
sys.path.insert(0, str(SRC))

import generate_singler_gex_gt as g  # noqa: E402

pytestmark = pytest.mark.unit

GENES = ["g1", "g2"]
CANON = ["spot_0", "spot_1", "spot_2", "spot_3"]


def _cells():
    return pd.DataFrame(
        {
            "spot_idx": [1, 1, 3],
            "singler_label": ["B cells", "B cells", "B cells"],
            "g1": [1.0, 2.0, 5.0],
            "g2": [0.0, 3.0, 1.0],
        }
    )


def test_layer_has_full_canonical_spot_set_in_order():
    cells = _cells()
    for label in g.CELL_TYPES:
        layer = g.build_type_layer(cells, label, CANON, GENES)
        assert list(layer.index) == CANON
        assert layer.index.name == "spot_id"
        assert len(layer) == len(CANON)
        assert list(layer.columns) == GENES


def test_absent_spots_are_zero_rows():
    layer = g.build_type_layer(_cells(), "B cells", CANON, GENES)
    assert (layer.loc[["spot_0", "spot_2"]] == 0).all().all()


def test_present_spots_match_plain_groupby_sum():
    cells = _cells()
    layer = g.build_type_layer(cells, "B cells", CANON, GENES)
    old = cells.groupby("spot_idx")[GENES].sum()
    old.index = [f"spot_{i}" for i in old.index]
    pd.testing.assert_frame_equal(layer.loc[old.index], old, check_names=False)


def test_type_absent_in_region_is_all_zero_layer():
    layer = g.build_type_layer(_cells(), "Macrophages", CANON, GENES)
    assert layer.shape == (4, 2)
    assert (layer == 0).all().all()


def test_cell_spot_not_in_canonical_raises():
    cells = _cells()
    cells.loc[0, "spot_idx"] = 99
    with pytest.raises(ValueError, match="spot_99"):
        g.build_type_layer(cells, "B cells", CANON, GENES)


def test_spot_name_matches_mapping_spot_id(tmp_path):
    mp = tmp_path / "map.csv"
    pd.DataFrame({"spot_idx": [3668, 12], "spot_id": ["spot_3668", "spot_12"]}, index=["a", "b"]).to_csv(mp)
    m = pd.read_csv(mp, index_col=0)
    assert all(f"spot_{i}" == s for i, s in zip(m["spot_idx"], m["spot_id"]))


def test_cli_requires_out_root_and_labels():
    r = subprocess.run([sys.executable, str(SRC / "generate_singler_gex_gt.py")], capture_output=True, text=True)
    assert r.returncode != 0
    assert "--out-root" in r.stderr


def test_end_to_end_synthetic(tmp_path):
    import numpy as np

    rng = np.random.default_rng(0)
    props = tmp_path / "props"
    props.mkdir()
    lab, gex, mp = [], [], []
    for r in range(5):
        pd.DataFrame(index=pd.Index(["spot_0", "spot_1", "spot_2"])).to_csv(props / f"Xenium_region_{r}_prop.csv")
        cid = f"c{r}"
        lab.append((cid, "B cells"))
        gex.append((cid, *rng.integers(1, 5, 2)))
        mp.append((cid, 1, "spot_1", r))
    pd.DataFrame(lab, columns=["cell_id", "singler_label"]).assign(cluster=0).to_csv(tmp_path / "l.csv", index=False)
    pd.DataFrame(gex, columns=["cell_id", "g1", "g2"]).to_csv(tmp_path / "gex.csv", index=False)
    pd.DataFrame(mp, columns=["cell_id", "spot_idx", "spot_id", "region_id"]).set_index("cell_id").to_csv(
        tmp_path / "map.csv"
    )
    out = tmp_path / "out"
    cmd = [sys.executable, str(SRC / "generate_singler_gex_gt.py"), "--labels-csv", str(tmp_path / "l.csv")]
    cmd += ["--out-root", str(out), "--props-dir", str(props), "--gex-csv", str(tmp_path / "gex.csv")]
    cmd += ["--mapping-csv", str(tmp_path / "map.csv")]
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert len(list(out.glob("Xenium_region_*/*_GT.csv"))) == 35
    mac = pd.read_csv(out / "Xenium_region_2" / "Macrophages_GT.csv", index_col=0)
    assert list(mac.index) == ["spot_0", "spot_1", "spot_2"] and (mac == 0).all().all()


def test_imports_without_gt_rebuild_on_path(tmp_path):
    import os

    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    r = subprocess.run(
        [sys.executable, "-I", str(SRC / "generate_singler_gex_gt.py"), "--help"],
        capture_output=True,
        text=True,
        env=env,
        cwd=tmp_path,
    )
    assert r.returncode == 0, r.stderr
    assert "gt_rebuild" not in (SRC / "generate_singler_gex_gt.py").read_text()


def _synthetic(tmp_path, bad_region=None):
    props = tmp_path / "props"
    props.mkdir()
    lab, gex, mp = [], [], []
    for r in range(5):
        pd.DataFrame(index=pd.Index(["spot_0", "spot_1"])).to_csv(props / f"Xenium_region_{r}_prop.csv")
        lab.append((f"c{r}", "B cells"))
        gex.append((f"c{r}", 3, 4))
        mp.append((f"c{r}", 99 if r == bad_region else 1, r))
    pd.DataFrame(lab, columns=["cell_id", "singler_label"]).to_csv(tmp_path / "l.csv", index=False)
    pd.DataFrame(gex, columns=["cell_id", "g1", "g2"]).to_csv(tmp_path / "gex.csv", index=False)
    m = pd.DataFrame(mp, columns=["cell_id", "spot_idx", "region_id"]).set_index("cell_id")
    m["spot_id"] = [f"spot_{i}" for i in m["spot_idx"]]
    m.to_csv(tmp_path / "map.csv")
    cmd = [sys.executable, str(SRC / "generate_singler_gex_gt.py"), "--labels-csv", str(tmp_path / "l.csv")]
    cmd += ["--out-root", str(tmp_path / "out"), "--props-dir", str(props), "--gex-csv", str(tmp_path / "gex.csv")]
    cmd += ["--mapping-csv", str(tmp_path / "map.csv")]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_bad_spot_in_last_region_writes_nothing(tmp_path):
    r = _synthetic(tmp_path, bad_region=4)
    assert r.returncode != 0 and "spot_99" in r.stderr
    assert not (tmp_path / "out").exists()


def test_counts_stay_integer_in_output(tmp_path):
    assert _synthetic(tmp_path).returncode == 0
    txt = (tmp_path / "out" / "Xenium_region_0" / "B_cells_GT.csv").read_text().splitlines()
    assert txt[1] == "spot_0,0,0" and txt[2] == "spot_1,3,4"
