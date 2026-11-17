"""Tests for Benchmarking/xenium_benchmarking/gt_rebuild/cluster_names.py (synthetic data only)."""

import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import cluster_names as cn  # noqa: E402
from gt_rebuild import gt_config as cfg  # noqa: E402

pytestmark = pytest.mark.unit


def _df():
    rows = []
    for i, lab in enumerate(["B cells", "T cells", "Epithelial"]):
        r = {"cluster": i, "label": lab, "delta_next": 0.1, "n_cells": 50 + i}
        r.update({f"score_{c}": 0.5 for c in cfg.CT6})
        rows.append(r)
    return pd.DataFrame(rows)


def _write(tmp_path, df):
    p = tmp_path / "cluster_names.csv"
    df.to_csv(p, index=False)
    return p


def test_valid_file_passes(tmp_path):
    out = cn.read_cluster_names(_write(tmp_path, _df()))
    assert len(out) == 3 and set(out["label"]) <= set(cfg.CT6)


def test_label_outside_ct6_raises(tmp_path):
    df = _df()
    df.loc[0, "label"] = "Neutrophils"
    with pytest.raises(ValueError, match="Neutrophils"):
        cn.read_cluster_names(_write(tmp_path, df))


def test_duplicate_cluster_row_raises(tmp_path):
    df = _df()
    df.loc[1, "cluster"] = 0
    with pytest.raises(ValueError, match="duplicate"):
        cn.read_cluster_names(_write(tmp_path, df))


def test_missing_score_column_raises(tmp_path):
    df = _df().drop(columns=["score_Endothelial"])
    with pytest.raises(ValueError, match="score_Endothelial"):
        cn.read_cluster_names(_write(tmp_path, df))


@pytest.mark.parametrize("col", ["label", "delta_next", "score_B cells"])
def test_nan_raises(tmp_path, col):
    df = _df()
    df.loc[0, col] = float("nan")
    with pytest.raises(ValueError, match="NaN"):
        cn.read_cluster_names(_write(tmp_path, df))


@pytest.mark.parametrize("bad", [0, -3, 2.5])
def test_bad_n_cells_raises(tmp_path, bad):
    df = _df()
    df["n_cells"] = df["n_cells"].astype(float)
    df.loc[0, "n_cells"] = bad
    with pytest.raises(ValueError, match="n_cells"):
        cn.read_cluster_names(_write(tmp_path, df))
