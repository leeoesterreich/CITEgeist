"""Tests for Benchmarking/xenium_benchmarking/gt_rebuild/cluster_cells.py (synthetic data only)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("leidenalg")
h5py = pytest.importorskip("h5py")
pytest.importorskip("pyarrow")

import anndata as ad  # noqa: E402
from scipy.sparse import csc_matrix, csr_matrix  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking"))

from gt_rebuild import cluster_cells as cc  # noqa: E402

pytestmark = pytest.mark.unit

GENES = [f"G{i}" for i in range(60)]


def _blobs(n_per=60, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for b in range(3):
        mu = np.full(len(GENES), 0.3)
        lo = b * 10
        mu[lo : lo + 10] = 8.0
        rows.append(rng.poisson(mu, size=(n_per, len(GENES))))
    X = np.vstack(rows).astype(np.float32)
    ids = [f"c{i}" for i in range(X.shape[0])]
    a = ad.AnnData(csr_matrix(X))
    a.obs_names = ids
    a.var_names = GENES
    return a, np.repeat([0, 1, 2], n_per)


def _write_h5(path, X, barcodes, names, ftypes):
    """Write a 10x-style cell_feature_matrix.h5 (features x cells, CSC)."""
    m = csc_matrix(X.T)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=m.data)
        g.create_dataset("indices", data=m.indices)
        g.create_dataset("indptr", data=m.indptr)
        g.create_dataset("shape", data=np.array(m.shape))
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        ft = g.create_group("features")
        ft.create_dataset("name", data=np.array(names, dtype="S"))
        ft.create_dataset("feature_type", data=np.array(ftypes, dtype="S"))


def _write_cells(path, ids):
    pd.DataFrame({"cell_id": ids, "x_centroid": 0.0}).to_parquet(path)


def _fixture(tmp_path, X, ids, names, ftypes, cell_ids=None):
    h5 = tmp_path / "m.h5"
    pq = tmp_path / "cells.parquet"
    _write_h5(h5, X, ids, names, ftypes)
    _write_cells(pq, ids if cell_ids is None else cell_ids)
    return h5, pq


def test_three_blobs_give_three_clusters():
    a, truth = _blobs()
    lab = cc.cluster_level1(a, resolution=0.8, seed=0)
    assert list(lab.index) == list(a.obs_names)
    assert lab.nunique() == 3
    assert pd.crosstab(lab.values, truth).gt(0).sum(axis=1).eq(1).all()


def test_same_seed_identical_labels():
    a, _ = _blobs()
    assert cc.cluster_level1(a, 0.8, 0).equals(cc.cluster_level1(a, 0.8, 0))


def test_input_counts_not_mutated():
    a, _ = _blobs()
    before = a.X.toarray().copy()
    cc.cluster_level1(a, 0.8, 0)
    assert np.array_equal(before, a.X.toarray())


def test_zero_count_cells_dropped_and_counted(tmp_path):
    rng = np.random.default_rng(1)
    X = rng.poisson(2.0, size=(10, 4)).astype(np.float32) + 1
    X[3] = 0
    X[7] = 0
    ids = [f"c{i}" for i in range(10)]
    names = ["A", "B", "C", "D"]
    h5, pq = _fixture(tmp_path, X, ids, names, ["Gene Expression"] * 4)
    a, qc = cc.load_xenium_rna(h5, pq)
    assert a.n_obs == 8 and qc["n_zero_count_dropped"] == 2
    assert qc["dropped_ids"] == ["c3", "c7"]
    assert qc["n_cells_total"] == 10 and qc["n_genes"] == 4
    assert "c3" not in a.obs_names


def test_control_probe_columns_removed(tmp_path):
    rng = np.random.default_rng(2)
    names = ["A", "B", "NegControlProbe_1", "BLANK_01", "antisense_X", "UnassignedCodeword_1", "Intergenic_1"]
    X = rng.poisson(2.0, size=(6, len(names))).astype(np.float32) + 1
    ids = [f"c{i}" for i in range(6)]
    h5, pq = _fixture(tmp_path, X, ids, names, ["Gene Expression"] * len(names))
    a, qc = cc.load_xenium_rna(h5, pq)
    assert list(a.var_names) == ["A", "B"]
    assert qc["n_control_removed"] == 5 and qc["n_genes"] == 2


def test_non_gene_expression_features_excluded(tmp_path):
    rng = np.random.default_rng(3)
    X = rng.poisson(2.0, size=(5, 3)).astype(np.float32) + 1
    h5, pq = _fixture(tmp_path, X, [f"c{i}" for i in range(5)], ["A", "B", "C"], ["Gene Expression"] * 2 + ["Protein"])
    a, _ = cc.load_xenium_rna(h5, pq)
    assert list(a.var_names) == ["A", "B"]


def test_barcode_order_mismatch_raises(tmp_path):
    X = np.ones((4, 2), dtype=np.float32)
    ids = ["c0", "c1", "c2", "c3"]
    h5, pq = _fixture(tmp_path, X, ids, ["A", "B"], ["Gene Expression"] * 2, cell_ids=["c1", "c0", "c2", "c3"])
    with pytest.raises(ValueError):
        cc.load_xenium_rna(h5, pq)


def test_barcode_content_mismatch_raises(tmp_path):
    X = np.ones((4, 2), dtype=np.float32)
    ids = ["c0", "c1", "c2", "c3"]
    h5, pq = _fixture(tmp_path, X, ids, ["A", "B"], ["Gene Expression"] * 2, cell_ids=["c0", "c1", "c2", "zz"])
    with pytest.raises(ValueError):
        cc.load_xenium_rna(h5, pq)
