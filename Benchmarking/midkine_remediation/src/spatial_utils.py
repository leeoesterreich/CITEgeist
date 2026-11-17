"""Lattice loading and spatial-field helpers shared by the G6/G7 synthetic gates.

Coordinates and QP proportions are the only real inputs (protocol C7: masks from QP
proportions and coordinates only; no SACE output exists before the freeze).
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import json
import os
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree

DATA_ROOT = Path(str(_cg.root("DATA") / "GSE289326" / "processed_files"))
REPO_ROOT = Path(os.environ.get("A12_REPO_ROOT") or Path(__file__).resolve().parents[3])
PROTOCOL_JSON = Path(__file__).resolve().parents[1] / "protocol.json"


def _protocol_lattice_dir() -> Path:
    """The folder the synthetic-gate lattice is recounted from: the protocol's inputs.synthetic_lattice_dir.
    (The old default output/module3_7type was wiped on 2026-10-22 and is not rebuilt.)"""
    return REPO_ROOT / json.loads(PROTOCOL_JSON.read_text())["inputs"]["synthetic_lattice_dir"]


QP_DIR = _protocol_lattice_dir()


def load_lattice(specimen: str, min_prop: float = 0.0, data_root: Path = DATA_ROOT, qp_dir: Path = QP_DIR) -> Tuple[np.ndarray, np.ndarray]:
    """(coords[n,2] in pixel units, cancer qp_prop[n]) for spots in the QP csv with cancer proportion >= min_prop."""
    qp = pd.read_csv(qp_dir / f"{specimen}_cell_prop_finetuned_results.csv", index_col=0)
    pos = pd.read_csv(data_root / specimen / "outs" / "spatial" / "tissue_positions.csv", header=None if _headerless(data_root / specimen / "outs" / "spatial" / "tissue_positions.csv") else 0, index_col=0)
    pos = pos.iloc[:, -2:]
    pos.columns = ["row", "col"]
    common = qp.index.intersection(pos.index)
    q = qp.loc[common, "Cancer_Cells"].to_numpy(float)
    xy = pos.loc[common, ["col", "row"]].to_numpy(float)
    keep = q >= min_prop
    return xy[keep], q[keep]


def _headerless(path: Path) -> bool:
    with open(path) as fh:
        first = fh.readline()
    return not first.lower().startswith("barcode")


def knn_weights(xy: np.ndarray, k: int = 6) -> sp.csr_matrix:
    """Row-standardised KNN weights (self excluded), the protocol's W."""
    n = len(xy)
    _, idx = cKDTree(xy).query(xy, k=k + 1)
    idx = idx[:, 1:]
    rows = np.repeat(np.arange(n), k)
    return sp.csr_matrix((np.full(n * k, 1.0 / k), (rows, idx.ravel())), shape=(n, n))


def gaussian_smoother(xy: np.ndarray, bw: float) -> sp.csr_matrix:
    """Sparse Gaussian-kernel smoother K with unit-variance output: field = K @ white_noise (bandwidth ``bw`` in coordinate units)."""
    tree = cKDTree(xy)
    D = tree.sparse_distance_matrix(tree, 3 * bw, output_type="coo_matrix")
    K = sp.csr_matrix((np.exp(-(D.data**2) / (2 * bw**2)), (D.row, D.col)), shape=(len(xy), len(xy)))
    K = K + sp.identity(len(xy), format="csr")           # sparse_distance_matrix drops the zero-distance diagonal
    norm = np.sqrt(np.asarray(K.multiply(K).sum(axis=1)).ravel())
    return sp.diags(1.0 / norm) @ K


def smooth_field(K: sp.csr_matrix, rng: np.random.Generator) -> np.ndarray:
    return K @ rng.standard_normal(K.shape[0])


def median_pitch(xy: np.ndarray) -> float:
    d, _ = cKDTree(xy).query(xy, k=2)
    return float(np.median(d[:, 1]))
