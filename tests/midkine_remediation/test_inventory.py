"""A12.2/A12.3: spot-UMI reader used for the pooled depth reference M (synthetic 10x-style h5)."""
import h5py
import numpy as np
import pytest
import scipy.sparse as sp

import inventory_depth_nuclei as INV

pytestmark = pytest.mark.unit


def _write(path, dense, barcodes):
    csc = sp.csc_matrix(dense)          # genes x barcodes, CSC by barcode like 10x
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g["barcodes"] = np.array([b.encode() for b in barcodes])
        g["data"], g["indices"], g["indptr"] = csc.data, csc.indices, csc.indptr
        g["shape"] = np.array(csc.shape)


def test_spot_umi_column_sums_including_empty_barcodes(tmp_path):
    dense = np.array([[1, 0, 0, 5], [2, 0, 0, 0], [0, 0, 0, 3]], float)   # barcode b (col 1) and c (col 2) are empty
    _write(tmp_path / "m.h5", dense, ["a", "b", "c", "d"])
    u = INV.spot_umi(tmp_path / "m.h5")
    assert u.to_dict() == {"a": 3.0, "b": 0.0, "c": 0.0, "d": 8.0}


def test_specimen_list_matches_protocol():
    import json
    from pathlib import Path

    proto = json.loads((Path(__file__).resolve().parents[2] / "Benchmarking/midkine_remediation/protocol.DRAFT.json").read_text())
    assert INV.SPECIMENS == proto["specimens"]["hurdle_specimens_all"]
