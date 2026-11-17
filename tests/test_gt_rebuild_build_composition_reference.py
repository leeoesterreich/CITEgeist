"""Tests for gt_rebuild/build_composition_reference.py (synthetic data only)."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")
import anndata as ad  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XB))

from gt_rebuild import build_composition_reference as bcr  # noqa: E402
from gt_rebuild import composition_arms as ca  # noqa: E402

CT = list(ca.CT6)
PER_TYPE = 12
SAMPLES = ["RCC1t", "RCC2t"]


@pytest.fixture
def adata():
    rng = np.random.default_rng(0)
    n = PER_TYPE * len(CT)
    labels = np.repeat(CT, PER_TYPE)
    obs = pd.DataFrame(
        {
            "cell_type": labels,
            "sample": [SAMPLES[i % 2] for i in range(n)],
            "batch": "b",
        },
        index=[f"cell{i}" for i in range(n)],
    )
    x = rng.poisson(2, size=(n, 7)).astype(np.float32)
    a = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(7)]))
    a.raw = a.copy()
    return a


def _keep(adata, counts):
    return ca.stratified_sample(adata.obs, counts, seed=0)


COUNTS = {"B cells": 2, "T cells": 3, "Macrophages": 4, "Fibroblasts": 1, "Epithelial": 5, "Endothelial": 6}


def test_counts_columns_equal_cell_types_index(adata, tmp_path):
    keep = _keep(adata, COUNTS)
    bcr.build_arm(adata, keep, tmp_path / "arm")
    for tool in ("rctd", "seurat"):
        d = tmp_path / "arm" / tool
        cols = pd.read_csv(d / "counts.csv", index_col=0, nrows=1).columns
        ct = pd.read_csv(d / "cell_types.csv", index_col=0)
        assert list(cols) == list(ct.index) == list(keep)
        assert not (d / "counts.csv").is_symlink()


def test_sample_column_present_for_card(adata, tmp_path):
    keep = _keep(adata, COUNTS)
    bcr.build_arm(adata, keep, tmp_path / "arm")
    ws = pd.read_csv(tmp_path / "arm" / "rctd" / "cell_types_with_sample.csv", index_col=0)
    assert list(ws.columns) == ["cell_type", "sample"]
    assert list(ws.index) == list(keep)
    assert set(ws["sample"]) <= set(SAMPLES)


def test_h5ad_exports_hold_only_kept_cells(adata, tmp_path):
    keep = _keep(adata, COUNTS)
    bcr.build_arm(adata, keep, tmp_path / "arm")
    for p in ("cell2location", "tangram"):
        a = ad.read_h5ad(tmp_path / "arm" / p / "reference.h5ad")
        assert list(a.obs_names) == list(keep)
    t = ad.read_h5ad(tmp_path / "arm" / "tangram" / "reference.h5ad")
    assert t.X.max() < 20  # log-normalised, not raw counts


def test_refuses_v40_processed_dirs():
    for d in ("processed", "processed_protein7", "processed_march", "processed_protein7_march"):
        with pytest.raises(ValueError):
            bcr.assert_not_v40(bcr.BASE / d / "x")
    bcr.assert_not_v40(bcr.BASE / "processed_comp" / "natural_s0")


def test_refuses_protected_dir(tmp_path):
    lst = tmp_path / "protected.txt"
    lst.write_text(str(tmp_path / "prot") + "\n")
    with pytest.raises(ValueError):
        bcr.assert_not_protected(tmp_path / "prot" / "sub", protected_file=lst)
    bcr.assert_not_protected(tmp_path / "other", protected_file=lst)


def _manifest(tmp_path, ref_csv, counts, sha=None):
    m = {
        "arm_set": ["natural"],
        "seeds": [0],
        "integer_counts": {"natural": counts},
        "reference_csv": str(ref_csv),
        "reference_sha256": sha or bcr.manifest.sha256_file(ref_csv),
    }
    p = tmp_path / "arms_manifest.json"
    p.write_text(json.dumps(m))
    return p


def test_label_counts_equal_manifest_counts(adata, tmp_path):
    keep = _keep(adata, COUNTS)
    bcr.build_arm(adata, keep, tmp_path / "arm")
    bcr.check_counts(adata.obs.loc[keep, "cell_type"], COUNTS)
    bad = dict(COUNTS, **{"B cells": 3})
    with pytest.raises(ValueError):
        bcr.check_counts(adata.obs.loc[keep, "cell_type"], bad)


def test_reference_sha_mismatch_refused(adata, tmp_path):
    ref = tmp_path / "ref.csv"
    adata.obs[["cell_type", "sample"]].to_csv(ref)
    good = bcr.load_manifest(_manifest(tmp_path, ref, COUNTS))
    bcr.verify_reference(good)
    ref.write_text(ref.read_text() + "\n")
    with pytest.raises(ValueError):
        bcr.verify_reference(good)


def test_export_refuses_without_signoff(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise bcr.manifest.SignoffMissing("no signoff")

    monkeypatch.setattr(bcr.manifest, "require_rnagt_signoff", boom)
    monkeypatch.setattr(bcr, "load_manifest", lambda *a, **k: pytest.fail("touched data before sign-off"))
    with pytest.raises(bcr.manifest.SignoffMissing):
        bcr.main(["--arm", "natural", "--seed", "0"])
    assert list(tmp_path.iterdir()) == []


def test_index_maps_to_arm_seed_from_manifest(tmp_path):
    m = {"arm_set": ["a", "b"], "seeds": [0, 1, 2]}
    pairs = [bcr.arm_seed_for_index(m, i) for i in range(6)]
    assert pairs == [("a", 0), ("a", 1), ("a", 2), ("b", 0), ("b", 1), ("b", 2)]
    with pytest.raises(IndexError):
        bcr.arm_seed_for_index(m, 6)


def test_report_written_with_counts_and_hashes(adata, tmp_path):
    keep = _keep(adata, COUNTS)
    rep = bcr.write_build_report(
        tmp_path / "build_report.json", adata.obs, keep, COUNTS, "natural", 0, {"reference_csv": "abc"}
    )
    on_disk = json.loads((tmp_path / "build_report.json").read_text())
    assert on_disk["label_counts"] == COUNTS == rep["label_counts"]
    assert on_disk["n_cells"] == sum(COUNTS.values())
    assert sum(on_disk["per_sample_counts"].values()) == sum(COUNTS.values())
    assert len(on_disk["kept_ids_sha256"]) == 64
    assert on_disk["inputs_sha256"] == {"reference_csv": "abc"}
    assert "_code_commit" in on_disk


# ------------------------------------------------------------------ full_size endpoint (Task 34)
import hashlib  # noqa: E402

pytestmark = pytest.mark.unit

FULL_COUNTS = {t: PER_TYPE for t in CT}


def test_endpoint_keeps_every_cell_in_obs_order(adata):
    keep = bcr.endpoint_keep(adata.obs, FULL_COUNTS)
    assert list(keep) == list(adata.obs.index)
    with pytest.raises(ValueError):
        bcr.endpoint_keep(adata.obs, dict(FULL_COUNTS, **{"B cells": PER_TYPE - 1}))


def _endpoint_manifest(tmp_path, adata):
    ref = tmp_path / "ref.csv"
    adata.obs[["cell_type", "sample"]].to_csv(ref)
    p = tmp_path / "full_size_manifest.json"
    p.write_text(
        json.dumps(
            {
                "label_counts": FULL_COUNTS,
                "reference_csv": str(ref),
                "reference_sha256": bcr.manifest.sha256_file(ref),
            }
        )
    )
    return p


def test_endpoint_refuses_without_signoff_and_writes_nothing(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise bcr.manifest.SignoffMissing("no signoff")

    monkeypatch.setattr(bcr.manifest, "require_rnagt_signoff", boom)
    monkeypatch.setattr(bcr, "load_manifest", lambda *a, **k: pytest.fail("touched data before sign-off"))
    monkeypatch.setenv("XEN_BASE", str(tmp_path))
    with pytest.raises(bcr.manifest.SignoffMissing):
        bcr.main(["--endpoint", "full_size"])
    assert list(tmp_path.iterdir()) == []


def test_endpoint_refuses_nonempty_out_root(adata, tmp_path, monkeypatch):
    monkeypatch.setattr(bcr.manifest, "require_rnagt_signoff", lambda **k: {})
    monkeypatch.setenv("XEN_BASE", str(tmp_path / "xb"))
    out = bcr.bp.comp_ref_dir("rctd", "full_size", 0).parent
    out.mkdir(parents=True)
    (out / "x").write_text("x")
    man = _endpoint_manifest(tmp_path, adata)
    with pytest.raises(FileExistsError):
        bcr.main(["--endpoint", "full_size", "--manifest", str(man)])


def test_endpoint_out_root_not_under_processed_comp_or_processed():
    out = bcr.bp.comp_ref_dir("rctd", "full_size", 0).parent
    for forbidden in (bcr.BASE / "processed", bcr.BASE / "processed_comp"):
        assert forbidden not in out.parents and out != forbidden
    bcr.assert_not_v40(out)


def test_endpoint_flags_mutually_exclusive():
    for extra in (["--arm", "natural"], ["--seed", "0"], ["--index", "0"]):
        with pytest.raises(SystemExit):
            bcr.main(["--endpoint", "full_size", *extra])
    with pytest.raises(SystemExit):
        bcr.main(["--verify-only"])


def test_build_report_accepts_endpoint_without_seed(adata, tmp_path):
    rep = bcr.write_build_report(
        tmp_path / "r.json", adata.obs, adata.obs.index, FULL_COUNTS, "full_size", None, {"x": "y"}
    )
    assert json.loads((tmp_path / "r.json").read_text())["seed"] is None
    assert rep["n_cells"] == len(adata.obs)


def _twin(adata, tmp_path):
    """Two identical full exports (stand-ins for out_root and processed/)."""
    a, b = tmp_path / "new", tmp_path / "proc"
    for d in (a, b):
        bcr.build_arm(adata, adata.obs.index, d)
    return a, b


def test_verify_identical_exports_passes(adata, tmp_path):
    a, b = _twin(adata, tmp_path)
    rep = bcr.verify_matches_processed(a, b, chunk_rows=5)
    assert set(rep) == set(bcr.H5AD_FILES + bcr.COUNTS_FILES + bcr.TABLE_FILES)
    assert rep["rctd/counts.csv"]["byte_equal"] is True


def _edit_last_row_csv(path):
    lines = path.read_text().splitlines()
    cells = lines[-1].split(",")
    cells[-1] = str(float(cells[-1]) + 1)
    lines[-1] = ",".join(cells)
    path.write_text("\n".join(lines) + "\n")


def _edit_last_row_h5ad(path):
    a = ad.read_h5ad(path)
    x = a.X.toarray() if hasattr(a.X, "toarray") else np.array(a.X)
    x[-1, -1] += 1
    a.X = x
    a.write_h5ad(path)


@pytest.mark.parametrize(
    "mutate",
    [
        lambda a: _edit_last_row_csv(a / "rctd" / "counts.csv"),
        lambda a: _edit_last_row_csv(a / "seurat" / "counts.csv"),
        lambda a: _edit_last_row_h5ad(a / "cell2location" / "reference.h5ad"),
        lambda a: _edit_last_row_h5ad(a / "tangram" / "reference.h5ad"),
    ],
    ids=["rctd_counts", "seurat_counts", "c2l_h5ad", "tangram_h5ad"],
)
def test_verify_detects_one_changed_count_in_last_row(adata, tmp_path, mutate):
    a, b = _twin(adata, tmp_path)
    mutate(a)
    with pytest.raises(bcr.EquivalenceError):
        bcr.verify_matches_processed(a, b, chunk_rows=5)


def test_verify_detects_swapped_obs_and_var_names_sample_raw_and_layer(adata, tmp_path):
    def swap(path, which):
        x = ad.read_h5ad(path)
        if which == "obs":
            names = list(x.obs_names)
            names[0], names[1] = names[1], names[0]
            x.obs_names = names
        elif which == "var":
            names = list(x.var_names)
            names[0], names[1] = names[1], names[0]
            x.var_names = names
        elif which == "sample":
            x.obs["sample"] = x.obs["sample"].astype(str)
            x.obs.iloc[0, x.obs.columns.get_loc("sample")] = "RCC9t"
        elif which == "raw":
            x.raw = x.copy()
        elif which == "layer":
            x.layers["extra"] = x.X.copy()
        x.write_h5ad(path)

    for which in ("obs", "var", "sample", "raw", "layer"):
        sub = tmp_path / which
        sub.mkdir()
        a, b = _twin(adata, sub)
        swap(a / "cell2location" / "reference.h5ad", which)
        with pytest.raises(bcr.EquivalenceError):
            bcr.verify_matches_processed(a, b, chunk_rows=5)


def test_verify_detects_cell_types_change(adata, tmp_path):
    a, b = _twin(adata, tmp_path)
    p = a / "rctd" / "cell_types_with_sample.csv"
    p.write_text(p.read_text().replace("RCC1t", "RCC9t", 1))
    with pytest.raises(bcr.EquivalenceError):
        bcr.verify_matches_processed(a, b)


def test_hash_is_chunk_size_independent_and_dense_sparse_equal(adata, tmp_path):
    import scipy.sparse as sp

    x = np.asarray(adata.X)
    assert bcr.hash_matrix(x, 3) == bcr.hash_matrix(x, 50) == bcr.hash_matrix(sp.csr_matrix(x), 7)
    y = x.copy()
    y[-1, -1] += 1
    assert bcr.hash_matrix(y, 3) != bcr.hash_matrix(x, 3)
    p = tmp_path / "c.csv"
    pd.DataFrame(x.T, index=[f"g{i}" for i in range(x.shape[1])], columns=adata.obs_names).to_csv(p)
    assert bcr.hash_counts_csv(p, 2) == bcr.hash_counts_csv(p, 100)


def test_byte_inequality_with_content_equality_passes(adata, tmp_path):
    a, b = _twin(adata, tmp_path)
    df = pd.read_csv(a / "rctd" / "counts.csv", index_col=0)
    df.to_csv(a / "rctd" / "counts.csv", float_format="%.3f")  # 3.000 vs 3.0: same values, different bytes
    assert (
        hashlib.sha256((a / "rctd/counts.csv").read_bytes()).digest()
        != hashlib.sha256((b / "rctd/counts.csv").read_bytes()).digest()
    )
    rep = bcr.verify_matches_processed(a, b, chunk_rows=5)
    assert rep["rctd/counts.csv"]["byte_equal"] is False
