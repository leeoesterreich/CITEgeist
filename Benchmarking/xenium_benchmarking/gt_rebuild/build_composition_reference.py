"""Per-arm reference exports for the composition experiment (Task 23).

Subsets the May reference (the one the full-size competitors were scored with) to the
cells drawn for one (arm, seed) and writes the C2L / Tangram / RCTD / Seurat inputs
under ``bp.comp_ref_dir``. Counts are physically subset (the R scripts read the whole
counts.csv against cell_types.csv). Modelled on
``reference_data/GSE156632/src/build_march_reference.py``.
"""

import argparse
import hashlib
import json
import logging
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, Iterator, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
sys.path.insert(0, str(XB))
sys.path.insert(0, str(REPO / "scripts"))

import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402

from . import composition_arms as ca  # noqa: E402
from . import manifest  # noqa: E402

logger = logging.getLogger(__name__)

BASE = XB / "reference_data" / "GSE156632"
SRC = BASE / "src"
COMBINED_H5AD = BASE / "processed" / "GSE156632_combined.h5ad"
V40_DIRS = ("processed", "processed_protein7", "processed_march", "processed_protein7_march")
PRODUCED_BY = "Benchmarking/xenium_benchmarking/gt_rebuild/build_composition_reference.py"
PROTECTED_LIST = Path(
    os.environ.get("CITEGEIST_PROTECTED_PATHS_FILE", str(Path.home() / ".config/citegeist-dvc/protected_paths.txt"))
)


def assert_not_v40(path, base: Path = BASE) -> None:
    """Raise if path resolves into (or is) a v40 / March reference directory."""
    p = Path(path).resolve()
    for d in V40_DIRS:
        root = (Path(base) / d).resolve()
        if p == root or root in p.parents:
            raise ValueError(f"refusing to write into v40/March reference dir: {p}")


def assert_not_protected(path, protected_file: Optional[Path] = None) -> None:
    """Raise if path equals, contains or is inside a DVC-protected folder (mirrors _refuse_protected_outdir.sh)."""
    lst = Path(protected_file) if protected_file is not None else PROTECTED_LIST
    if not lst.is_file():
        logger.warning("no protected-path list at %s; guard skipped", lst)
        return
    t = Path(os.path.realpath(path))
    for line in lst.read_text().splitlines():
        if not line.strip():
            continue
        p = Path(os.path.realpath(line.strip()))
        if t == p or p in t.parents or t in p.parents:
            raise ValueError(f"{path} overlaps DVC-protected folder {line.strip()}")


def load_manifest(path: Path = ca.DEFAULT_MANIFEST) -> dict:
    return json.loads(Path(path).read_text())


def verify_reference(m: dict) -> str:
    """The reference labels file must be byte-identical to the one the arms were frozen against."""
    got = manifest.sha256_file(m["reference_csv"])
    if got != m["reference_sha256"]:
        raise ValueError(f"reference sha256 {got} != manifest {m['reference_sha256']}; refusing to sample")
    return got


def arm_seed_for_index(m: dict, index: int) -> Tuple[str, int]:
    pairs = [(a, s) for a in m["arm_set"] for s in m["seeds"]]
    if not 0 <= index < len(pairs):
        raise IndexError(f"array index {index} outside 0..{len(pairs) - 1}")
    return pairs[index]


def check_counts(cell_types: pd.Series, expected: Dict[str, int]) -> Dict[str, int]:
    got = {t: int((cell_types == t).sum()) for t in ca.CT6}
    want = {t: int(expected.get(t, 0)) for t in ca.CT6}
    if got != want or len(cell_types) != sum(want.values()):
        raise ValueError(f"kept label counts {got} != manifest counts {want}")
    return got


def build_arm(adata_combined, keep: pd.Index, out_dir: Path) -> None:
    """Write all tool exports for the kept cells into out_dir (counts physically subset)."""
    out_dir = Path(out_dir)
    assert_not_v40(out_dir)
    assert_not_protected(out_dir)
    sys.path.insert(0, str(SRC))
    import split_tcells_to_protein7 as split
    from process_reference import export_for_cell2location, export_for_tangram

    sub = adata_combined[list(keep)].copy()
    out_dir.mkdir(parents=True, exist_ok=True)
    export_for_cell2location(sub, out_dir)
    export_for_tangram(sub, out_dir)
    src = sub.raw.to_adata() if sub.raw is not None else sub
    x = src.X.toarray() if hasattr(src.X, "toarray") else src.X
    rctd = out_dir / "rctd"
    seurat = out_dir / "seurat"
    rctd.mkdir(exist_ok=True)
    seurat.mkdir(exist_ok=True)
    pd.DataFrame(x.T, index=src.var_names, columns=sub.obs_names).to_csv(rctd / "counts.csv")
    shutil.copyfile(rctd / "counts.csv", seurat / "counts.csv")
    for d in (rctd, seurat):
        sub.obs[["cell_type"]].to_csv(d / "cell_types.csv")
    split.write_cell_types_with_sample(sub, rctd)


def write_build_report(
    path: Path, obs: pd.DataFrame, keep: pd.Index, expected: Dict[str, int], arm: str, seed: Optional[int], inputs: dict
) -> dict:
    kept = obs.loc[list(keep)]
    payload = {
        "arm": arm,
        "seed": None if seed is None else int(seed),
        "n_cells": int(len(kept)),
        "label_counts": check_counts(kept["cell_type"], expected),
        "kept_ids_sha256": hashlib.sha256("\n".join(map(str, keep)).encode()).hexdigest(),
        "per_sample_counts": {str(k): int(v) for k, v in kept["sample"].value_counts().sort_index().items()},
        "per_sample_label_counts": {
            str(s): {t: int(v) for t, v in g["cell_type"].value_counts().items()} for s, g in kept.groupby("sample")
        },
        "inputs_sha256": inputs,
    }
    return provenance.write_stamped_json(path, payload, PRODUCED_BY)


class EquivalenceError(ValueError):
    """A full_size export differs in content from processed/."""


def _h(*parts: bytes) -> "hashlib._Hash":
    h = hashlib.sha256()
    for p in parts:
        h.update(p)
    return h


def hash_matrix(X, chunk_rows: int = 2000) -> str:
    """Chunk-size-independent, representation-independent (dense == CSR) sha256 over a whole matrix.

    Three streams (row lengths, column indices, values as float64) are hashed separately so the digest does
    not depend on where chunk boundaries fall.
    """
    import scipy.sparse as sp

    hl, hi, hd = hashlib.sha256(), hashlib.sha256(), hashlib.sha256()
    n = X.shape[0]
    for a in range(0, n, chunk_rows):
        blk = X[a : a + chunk_rows]
        blk = sp.csr_matrix(blk.toarray() if hasattr(blk, "toarray") and not sp.issparse(blk) else blk)
        blk.sum_duplicates()
        blk.eliminate_zeros()
        blk.sort_indices()
        hl.update(np.diff(blk.indptr).astype(np.int64).tobytes())
        hi.update(blk.indices.astype(np.int64).tobytes())
        hd.update(blk.data.astype(np.float64).tobytes())
    return _h(str(X.shape).encode(), hl.digest(), hi.digest(), hd.digest()).hexdigest()


def _csv_chunks(path: Path, chunk_rows: int) -> Iterator[pd.DataFrame]:
    return pd.read_csv(path, index_col=0, chunksize=chunk_rows)


def hash_counts_csv(path: Path, chunk_rows: int = 500) -> Tuple[list, str, str]:
    """(ordered header, sha256 of ordered row ids, sha256 of all parsed float64 values in file order)."""
    header = list(pd.read_csv(path, index_col=0, nrows=0).columns)
    hi, hv = hashlib.sha256(), hashlib.sha256()
    for chunk in _csv_chunks(path, chunk_rows):
        hi.update("\n".join(map(str, chunk.index)).encode() + b"\n")
        hv.update(np.ascontiguousarray(chunk.to_numpy(dtype=np.float64)).tobytes())
    return header, hi.hexdigest(), hv.hexdigest()


def sha256_stream(path: Path) -> str:
    return manifest.sha256_file(path)


def _fail(what: str, a, b) -> None:
    raise EquivalenceError(f"{what}: {a!r} != {b!r}")


def compare_h5ad(new: Path, ref: Path, chunk_rows: int = 2000) -> dict:
    import anndata as ad

    a, b = ad.read_h5ad(new), ad.read_h5ad(ref)
    if list(a.obs_names) != list(b.obs_names):
        _fail(f"{new.name} obs_names order", len(a.obs_names), len(b.obs_names))
    if list(a.var_names) != list(b.var_names):
        _fail(f"{new.name} var_names order", len(a.var_names), len(b.var_names))
    if list(a.obs.columns) != list(b.obs.columns):
        _fail(f"{new.name} obs columns", list(a.obs.columns), list(b.obs.columns))
    for c in a.obs.columns:
        if not a.obs[c].astype(object).equals(b.obs[c].astype(object)):
            _fail(f"{new.name} obs column {c}", "differs", "reference")
    if a.X.shape != b.X.shape:
        _fail(f"{new.name} X.shape", a.X.shape, b.X.shape)
    if (a.raw is None) != (b.raw is None):
        _fail(f"{new.name} raw presence", a.raw is not None, b.raw is not None)
    if sorted(a.layers.keys()) != sorted(b.layers.keys()):
        _fail(f"{new.name} layers", sorted(a.layers.keys()), sorted(b.layers.keys()))
    pairs = [("X", a.X, b.X)]
    if a.raw is not None:
        if a.raw.X.shape != b.raw.X.shape:
            _fail(f"{new.name} raw.X.shape", a.raw.X.shape, b.raw.X.shape)
        pairs.append(("raw.X", a.raw.X, b.raw.X))
    pairs += [(f"layers[{k}]", a.layers[k], b.layers[k]) for k in sorted(a.layers.keys())]
    hashes = {}
    for name, x, y in pairs:
        hx, hy = hash_matrix(x, chunk_rows), hash_matrix(y, chunk_rows)
        if hx != hy:
            _fail(f"{new.name} {name} matrix hash", hx, hy)
        hashes[name] = hx
    return {"matrix_hashes": hashes, "byte_equal": sha256_stream(new) == sha256_stream(ref)}


def compare_counts_csv(new: Path, ref: Path, chunk_rows: int = 500) -> dict:
    ha, ia, va = hash_counts_csv(new, chunk_rows)
    hb, ib, vb = hash_counts_csv(ref, chunk_rows)
    if ha != hb:
        _fail(f"{new.parent.name}/{new.name} header (cell ids)", len(ha), len(hb))
    if ia != ib:
        _fail(f"{new.parent.name}/{new.name} first column (gene ids)", ia, ib)
    if va != vb:
        _fail(f"{new.parent.name}/{new.name} values hash", va, vb)
    return {"values_hash": va, "byte_equal": sha256_stream(new) == sha256_stream(ref)}


def compare_table_csv(new: Path, ref: Path) -> dict:
    a, b = pd.read_csv(new, index_col=0), pd.read_csv(ref, index_col=0)
    if list(a.index) != list(b.index):
        _fail(f"{new.parent.name}/{new.name} index order", len(a), len(b))
    if list(a.columns) != list(b.columns) or not a.astype(object).equals(b.astype(object)):
        _fail(f"{new.parent.name}/{new.name} content", "differs", "reference")
    return {"byte_equal": sha256_stream(new) == sha256_stream(ref)}


H5AD_FILES = ("cell2location/reference.h5ad", "tangram/reference.h5ad")
COUNTS_FILES = ("rctd/counts.csv", "seurat/counts.csv")
TABLE_FILES = ("rctd/cell_types.csv", "rctd/cell_types_with_sample.csv", "seurat/cell_types.csv")


def verify_matches_processed(out_root: Path, processed_root: Path, chunk_rows: int = 2000) -> dict:
    """Semantic, full-content equivalence of every export a method reads. Byte equality is reported only.

    Raises EquivalenceError naming the first differing item.
    """
    out_root, processed_root = Path(out_root), Path(processed_root)
    report = {}
    for rel in H5AD_FILES:
        report[rel] = compare_h5ad(out_root / rel, processed_root / rel, chunk_rows)
    for rel in COUNTS_FILES:
        report[rel] = compare_counts_csv(out_root / rel, processed_root / rel)
    for rel in TABLE_FILES:
        report[rel] = compare_table_csv(out_root / rel, processed_root / rel)
    return report


def endpoint_keep(obs: pd.DataFrame, expected: Dict[str, int]) -> pd.Index:
    """full_size keeps EVERY reference cell, in reference-label-file order, no sampling."""
    check_counts(obs["cell_type"], expected)
    return obs.index


def _main_endpoint(a) -> None:
    fm = load_manifest(a.manifest or ca.FULL_SIZE_MANIFEST)
    out_root = bp.comp_ref_dir("rctd", "full_size", 0).parent
    assert_not_v40(out_root)
    assert_not_protected(out_root)
    if a.verify_only:
        rep = verify_matches_processed(out_root, BASE / "processed")
        (out_root / "equivalence_report.json").write_text(json.dumps(rep, indent=2))
        print(json.dumps(rep, indent=2))
        return
    if out_root.exists() and any(out_root.iterdir()):
        raise FileExistsError(f"REFUSED: {out_root} is not empty; builds write a new dir")
    ref_sha = verify_reference(fm)

    import scanpy as sc

    obs = pd.read_csv(fm["reference_csv"], index_col=0)
    obs["cell_type"], unmapped = ca.map_reference_labels(obs["cell_type"])
    if unmapped:
        raise ValueError(f"unmapped reference labels: {unmapped}")
    keep = endpoint_keep(obs, fm["label_counts"])
    adata = sc.read_h5ad(COMBINED_H5AD)
    if set(adata.obs_names) != set(obs.index):
        raise ValueError("combined h5ad barcodes differ from the frozen reference labels file")
    adata.obs["cell_type"] = obs["cell_type"].reindex(adata.obs_names).values
    adata.obs["sample"] = obs["sample"].reindex(adata.obs_names).values
    build_arm(adata, keep, out_root)
    write_build_report(
        out_root / "build_report.json",
        obs,
        keep,
        fm["label_counts"],
        "full_size",
        None,
        {
            "reference_csv": ref_sha,
            "combined_h5ad": manifest.sha256_file(COMBINED_H5AD),
            "full_size_manifest": manifest.sha256_file(a.manifest or ca.FULL_SIZE_MANIFEST),
        },
    )
    logger.info("built full_size: %d cells -> %s", len(keep), out_root)


def main(argv: Optional[Sequence[str]] = None) -> None:
    ap = argparse.ArgumentParser(description="Build one (arm, seed) reference export.")
    ap.add_argument("--arm")
    ap.add_argument("--seed", type=int)
    ap.add_argument("--index", type=int, help="array index into arm_set x seeds (manifest order)")
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--endpoint", choices=list(ca.ENDPOINTS), help="build the all-cells endpoint instead of an arm")
    ap.add_argument("--verify-only", action="store_true", help="with --endpoint: compare the export to processed/")
    a = ap.parse_args(argv)
    if a.endpoint and (a.arm is not None or a.seed is not None or a.index is not None):
        ap.error("--endpoint is mutually exclusive with --arm/--seed/--index")
    if a.verify_only and not a.endpoint:
        ap.error("--verify-only needs --endpoint")
    manifest.require_rnagt_signoff(labelset="rnaprot")
    if a.endpoint:
        return _main_endpoint(a)
    a.manifest = a.manifest or ca.DEFAULT_MANIFEST
    m = load_manifest(a.manifest)
    if a.index is not None:
        arm, seed = arm_seed_for_index(m, a.index)
    elif a.arm is not None and a.seed is not None:
        arm, seed = a.arm, a.seed
    else:
        ap.error("give --index or both --arm and --seed")
    if arm not in m["arm_set"] or seed not in m["seeds"]:
        raise ValueError(f"({arm}, {seed}) not in the frozen manifest")
    out_root = bp.comp_ref_dir("rctd", arm, seed).parent
    assert_not_v40(out_root)
    assert_not_protected(out_root)
    if out_root.exists() and any(out_root.iterdir()):
        raise FileExistsError(f"REFUSED: {out_root} is not empty; builds write a new dir")
    ref_sha = verify_reference(m)

    import scanpy as sc

    obs = pd.read_csv(m["reference_csv"], index_col=0)
    obs["cell_type"], unmapped = ca.map_reference_labels(obs["cell_type"])
    if unmapped:
        raise ValueError(f"unmapped reference labels: {unmapped}")
    keep = ca.stratified_sample(obs, m["integer_counts"][arm], seed)
    check_counts(obs.loc[keep, "cell_type"], m["integer_counts"][arm])
    adata = sc.read_h5ad(COMBINED_H5AD)
    if set(adata.obs_names) != set(obs.index):
        raise ValueError("combined h5ad barcodes differ from the frozen reference labels file")
    adata.obs["cell_type"] = obs["cell_type"].reindex(adata.obs_names).values
    adata.obs["sample"] = obs["sample"].reindex(adata.obs_names).values
    build_arm(adata, keep, out_root)
    write_build_report(
        out_root / "build_report.json",
        obs,
        keep,
        m["integer_counts"][arm],
        arm,
        seed,
        {
            "reference_csv": ref_sha,
            "combined_h5ad": manifest.sha256_file(COMBINED_H5AD),
            "arms_manifest": manifest.sha256_file(a.manifest),
        },
    )
    logger.info("built %s seed %s: %d cells -> %s", arm, seed, len(keep), out_root)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
