"""A10 Step 2: reconstruct Tangram proportions + per-type GEX layers from the saved mappings.

Ruling B3 (amendment-1b review): Tangram per-type GEX = Tangram's native mapping projection
(for type t: sum over mapped cells of type t of mapping weight x single-cell RAW expression;
no log, no per-spot weight normalisation), using the INNER 2026-04-28 mapping
(tangram_map_N/tangram_map_N/...), pinned by md5. Proportions = tg.project_cell_annotations
(mapping.T @ one-hot(cell types)), row-normalised, i.e. what Tangram/src/run_benchmark.py does.

Outputs (never overwrites existing predictions):
  regen/results/tangram/<cond>/Wu_rep_<r>/cell_type_proportions.csv
  regen/results/tangram/<cond>/Wu_rep_<r>/layers/<type>_layer.csv   (spots x GT-genes, native raw projection)
  regen/results/tangram/<cond>/Wu_rep_<r>/provenance.json
Diagnostic only (rep 0 of each condition): layers_legacy/ = run_benchmark.extract_gex_layers
(per-spot-normalised weights @ log-normalised reference), to explain any mismatch vs the cached values.
"""

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import argparse
import hashlib
import json
import os
import shutil
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData

try:
    from fast_matrix_market import mmread
except ImportError:
    from scipy.io import mmread

B = Path(str(_cg.REPO_ROOT / "Benchmarking" / "simulation_benchmarking"))
WU_REF = str(_cg.root("DATA") / "BC-Datasets" / "Wu_2021" / "Wu_etal_2021_BRCA_scRNASeq")

# Inner (2026-04-28) mapping md5s, pinned (measured 2026-10-14).
PINNED_MD5 = {
    ("high_seg", 0): "fe4451d5fb303ca3dca58723e0746dc5",
    ("high_seg", 1): "301f2be03523bfda1e0b2dbddaa7d3e9",
    ("high_seg", 2): "667e0c16606a2aaaf9721e3e8bfff71e",
    ("high_seg", 3): "f126f3a78bb68ea695349d447cbbfdff",
    ("high_seg", 4): "da53565b48855dd09b201779b1436509",
    ("mixed", 0): "96281b77f9ec7118cf589b0a17e86ea9",
    ("mixed", 1): "ac79f214fc1d69794354b9d0f8e0fd78",
    ("mixed", 2): "c4a09d2ebb530125047e7e47f4396077",
    ("mixed", 3): "295d44f99cb59f45d376df56b6bbee3b",
    ("mixed", 4): "4aa91b79fdbd57b512a1530f2eebf5f1",
}


# OUTER (2026-04-25) mapping md5s (diagnostic generation; reproduces the cached prop_comparison.json Tangram values)
PINNED_MD5_OUTER = {
    ("high_seg", 0): "0435d415e3b9b0e1055e53db0c51e1e0",
    ("high_seg", 1): "7135ad1c67e28ad95db22340cd9b4009",
    ("high_seg", 2): "6b80411361122cc8317305588caeb8de",
    ("high_seg", 3): "9e08bd8d4a31a890c9d4f76260f03c85",
    ("high_seg", 4): "20d78912ef2ae7c5158609e832818517",
    ("mixed", 0): "d8efc96cb5fd7ea000b7af303f1a2d2c",
    ("mixed", 1): "5783e845b35a9cd2762b7b0536b79c27",
    ("mixed", 2): "e1ff01da1a93932c6e699d999527880f",
    ("mixed", 3): "4506b05af9b2216b2e05a139928f41c0",
    ("mixed", 4): "98089ca5d828d26f31941e3e7af6569a",
}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def load_reference_raw():
    """Exactly Tangram/src/run_benchmark.load_wu_reference, minus normalize/log1p (raw counts kept)."""
    genes_df = pd.read_csv(os.path.join(WU_REF, "genes.tsv"), header=None, sep="\t")
    barcodes = pd.read_csv(os.path.join(WU_REF, "barcodes.tsv"), header=None)[0]
    matrix = mmread(os.path.join(WU_REF, "matrix.mtx")).T.tocsc().tocsr().astype("float32")
    adata = AnnData(matrix, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=genes_df[0]))
    adata.obs.index.name = None
    meta = pd.read_csv(os.path.join(WU_REF, "metadata.csv"), index_col=0)
    adata.obs["cell_type"] = meta["celltype_major"].values
    adata.obs["subtype"] = meta["subtype"].values
    adata = adata[adata.obs["subtype"] == "ER+"].copy()
    sc.pp.filter_cells(adata, min_genes=10)
    sc.pp.filter_genes(adata, min_cells=3)
    return adata


def write_csv_safe(df, path):
    """Write via node-local tmp then move (NFS wide-CSV corruption)."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, dir="/tmp") as tmp:
        tp = Path(tmp.name)
    df.to_csv(tp, float_format="%.7g")
    shutil.move(str(tp), str(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", default="0-9", help="array-style task ids 0-9 (t<5 high_seg, else mixed)")
    ap.add_argument("--out-root", default=str(B / "regen/results/tangram"))
    ap.add_argument("--gen", default="inner", choices=["inner", "outer"], help="mapping generation (inner = ruling B3)")
    ap.add_argument("--legacy-all", action="store_true", help="write layers_legacy for every rep (default: rep 0 only)")
    args = ap.parse_args()
    lo, hi = [int(x) for x in args.tasks.split("-")]

    t0 = time.time()
    ref = load_reference_raw()
    print(f"reference (ER+, raw counts) {ref.shape}  load {time.time()-t0:.0f}s", flush=True)
    ref_up = ref.var_names.str.upper()
    ref_logn = None  # built lazily for the diagnostic

    for t in range(lo, hi + 1):
        cond, rep = ("high_seg", t) if t < 5 else ("mixed", t - 5)
        ts = time.time()
        if args.gen == "inner":
            mp = B / f"Tangram/{cond}/tangram_map_{rep}/tangram_map_{rep}/Tangram_{cond}_rep_{rep}_mapping.h5ad"
            pin = PINNED_MD5[(cond, rep)]
        else:
            mp = B / f"Tangram/{cond}/tangram_map_{rep}/Tangram_{cond}_rep_{rep}_mapping.h5ad"
            pin = PINNED_MD5_OUTER[(cond, rep)]
        got = md5(mp)
        assert got == pin, f"mapping md5 mismatch {mp}: {got}"
        ad_map = sc.read_h5ad(mp)
        assert list(ad_map.obs_names) == list(ref.obs_names), "mapping cells != reference cells (order/ids)"
        assert (ad_map.obs["cell_type"].astype(str).values == ref.obs["cell_type"].astype(str).values).all()
        M = ad_map.X.toarray() if sp.issparse(ad_map.X) else np.asarray(ad_map.X)  # cells x spots
        spots = np.asarray(ad_map.var_names)
        types = ref.obs["cell_type"].astype(str).values
        type_list = sorted(set(types))
        assert len(type_list) == 9, type_list

        gt_dir = B / f"replicates/{cond}/ST_GEX_sim/sample_{rep}/layers"
        gt_cols = pd.read_csv(next(gt_dir.glob("*_GT.csv")), index_col=0, nrows=3)
        gt_idx = pd.read_csv(next(gt_dir.glob("*_GT.csv")), index_col=0, usecols=[0]).index
        assert set(spots) == set(gt_idx) and len(spots) == len(gt_idx), "mapping spots != GT spots"
        gt_genes = set(gt_cols.columns.str.upper())
        keep = np.array([g in gt_genes for g in ref_up])
        genes_out = ref_up[keep]

        out = Path(args.out_root) / cond / f"Wu_rep_{rep}"
        (out / "layers").mkdir(parents=True, exist_ok=True)

        # proportions: tg.project_cell_annotations = mapping.T @ one-hot(types), then row-normalise
        onehot = np.stack([(types == c).astype(np.float64) for c in type_list], axis=1)  # cells x types
        ct_pred = M.T @ onehot  # spots x types
        props = pd.DataFrame(ct_pred, index=spots, columns=type_list)
        props = props.div(props.sum(axis=1), axis=0).fillna(0)
        write_csv_safe(props, out / "cell_type_proportions.csv")

        # native raw projection layers
        Xk = ref.X[:, np.where(keep)[0]]
        Xk = Xk.tocsr() if sp.issparse(Xk) else sp.csr_matrix(Xk)
        for c in type_list:
            m = types == c
            L = (sp.csr_matrix(M[m]).T @ Xk[m]).toarray()  # spots x genes
            df = pd.DataFrame(L, index=spots, columns=genes_out)
            write_csv_safe(df, out / "layers" / f"{c.replace(' ', '_')}_layer.csv")

        diag = None
        if rep == 0 or args.legacy_all:  # diagnostic: legacy extract_gex_layers (normalised weights @ log-normalised reference)
            if ref_logn is None:
                tmp = ref.copy()
                sc.pp.normalize_total(tmp)
                sc.pp.log1p(tmp)
                ref_logn = tmp.X
            Xl = ref_logn[:, np.where(keep)[0]]
            Xl = Xl.toarray() if sp.issparse(Xl) else Xl
            (out / "layers_legacy").mkdir(exist_ok=True)
            for c in type_list:
                m = types == c
                W = M[m] / (M[m].sum(axis=0, keepdims=True) + 1e-10)
                df = pd.DataFrame(W.T @ Xl[m], index=spots, columns=genes_out)
                write_csv_safe(df, out / "layers_legacy" / f"{c.replace(' ', '_')}_layer.csv")
            diag = "layers_legacy written (diagnostic only)"

        prov = {
            "condition": cond, "replicate": rep, "mapping_path": str(mp), "mapping_md5": got,
            "n_cells": int(M.shape[0]), "n_spots": int(M.shape[1]), "n_types": len(type_list), "types": type_list,
            "n_layer_genes": int(len(genes_out)), "n_dup_upper_genes": int(pd.Index(genes_out).duplicated().sum()),
            "gex_method": "native mapping projection: sum_c in type t M[c,s]*raw_X[c,g]",
            "reference": WU_REF, "diagnostic": diag, "seconds": round(time.time() - ts, 1),
        }
        (out / "provenance.json").write_text(json.dumps(prov, indent=2))
        print(f"[{cond} rep{rep}] done in {time.time()-ts:.0f}s ({len(genes_out)} genes)", flush=True)
    print(f"ALL DONE {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
