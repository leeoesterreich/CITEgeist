"""A12.6 probe (outcome-blind): per specimen, are the six basal keratins and MDK present in the raw GEX, and which survive filter_gex?

Reads only the Visium h5 (counts), the QP csv barcodes and the module3 filter rule (nonzero >= 1% of spots AND mean over nonzero spots > 1.0,
computed on the spots kept by the spot filter min_counts 100 (-S1) / 25). Writes a JSON; no outcome data are touched.
"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

GENES = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17", "MDK"]
ROOT = Path(str(_cg.root("DATA") / "GSE289326" / "processed_files"))
SP = ["HCC22-088-P1-S1", "HCC22-088-P1-S2", "HCC22-088-P2-S1", "HCC22-088-P2-S2", "HCC22-088-P3-S1_A", "HCC22-088-P3-S2", "HCC22-088-P4-S1",
      "HCC22-088-P4-S2_1i_rep", "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep", "HCC22-088-P6-S1", "HCC22-088-P6-S2_D"]
out = {}
for s in SP:
    with h5py.File(ROOT / s / "outs" / "filtered_feature_bc_matrix.h5", "r") as f:
        m = f["matrix"]
        bc = np.array([b.decode() for b in m["barcodes"][:]])
        ft = m["features"]["feature_type"][:]
        names = np.array([b.decode() for b in m["features"]["name"][:]])
        shape = m["shape"][:]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(shape)).tocsr()   # genes x spots
    ge = np.where(ft == b"Gene Expression")[0]
    X, names = X[ge], names[ge]
    qp = pd.read_csv(f"output/module3_7type/{s}_cell_prop_finetuned_results.csv", index_col=0)
    keep = np.isin(bc, qp.index)
    X = X[:, keep]
    nz = np.asarray((X > 0).sum(axis=1)).ravel()
    tot = np.asarray(X.sum(axis=1)).ravel()
    n = X.shape[1]
    rec = {"n_spots": int(n), "n_genes_all": int(len(names)), "total_counts_all_genes": float(tot.sum())}
    mean_nz = np.where(nz > 0, tot / np.maximum(nz, 1), 0)
    passes = (nz >= 0.01 * n) & (mean_nz > 1.0)
    rec["n_genes_pass_module3_filter"] = int(passes.sum())
    rec["total_counts_pass_module3_filter"] = float(tot[passes].sum())
    rec["n_genes_nonzero_ge_1pct"] = int((nz >= 0.01 * n).sum())
    for g in GENES:
        i = np.where(names == g)[0]
        if len(i) != 1:
            rec[g] = {"present": False, "n_matches": int(len(i))}
            continue
        i = int(i[0])
        rec[g] = {"present": True, "frac_spots_nonzero": float(nz[i] / n), "mean_over_nonzero": float(mean_nz[i]), "total": float(tot[i]),
                  "pass_nonzero_1pct": bool(nz[i] >= 0.01 * n), "pass_module3_filter": bool(passes[i])}
    out[s] = rec
    print(s, {g: (rec[g].get("pass_module3_filter"), round(rec[g].get("frac_spots_nonzero", 0), 3)) for g in GENES}, flush=True)
Path(sys.argv[1]).write_text(json.dumps(out, indent=1))
