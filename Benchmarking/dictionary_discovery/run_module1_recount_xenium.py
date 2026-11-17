"""Module 1 (marker interest) retained-marker count for one Xenium RCC region (v42 recount, Xenium analogue of run_module1_recount.py).

Same Module 1 call and settings as the breast recount (morans_k 8, 99 permutations, seed 42, raw protein counts).
One SLURM array task = one region (0..4). MODE=spot (default) reads the pseudo-Visium spot protein matrix that the Xenium
benchmark uses (xenium_pseudovisium/data_protein_gt/h5ad_objects/Xenium_region_<r>_CITE.h5ad). MODE=cell reads the
single-cell-resolution protein matrix the April 2026 Module 1/2 runs used (load_xenium_singlecell), for comparison only.
Writes <out_dir>/region_<r>.json with the retained markers (same fields as the breast recount files, plus per-marker scores).

Usage (via run_module1_recount_xenium.sbatch):  run_module1_recount_xenium.py OUT_DIR   (env SLURM_ARRAY_TASK_ID = 0..4)
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy.sparse

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
from CITEgeist.model import identify_interesting_markers  # noqa: E402

region = int(os.environ["SLURM_ARRAY_TASK_ID"])
mode = os.environ.get("MODE", "spot")
if mode == "spot":
    import scanpy as sc

    ad = sc.read_h5ad(
        REPO / "Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects" / f"Xenium_region_{region}_CITE.h5ad"
    )
elif mode == "cell":
    from load_xenium_singlecell import load_xenium_singlecell

    _, ad = load_xenium_singlecell(region_id=region, seed=region)
else:
    raise ValueError(mode)
X = ad.X.toarray() if scipy.sparse.issparse(ad.X) else np.asarray(ad.X)
X = np.asarray(X, dtype=np.float64)
coords = np.asarray(ad.obsm["spatial"], dtype=np.float64)
markers = list(ad.var_names)
finite = np.all(np.isfinite(coords), axis=1)
X, coords = X[finite], coords[finite]
print(f"region {region} mode {mode}: {X.shape[0]} obs x {X.shape[1]} markers, dropped {int((~finite).sum())} non-finite coords", flush=True)

m = identify_interesting_markers(
    X=X, coords=coords, marker_names=markers, morans_k=8, morans_n_perm=99, seed=42, verbose=True
)
out = dict(
    sample=f"Xenium_region_{region}",
    mode=mode,
    n_spots=int(X.shape[0]),
    n_markers=len(markers),
    n_retained=len(m.interesting_markers),
    retained=list(m.interesting_markers),
    scores={
        k.name: dict(
            kurtosis=float(k.kurtosis),
            gmm_snr=float(k.gmm_snr),
            morans_i=float(k.morans_i),
            morans_pvalue=float(k.morans_i_pvalue),
            passed_kurtosis=bool(k.passed_kurtosis),
            passed_gmm=bool(k.passed_gmm),
            passed_morans=bool(k.passed_morans),
        )
        for k in m.markers
    },
    commit=subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
)
Path(sys.argv[1], f"region_{region}.json").write_text(json.dumps(out, indent=1))
print("DONE region", region, out["n_retained"], "of", out["n_markers"])
