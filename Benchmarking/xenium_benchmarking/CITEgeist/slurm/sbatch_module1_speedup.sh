#!/bin/bash
#SBATCH --job-name=m1_speedup
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Speedup benchmark for Module 1 (marker_interest.py) parallelization.
# Runs identify_interesting_markers twice on the same data:
#   pass 1: CITEGEIST_N_JOBS=1  (serial baseline)
#   pass 2: CITEGEIST_N_JOBS=16 (parallel, uses joblib)
# Reports wall-clock seconds and speedup factor.

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
set -eo pipefail

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
OUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_module1_speedup"
mkdir -p "${OUT_DIR}"
mkdir -p "${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log"

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

export NUMBA_CACHE_DIR=/tmp/numba_cache

cd "${REPO_ROOT}"

REGION_ID="${SLURM_ARRAY_TASK_ID}"
export OUT_JSON="${OUT_DIR}/region_${REGION_ID}_speedup.json"

python - <<'PYEOF'
import json, os, sys, time
from pathlib import Path

REPO_ROOT = Path(__import__("os").environ["CITEGEIST_REPO_ROOT"])
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

import numpy as np
from load_xenium_singlecell import load_xenium_singlecell
from CITEgeist.model.marker_interest import identify_interesting_markers

region_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
out_json = Path(os.environ["OUT_JSON"])
n_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "16"))

print(f"Loading region {region_id}...", flush=True)
_, adata_protein = load_xenium_singlecell(region_id=region_id)
X = np.asarray(adata_protein.X.todense() if hasattr(adata_protein.X, "todense") else adata_protein.X, dtype=np.float64)
coords = np.asarray(adata_protein.obsm["spatial"], dtype=np.float64)
marker_names = list(adata_protein.var_names)
print(f"Data: {X.shape[0]} cells × {X.shape[1]} markers", flush=True)

# --- serial pass ---
os.environ["CITEGEIST_N_JOBS"] = "1"
print("Running M1 serial (n_jobs=1)...", flush=True)
t0 = time.perf_counter()
identify_interesting_markers(X=X, coords=coords, marker_names=marker_names,
                              morans_k=15, morans_n_perm=199, verbose=False)
serial_seconds = time.perf_counter() - t0
print(f"  serial: {serial_seconds:.1f}s", flush=True)

# --- parallel pass ---
os.environ["CITEGEIST_N_JOBS"] = str(n_cpus)
print(f"Running M1 parallel (n_jobs={n_cpus})...", flush=True)
t0 = time.perf_counter()
identify_interesting_markers(X=X, coords=coords, marker_names=marker_names,
                              morans_k=15, morans_n_perm=199, verbose=False)
parallel_seconds = time.perf_counter() - t0
print(f"  parallel: {parallel_seconds:.1f}s", flush=True)

payload = {
    "region_id": region_id,
    "n_cells": int(X.shape[0]),
    "n_markers": int(X.shape[1]),
    "n_cpus": n_cpus,
    "serial_seconds": serial_seconds,
    "parallel_seconds": parallel_seconds,
    "speedup_x": serial_seconds / parallel_seconds if parallel_seconds > 0 else None,
}
out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps(payload, indent=2, sort_keys=True), flush=True)
PYEOF

echo "M1 speedup result for region ${REGION_ID} written to ${OUT_JSON}"
