#!/bin/bash
#SBATCH --job-name=a12_11_nmf
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=112G
#SBATCH --time=08:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_11_%j.log
# A12.11 Module 4 (NMF) program counts from the frozen path, over the 12 module3_7type specimens (htc CPU:
# the NMF / Moran path is numpy/sklearn-only; no GPU code path). Single job -> one merged nmf_counts.json.
# Env: RUN_ID, OUT_ROOT (= output/nmf_counts_module3_7type/<protocol_sha12>/<RUN_ID>).
# Required env MODULE3_DIR (stamped B1 dir) and LAYERS_RUN_ID (the B4 id): neither has a default (the old defaults output/module3_7type and r1 are v40/v41).
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../_paths.sh"
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
set -euo pipefail
cd "${CODE_DIR:-${CITEGEIST_REPO_ROOT}}"
MAIN=${CITEGEIST_REPO_ROOT}   # data/outputs stay in the main tree; CODE_DIR (if set) is a clean detached worktree
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
M3=(--module3-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}" --layers-root "$MAIN/output/sace_layers_module3_7type" --layers-run-id "${LAYERS_RUN_ID:?LAYERS_RUN_ID (the B4 id) is required}")
echo "run_id=$RUN_ID out=$OUT_ROOT host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/nmf_module4.py \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" --protocol Benchmarking/midkine_remediation/protocol.json "${M3[@]}"
echo "done end=$(date -Is)"
