#!/bin/bash
#SBATCH --job-name=a12_6_sace
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-11
#SBATCH --cpus-per-task=8
#SBATCH --mem=112G
#SBATCH --time=06:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_6_%A_%a.log
# v2 B4 SACE layers for the stamped B1 proportions, one specimen per array task (htc CPU: run_sace_layers is numpy-only, no GPU code path).
# Env: CODE_DIR (worktree to run from; default main tree), RUN_ID, OUT_ROOT (= output/sace_layers_module3_7type/<protocol_sha12>/<RUN_ID>), DATA_ROOT.
# Required env MODULE3_DIR (the B1 output dir) and NUCLEI_DIR (the stamped B2 artifact dir, <root>/v2b2); neither has a default.
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
SPECS=(HCC22-088-P1-S1 HCC22-088-P1-S2 HCC22-088-P2-S1 HCC22-088-P2-S2 HCC22-088-P3-S1_A HCC22-088-P3-S2 HCC22-088-P4-S1 HCC22-088-P4-S2_1i_rep HCC22-088-P5-S1 HCC22-088-P5-S2_F_rep HCC22-088-P6-S1 HCC22-088-P6-S2_D)
SPECIMEN=${SPECS[$SLURM_ARRAY_TASK_ID]}
M3=(--module3-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}" --nuclei-dir "${NUCLEI_DIR:?NUCLEI_DIR (the B2 artifact dir) is required}")
echo "specimen=$SPECIMEN out=$OUT_ROOT/$SPECIMEN host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/sace_layers.py --specimen "$SPECIMEN" --data-root "$DATA_ROOT" \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" --protocol Benchmarking/midkine_remediation/protocol.json "${M3[@]}"
echo "done $SPECIMEN end=$(date -Is)"
