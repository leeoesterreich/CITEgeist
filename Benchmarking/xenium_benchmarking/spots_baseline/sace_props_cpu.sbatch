#!/bin/bash
#SBATCH --job-name=sace_props_cpu
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
# SACE with SPOTS (a) and ground-truth proportions (CPU arms), one region per array task. Env: CODE_DIR (clean detached worktree), OUT_DIR; optional PROPS_DIR + PROFILE (marker_share arm: props tree and achievable|ms_v41).
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
cd "$CODE_DIR"
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
echo "code=$CODE_DIR head=$(git rev-parse HEAD) out=$OUT_DIR region=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is)"
for ARM in ${ARMS:-qp_saved}; do "$PY" -u Benchmarking/xenium_benchmarking/spots_baseline/sace_props_ablation.py "$ARM" "$SLURM_ARRAY_TASK_ID" "$OUT_DIR" ${PROPS_DIR:+"$PROPS_DIR" "${PROFILE:-achievable}"}; done
echo "done end=$(date -Is)"
