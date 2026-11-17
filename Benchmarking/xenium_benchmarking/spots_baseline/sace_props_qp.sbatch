#!/bin/bash
#SBATCH --job-name=sace_props_qp
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:1
#SBATCH --array=0-4
# SACE with canonical cuOPT QP proportions (GPU arm; cuOPT has no CPU fallback). Submit via gpu-race-submit.sh. Env: CODE_DIR, OUT_DIR.
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
echo "race namespace=${GPU_RACE_NAMESPACE:-none} partition=${GPU_RACE_PARTITION:-none} gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo none)"
echo "code=$CODE_DIR head=$(git rev-parse HEAD) out=$OUT_DIR region=$SLURM_ARRAY_TASK_ID host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/xenium_benchmarking/spots_baseline/sace_props_ablation.py qp "$SLURM_ARRAY_TASK_ID" "$OUT_DIR"
echo "done end=$(date -Is)"
