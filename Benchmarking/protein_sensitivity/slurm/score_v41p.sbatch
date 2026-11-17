#!/bin/bash
#SBATCH --job-name=psv41p_score
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# v41p scoring on htc. SCORE_MODE selects the module: gate (gate-score stage) or score (re-gate the final root, then strict scoring). Env: SCORE_MODE, OUT_ROOT, CODE_DIR.
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
cd "${CODE_DIR:?CODE_DIR is required (clean worktree)}"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID
export TMPDIR=${SLURM_TMPDIR:-/tmp}/psv41p_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}; mkdir -p "$TMPDIR"
PY="$(cg_env_python CITEgeist_env)"
echo "score_mode=${SCORE_MODE:?} head=$(git rev-parse HEAD) tmpdir=$TMPDIR"
case "$SCORE_MODE" in
  gate) PYTHONPATH=Benchmarking $PY -u -m protein_sensitivity.gate --results "${OUT_ROOT:?}" --citegeist-subdir citegeist --n-types 6 7 ;;
  score) PYTHONPATH=Benchmarking $PY -u -m protein_sensitivity.aggregate_v41p --root "${OUT_ROOT:?}" --n-types 6 7 ;;
  *) echo "SCORE_MODE=$SCORE_MODE not implemented yet" >&2; exit 2 ;;
esac
