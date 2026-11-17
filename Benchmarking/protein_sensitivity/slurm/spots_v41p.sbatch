#!/bin/bash
#SBATCH --job-name=psv41p_spots
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# v41p SPOTS (reading a) through the SAME perturbed matrices as the CITEgeist sweep. Array 0-4 = regions; loops N in 6 7.
# Env: OUT_ROOT, CODE_DIR (clean worktree). Reuses $OUT_ROOT/conditions_xenium_${N}.json from the sweep.
# Unit dir pattern (must match strict_score.unit_dir): $OUT_ROOT/spots/spots_<n>_<region>
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
export TMPDIR=${SLURM_TMPDIR:-/tmp}/psv41p_spots_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}; mkdir -p "$TMPDIR"
PY="$(cg_env_python CITEgeist_env)"
UNIT=${SLURM_ARRAY_TASK_ID:?array task id required}
echo "profile=ms_v41 zero_markers=fail head=$(git rev-parse HEAD) tmpdir=$TMPDIR"
for N in 6 7; do
  OUT="${OUT_ROOT:?}/spots/spots_${N}_${UNIT}"
  echo "n_types=$N unit=$UNIT out=$OUT"
  mkdir -p "$OUT"
  PYTHONPATH=Benchmarking $PY -u Benchmarking/protein_sensitivity/spots_sensitivity.py run-unit "$OUT" "$N" "$UNIT" \
    --conditions "$OUT_ROOT/conditions_xenium_${N}.json" --profile ms_v41
done
