#!/bin/bash
#SBATCH --job-name=psv41p
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# v41p protein-sensitivity run (marker share + GMM brightness refit + ms_v41 profile) on htc. Copy of sensitivity_cpu.sbatch.
# Inputs via env: UNITS (space-separated regions), N_TYPES, CONDITIONS_JSON, OUT_ROOT, CODE_DIR (clean worktree).
# Unit dir pattern (must match strict_score.unit_dir): $OUT_ROOT/citegeist/xenium_<n>_<region>
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
# node-local temp: baseline.run_xenium_qp uses tempfile.TemporaryDirectory(); on NFS its cleanup raised OSError(39)
export TMPDIR=${SLURM_TMPDIR:-/tmp}/psv41p_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}; mkdir -p "$TMPDIR"
export PHASE1=marker_share
PY="$(cg_env_python CITEgeist_env)"
UNITS_ARR=($UNITS)
UNIT=${UNITS_ARR[$SLURM_ARRAY_TASK_ID]}
OUT="${OUT_ROOT:?}/citegeist/xenium_${N_TYPES}_${UNIT}"
echo "profile=ms_v41 brightness=1 phase1=$PHASE1 head=$(git rev-parse HEAD) tmpdir=$TMPDIR out=$OUT"
mkdir -p "$OUT"
PYTHONPATH=Benchmarking $PY -u Benchmarking/protein_sensitivity/run_sensitivity.py \
  --dataset xenium --unit "$UNIT" --n-types "$N_TYPES" --conditions "$CONDITIONS_JSON" --out "$OUT" \
  --profile ms_v41 --brightness 1
