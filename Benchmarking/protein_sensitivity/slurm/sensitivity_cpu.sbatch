#!/bin/bash
#SBATCH --job-name=protein_sens_cpu
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --output=Benchmarking/protein_sensitivity/logs/%x_%A_%a.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
# CPU (htc) copy of sensitivity.sbatch for PHASE1=marker_share (no cuOPT, no GPU).
# Inputs via env (sbatch exports the caller's env): DATASET (xenium|sim), UNITS (space-separated), N_TYPES, CONDITIONS_JSON, OUT_ROOT, PHASE1, optional SIM_ROOT (sim replicate root; read by baseline.resolve_sim_root, unset = canonical replicates)
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
UNITS_ARR=($UNITS)
UNIT=${UNITS_ARR[$SLURM_ARRAY_TASK_ID]}
OUT="${OUT_ROOT}/${DATASET}_${N_TYPES}_${UNIT//:/_}"
echo "phase1=${PHASE1:-qp} out=$OUT sim_root=${SIM_ROOT:-default}"
mkdir -p "$OUT"
EXTRA=()
if [ "${SAVE_LAYERS:-0}" = 1 ]; then EXTRA=(--save-layers "$OUT/${UNIT%%:*}/Wu_rep_${UNIT##*:}"); fi
PYTHONPATH=Benchmarking $PY -u Benchmarking/protein_sensitivity/run_sensitivity.py \
  --dataset "$DATASET" --unit "$UNIT" --n-types "$N_TYPES" --conditions "$CONDITIONS_JSON" --out "$OUT" "${EXTRA[@]}"
