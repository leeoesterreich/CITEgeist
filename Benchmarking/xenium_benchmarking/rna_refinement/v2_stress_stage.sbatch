#!/bin/bash
#SBATCH --job-name=v2stress
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --exclude=htc-1024-n2,htc-n37
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# v42 plan Task 9 stages. Env: STAGE (x7 x7score x8 wh whscore), SHA (the commit that holds this code), WH_OUT (wh, whscore: the NEW run dir).
# Submitted by .commands/submit-v2-stress.sh. x7 and wh are arrays (SLURM_ARRAY_TASK_ID = unit: x7 0-9 = 6-type regions 0-4 then 7-type
# regions 0-4; wh 0-4 = region). Same code proof as v2_official_stage.sbatch (x7 stages also prove Benchmarking/protein_sensitivity).
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
: "${STAGE:?set STAGE}"; : "${SHA:?set SHA to the commit that holds this code}"
PY="$(cg_env_python CITEgeist_env)"
REPO=${CITEGEIST_REPO_ROOT}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID PYTHONPATH=$REPO:$REPO/Benchmarking
cd "$REPO"
RF=Benchmarking/xenium_benchmarking/rna_refinement
STRESS=
case "$STAGE" in x7|x7score) STRESS=--stress ;; esac
"$PY" -u $RF/v2_official.py check-code "$SHA" $STRESS
echo "stage=$STAGE sha=$SHA host=$(hostname) task=${SLURM_ARRAY_TASK_ID:-none} start=$(date -Is)"
case "$STAGE" in
  x7)       "$PY" -u $RF/v2_stress.py sweep "${SLURM_ARRAY_TASK_ID:?x7 is an array job}" ;;
  x7score)  "$PY" -u $RF/v2_stress.py score ;;
  x8)       "$PY" -u $RF/v2_heldout.py x8
            "$PY" -u $RF/v2_stress.py complete x8 ;;
  wh)       "$PY" -u $RF/v2_heldout.py wh-region "${WH_OUT:?set WH_OUT}" "${SLURM_ARRAY_TASK_ID:?wh is an array job}" ;;
  whscore)  "$PY" -u $RF/v2_heldout.py wh-score "${WH_OUT:?set WH_OUT}" ;;
  *)        echo "unknown STAGE $STAGE" >&2; exit 2 ;;
esac
echo "stage=$STAGE done $(date -Is)"
