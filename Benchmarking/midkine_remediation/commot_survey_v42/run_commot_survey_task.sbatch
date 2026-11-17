#!/bin/bash
#SBATCH --job-name=cmsurvey_b26_task
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
# Array task of stage B26: task 0-11 = level 1 (specimens 0-11), task 12-23 = level 2 (specimens 0-11). P3-S2 level 1 took 4h47m in the survey, hence 8 h.
# Submitted by .commands/submit-commot-survey-stage.sh with env SHA, CODE_DIR (clean detached worktree), MAIN. COMMOT env python (commot 0.0.3 needs numpy 1.x).
# Data are read from MAIN (V2_NAMESPACE_ROOT); the out dir is results/<B26 run id>; the job refuses to start from a dirty worktree.
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
MAIN=${CITEGEIST_REPO_ROOT}
cd "$MAIN"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_${SLURM_JOB_ID} NUMBA_CACHE_DIR=/tmp/numba_${SLURM_JOB_ID} OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
PY="$(cg_env_python COMMOT)"
: "${CODE_DIR:?}" "${SHA:?}"
export V2_NAMESPACE_ROOT=$MAIN PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$CODE_DIR:$CODE_DIR/scripts
cd "$CODE_DIR"
[ "$(git rev-parse HEAD)" = "$SHA" ] && [ -z "$(git status --porcelain -uno)" ] || { echo "worktree is not clean at $SHA" >&2; exit 1; }
T=${SLURM_ARRAY_TASK_ID}
if [ "$T" -lt 12 ]; then LEVEL=l1; else LEVEL=l2; fi
echo "stage=B26 task=$T level=$LEVEL specimen=$((T % 12)) sha=$SHA host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/commot_survey_v42/commot_survey.py --level "$LEVEL" --specimen "$((T % 12))"
echo "done end=$(date -Is)"
