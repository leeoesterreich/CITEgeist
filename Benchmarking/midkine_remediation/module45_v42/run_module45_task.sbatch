#!/bin/bash
#SBATCH --job-name=m45_v2b251_task
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=3:00:00
# Array task: Module 4 + 4b for ONE specimen (SLURM_ARRAY_TASK_ID indexes run_module45.SPECIMENS). Max RSS of the biggest specimen in the serial run (job 12148263) was 8.1 GB.
# Stage B25 (stamped), submitted by .commands/submit-m45-stage.sh with env SHA, CODE_DIR (clean detached worktree), STAGE_STAMP=1.
# Data are read from MAIN (V2_NAMESPACE_ROOT); the out dir is results/<B25 run id>; the job refuses to start from a dirty worktree.
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_${SLURM_JOB_ID} NUMBA_CACHE_DIR=/tmp/numba_${SLURM_JOB_ID}
PY="$(cg_env_python CITEgeist_env)"
: "${CODE_DIR:?}" "${SHA:?}"
export V2_NAMESPACE_ROOT=$MAIN PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$CODE_DIR:$CODE_DIR/scripts
cd "$CODE_DIR"
[ "$(git rev-parse HEAD)" = "$SHA" ] && [ -z "$(git status --porcelain -uno)" ] || { echo "worktree is not clean at $SHA" >&2; exit 1; }
echo "stage=B25 task=${SLURM_ARRAY_TASK_ID} sha=$SHA host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/module45_v42/run_module45.py --stamp --specimen "${SLURM_ARRAY_TASK_ID}"
echo "done end=$(date -Is)"
