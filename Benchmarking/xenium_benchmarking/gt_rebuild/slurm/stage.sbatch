#!/usr/bin/env bash
#SBATCH --job-name=gt_rebuild
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=06:00:00
# One GT-rebuild stage per job (STAGE=stability_variant: a 3-task array, one variant per task). Submit via .commands/submit-gt-rebuild.sh (it sets --export=ALL,STAGE=<stage>).
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
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
: "${STAGE:?STAGE not set; submit with .commands/submit-gt-rebuild.sh}"
export MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID
export HDF5_USE_FILE_LOCKING=FALSE
export PY="$(cg_env_python CITEgeist_env)"
export RSCRIPT=${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env/bin/Rscript
mkdir -p "$MPLCONFIGDIR"
echo "gt_rebuild stage=$STAGE job=$SLURM_JOB_ID task=${SLURM_ARRAY_TASK_ID:-none} root=${GT_REBUILD_OUT_ROOT:-default} $(date -Is)"
bash ${CITEGEIST_REPO_ROOT}/.commands/gt-rebuild-run.sh "$STAGE"
echo "gt_rebuild stage=$STAGE done $(date -Is)"
