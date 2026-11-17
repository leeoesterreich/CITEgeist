#!/usr/bin/env bash
#SBATCH --job-name=rnaprot_gt
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=03:00:00
# rnaprot GT build. Submit via .commands/submit-rnaprot-gt.sh (sets OUT).
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
: "${OUT:?OUT not set; submit with .commands/submit-rnaprot-gt.sh}"
PY=${PY:-"$(cg_env_python CITEgeist_env)"}
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONPATH=${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking
cd "$PYTHONPATH"
echo "rnaprot build job=$SLURM_JOB_ID out=$OUT $(date -Is)"
$PY -u -m gt_rebuild.rnaprot build --out-root "$OUT"
echo "rnaprot build done $(date -Is)"
