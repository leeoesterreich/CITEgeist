#!/bin/bash
#SBATCH --job-name=prep_erpos_ref
#SBATCH --output=Benchmarking/competitor_gex_validation/logs/01_prepare_reference_%j.out
#SBATCH --error=Benchmarking/competitor_gex_validation/logs/01_prepare_reference_%j.err
#SBATCH --partition=htc
#SBATCH --cluster=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1:00:00

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
set -eo pipefail

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
cd "${REPO_ROOT}"
export MPLCONFIGDIR=/tmp/mpl_${SLURM_JOB_ID}
export HDF5_USE_FILE_LOCKING=FALSE


echo "=== Preparing ER+ Wu 2021 reference ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "head=$(git rev-parse HEAD)"
echo "Start: $(date)"

"$(cg_env_python CITEgeist_env)" -u Benchmarking/competitor_gex_validation/prepare_erpos_reference.py

echo "=== Done: $(date) ==="
