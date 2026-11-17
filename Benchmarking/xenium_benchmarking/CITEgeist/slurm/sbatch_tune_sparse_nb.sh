#!/bin/bash
#SBATCH --job-name=tune_snb
#SBATCH --output=tune_sparse_nb_%j.log
#SBATCH --error=tune_sparse_nb_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc

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
set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

export OMP_NUM_THREADS=1

cd ${CITEGEIST_REPO_ROOT}

python Benchmarking/xenium_benchmarking/CITEgeist/src/tune_sparse_nb_grid.py
