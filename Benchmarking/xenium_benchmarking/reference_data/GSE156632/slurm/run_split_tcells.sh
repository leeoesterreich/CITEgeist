#!/bin/bash
#SBATCH --job-name=split_tcells
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../../_paths.sh"
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
set -e

SRC_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/reference_data/GSE156632/src"
ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"

echo "=============================================="
echo "Splitting T cells to derive 7-type reference"
echo "=============================================="
echo "Start time: $(date)"

mkdir -p "${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/slurm_log"

export PATH="${ENV}/bin:$PATH"
source activate "${ENV}"

python "${SRC_DIR}/split_tcells_to_protein7.py"

echo "=============================================="
echo "Done. End time: $(date)"
