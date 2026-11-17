#!/bin/bash
#SBATCH --job-name=m2_candidates
#SBATCH --array=3-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
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

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
SRC_SCRIPT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/generate_module2_candidate_reviews.py"
OUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_module2_region_candidates_fullres"

mkdir -p "${OUT_DIR}"
mkdir -p "${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log"

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

export CITEGEIST_N_JOBS="${SLURM_CPUS_PER_TASK}"
export NUMBA_CACHE_DIR=/tmp/numba_cache

cd "${REPO_ROOT}"

python "${SRC_SCRIPT}" \
  --output-dir "${OUT_DIR}" \
  --region-id "${SLURM_ARRAY_TASK_ID}"

echo "Module 2 candidates complete for region ${SLURM_ARRAY_TASK_ID}"
