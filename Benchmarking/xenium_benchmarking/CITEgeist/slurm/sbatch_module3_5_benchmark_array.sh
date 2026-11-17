#!/bin/bash
#SBATCH --job-name=module3_5_bench
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --gres=gpu:1

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
SRC_SCRIPT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/run_module3_5_xenium_benchmark.py"
OUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_module3_5_benchmark_array"
LOG_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log"

mkdir -p "${OUT_DIR}"
mkdir -p "${LOG_DIR}"

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

REGION_ID="${SLURM_ARRAY_TASK_ID}"

cd "${REPO_ROOT}"

NUMBA_CACHE_DIR=/tmp/numba_cache \
python "${SRC_SCRIPT}" \
  --region-id "${REGION_ID}" \
  --output-dir "${OUT_DIR}"

echo "Module 3.5 benchmark complete for region ${REGION_ID}"
