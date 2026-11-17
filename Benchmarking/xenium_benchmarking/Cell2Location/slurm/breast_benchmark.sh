#!/bin/bash
#SBATCH --job-name=c2l_breast
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --array=0-4
#SBATCH --output=slurm/logs/c2l_breast_%A_%a.log

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
module purge
module load gcc/12.2.0

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
BASE_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
REGION_ID=${SLURM_ARRAY_TASK_ID}

source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/cell2location_env

python "${BASE_DIR}/Cell2Location/src/run_benchmark.py" \
    --region-id "${REGION_ID}" \
    --input-dir "${BASE_DIR}/data_breast" \
    --ref-dir "${BASE_DIR}/reference_data/Wu_2021/cell2location" \
    --output-dir "${BASE_DIR}/Cell2Location/output_breast" \
    --prefix "Xenium"

echo "C2L breast region ${REGION_ID} complete"
