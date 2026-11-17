#!/bin/bash
#SBATCH --job-name=citegeist_benchmark
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --array=0-4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=slurm/logs/benchmark_%A_%a.out
#SBATCH --error=slurm/logs/benchmark_%A_%a.err

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
set -u  # Re-enable nounset after conda activation

module load gurobi/12.0.3

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
INPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
PREPROCESS_DIR="${SCRIPT_DIR}/output/preprocessing"
REGION_ID=${SLURM_ARRAY_TASK_ID}

echo "Running region ${REGION_ID} at $(date)"

# Run 6-type
echo "=== 6-type ==="
python "${SCRIPT_DIR}/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --n-types 6 \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${SCRIPT_DIR}/output_singler_6type" \
    --nuclei-counts-dir "${PREPROCESS_DIR}/nuclei_counts" \
    --features-dir "${PREPROCESS_DIR}/features" \
    --run-gex

# Run 7-type
echo "=== 7-type ==="
python "${SCRIPT_DIR}/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --n-types 7 \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${SCRIPT_DIR}/output_singler_7type" \
    --nuclei-counts-dir "${PREPROCESS_DIR}/nuclei_counts" \
    --features-dir "${PREPROCESS_DIR}/features" \
    --run-gex

echo "Region ${REGION_ID} complete at $(date)"
