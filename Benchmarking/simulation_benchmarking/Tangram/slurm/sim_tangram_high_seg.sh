#!/bin/bash
#SBATCH --job-name=tg_sim_high_seg
#SBATCH --array=0-4
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
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
set -e

REPLICATE_ID=${SLURM_ARRAY_TASK_ID}
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/Tangram"
SRC_DIR="${BASE_DIR}/src"
INPUT_PATH="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/replicates/high_seg/h5ad_objects/Wu_rep_${REPLICATE_ID}_GEX.h5ad"
WU_REF="${CITEGEIST_DATA_ROOT}/BC-Datasets/Wu_2021/Wu_etal_2021_BRCA_scRNASeq"
# Use GPU_RACE_NAMESPACE for safe writes; rename to canonical path at end
RACE_NS="${GPU_RACE_NAMESPACE:-default}"
TEMP_OUTPUT_DIR="${BASE_DIR}/high_seg/.race_tmp/${RACE_NS}/tangram_map_${REPLICATE_ID}"
FINAL_OUTPUT_DIR="${BASE_DIR}/high_seg/tangram_map_${REPLICATE_ID}"

echo "=============================================="
echo "Tangram Simulation Benchmark - high_seg rep ${REPLICATE_ID}"
echo "=============================================="
echo "Start time: $(date)"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID} | Task: ${SLURM_ARRAY_TASK_ID}"
echo "Race namespace: ${RACE_NS}"
echo ""
echo "Input:     ${INPUT_PATH}"
echo "Reference: ${WU_REF}"
echo "Output:    ${FINAL_OUTPUT_DIR}"
echo ""

if [ ! -f "${INPUT_PATH}" ]; then
    echo "ERROR: Input data not found at ${INPUT_PATH}"
    exit 1
fi

if [ ! -d "${WU_REF}" ]; then
    echo "ERROR: Wu reference not found at ${WU_REF}"
    exit 1
fi

echo "Activating tangram environment..."
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/tangram_env
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate tangram-env" >&2
    exit 1
fi

echo "Python: $(which python)"
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}')"
echo ""

python "${SRC_DIR}/run_benchmark.py" \
    --replicate-id ${REPLICATE_ID} \
    --condition high_seg \
    --input-path "${INPUT_PATH}" \
    --wu-reference-path "${WU_REF}" \
    --output-dir "${TEMP_OUTPUT_DIR}" \
    --n-markers 100 \
    --num-epochs 500

# Move to canonical output path (safe: only winner survives to this point)
mv "${TEMP_OUTPUT_DIR}" "${FINAL_OUTPUT_DIR}"

echo ""
echo "=============================================="
echo "Rep ${REPLICATE_ID} Complete! End time: $(date)"
echo "=============================================="
ls -la "${FINAL_OUTPUT_DIR}/"
