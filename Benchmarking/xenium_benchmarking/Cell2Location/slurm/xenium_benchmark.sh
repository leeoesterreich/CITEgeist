#!/bin/bash
#SBATCH --job-name=c2l_xenium
#SBATCH --array=0-4
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100

# Run Cell2Location deconvolution on Xenium pseudo-Visium regions
# Array job: processes regions 0-4

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

# Load CUDA module for GPU support
module load cuda/11.8.0
echo "Loaded CUDA: $(nvcc --version 2>/dev/null | grep release || echo 'CUDA loaded')"

# Get region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

# Directories (use absolute paths - SCRIPT_DIR resolves incorrectly in SLURM)
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
C2L_DIR="${BASE_DIR}/Cell2Location"
SRC_DIR="${C2L_DIR}/src"
INPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
REF_DIR="${C2L_DIR}/reference_model_protein7"
OUTPUT_DIR="${C2L_DIR}/output_protein_gt"

echo "=============================================="
echo "Cell2Location Xenium Benchmark - Region ${REGION_ID}"
echo "=============================================="
echo "Start time: $(date)"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Reference directory: ${REF_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if reference model exists
if [ ! -d "${REF_DIR}/model" ]; then
    echo "ERROR: Reference model not found at ${REF_DIR}/model"
    echo "Please run train_reference.sh first"
    exit 1
fi

# Check if input data exists
if [ ! -f "${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad" ]; then
    echo "ERROR: Input data not found"
    echo "Expected: ${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad"
    exit 1
fi

# Activate Cell2Location environment (use export PATH + source activate pattern)
ENV_PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/cell2location_env"
MINICONDA_PATH="$(conda info --base)"
echo "Activating cell2location_env from ${ENV_PATH}..."
export PATH="${ENV_PATH}/bin:${MINICONDA_PATH}/bin:$PATH"
source activate "${ENV_PATH}"

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Check GPU
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}'); print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"
echo ""

# Run deconvolution
python "${SRC_DIR}/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --ref-dir "${REF_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-cells-per-location 5 \
    --max-epochs 30000 \
    --prefix "Xenium"

echo ""
echo "=============================================="
echo "Region ${REGION_ID} Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/Xenium_region_${REGION_ID}/"
