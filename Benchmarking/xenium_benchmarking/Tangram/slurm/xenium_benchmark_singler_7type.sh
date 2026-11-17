#!/bin/bash
#SBATCH --job-name=tg_xenium_7t
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s

# Run Tangram deconvolution on Xenium pseudo-Visium regions
# Array job: processes regions 0-4
# SingleR 7-type reference (separate CD4/CD8 T cells)

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

# Get region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

# Directories (use absolute paths - SCRIPT_DIR resolves incorrectly in SLURM)
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
TG_DIR="${BASE_DIR}/Tangram"
SRC_DIR="${TG_DIR}/src"
INPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
REF_PATH="${REF_PATH:-${BASE_DIR}/reference_data/GSE156632/processed_protein7/tangram/reference.h5ad}"
OUTPUT_DIR="${OUTPUT_DIR:-${TG_DIR}/output_singler_7type}"

echo "=============================================="
echo "Tangram Xenium Benchmark - Region ${REGION_ID}"
echo "SingleR 7-type reference (separate CD4/CD8 T cells)"
echo "=============================================="
echo "Start time: $(date)"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Reference: ${REF_PATH}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if reference exists
if [ ! -f "${REF_PATH}" ]; then
    echo "ERROR: Reference data not found at ${REF_PATH}"
    echo "Please run reference_data/GSE156632/slurm/process_reference.sh first"
    exit 1
fi

# Check if input data exists
if [ ! -f "${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad" ]; then
    echo "ERROR: Input data not found"
    echo "Expected: ${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad"
    exit 1
fi

# Activate Tangram environment
# Using the Tangram environment (verified to work with GPU)
echo "Activating tangram environment..."
export PATH=$(conda info --base)/bin:$PATH
source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/tangram-env
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to activate tangram environment"
    exit 1
fi

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
    --ref-path "${REF_PATH}" \
    --output-dir "${OUTPUT_DIR}" \
    --n-markers 100 \
    --num-epochs 500 \
    --prefix "Xenium"

echo ""
echo "=============================================="
echo "Region ${REGION_ID} Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/Xenium_region_${REGION_ID}/"
