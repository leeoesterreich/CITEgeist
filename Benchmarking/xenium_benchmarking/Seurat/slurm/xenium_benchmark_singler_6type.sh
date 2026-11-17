#!/bin/bash
#SBATCH --job-name=seurat_xenium_6t
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Run Seurat TransferData deconvolution on Xenium pseudo-Visium regions
# Array job: processes regions 0-4
# SingleR 6-type reference (merged T cells)
# Note: Seurat only outputs proportions, not GEX layers

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

# Directories - use absolute paths
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking"
SEURAT_DIR="${BASE_DIR}/Seurat"
SRC_DIR="${SEURAT_DIR}/src"
INPUT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
REF_DIR="${REF_DIR:-${BASE_DIR}/reference_data/GSE156632/processed/seurat}"
OUTPUT_DIR="${OUTPUT_DIR:-${SEURAT_DIR}/output_singler_6type}"
TEMP_DIR="${OUTPUT_DIR}/temp_csvs"

# Use environment binaries directly (no conda activation needed)
PYTHON="$(cg_env_python CITEgeist_env)"
RSCRIPT="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/r_deconv_env/bin/Rscript"

echo "=============================================="
echo "Seurat Xenium Benchmark - Region ${REGION_ID}"
echo "SingleR 6-type reference (merged T cells)"
echo "=============================================="
echo "Start time: $(date)"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Reference directory: ${REF_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"

# Check if reference exists
if [ ! -f "${REF_DIR}/counts.csv" ]; then
    echo "ERROR: Reference counts not found at ${REF_DIR}/counts.csv"
    echo "Please run reference_data/GSE156632/slurm/process_reference.sh first"
    exit 1
fi

# Check if input data exists
if [ ! -f "${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad" ]; then
    echo "ERROR: Input data not found"
    echo "Expected: ${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad"
    exit 1
fi

# Step 1: Convert spatial h5ad to CSV (using Python)
echo "=============================================="
echo "Step 1: Converting spatial h5ad to CSV"
echo "=============================================="

echo "Python: ${PYTHON}"

${PYTHON} "${SRC_DIR}/convert_h5ad.py" \
    --input "${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad" \
    --output-dir "${TEMP_DIR}" \
    --mode spatial \
    --prefix "Xenium_region_${REGION_ID}"

# Step 2: Run Seurat (using R)
echo ""
echo "=============================================="
echo "Step 2: Running Seurat TransferData"
echo "=============================================="

echo "Rscript: ${RSCRIPT}"
echo ""

${RSCRIPT} "${SRC_DIR}/run_benchmark.R" \
    --region-id ${REGION_ID} \
    --ref-counts "${REF_DIR}/counts.csv" \
    --ref-cell-types "${REF_DIR}/cell_types.csv" \
    --spatial-counts "${TEMP_DIR}/Xenium_region_${REGION_ID}_counts.csv" \
    --output-dir "${OUTPUT_DIR}" \
    --prefix "Xenium" \
    --n-features 2000

echo ""
echo "=============================================="
echo "Region ${REGION_ID} Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/Xenium_region_${REGION_ID}/"
