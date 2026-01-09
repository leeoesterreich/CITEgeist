#!/bin/bash
#SBATCH --job-name=rctd_xenium
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/RCTD/slurm/slurm_log/xenium_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/RCTD/slurm/slurm_log/xenium_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Run RCTD deconvolution on Xenium pseudo-Visium regions
# Array job: processes regions 0-4
# Note: RCTD only outputs proportions, not GEX layers

set -e

# Get region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

# Directories - use absolute paths to avoid SLURM symlink issues
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
RCTD_DIR="${BASE_DIR}/RCTD"
SRC_DIR="${RCTD_DIR}/src"
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt"
REF_DIR="${BASE_DIR}/reference_data/GSE156632/processed_granular/rctd"
OUTPUT_DIR="${RCTD_DIR}/output_granular"
TEMP_DIR="${OUTPUT_DIR}/temp_csvs"

echo "=============================================="
echo "RCTD Xenium Benchmark - Region ${REGION_ID}"
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

# Environment paths (use source activate with full path)
CITEGEIST_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
R_ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/r_deconv_env"

# Step 1: Convert spatial h5ad to CSV (using Python)
echo "=============================================="
echo "Step 1: Converting spatial h5ad to CSV"
echo "=============================================="

# Activate CITEgeist_env for Python conversion (has scanpy)
echo "Activating CITEgeist_env..."
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

python "${SRC_DIR}/convert_h5ad.py" \
    --input "${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad" \
    --output-dir "${TEMP_DIR}" \
    --mode spatial \
    --prefix "Xenium_region_${REGION_ID}"

# Step 2: Run RCTD (using R)
echo ""
echo "=============================================="
echo "Step 2: Running RCTD"
echo "=============================================="

# Activate R environment (using full path as specified)
echo "Activating r_deconv_env from ${R_ENV_PATH}..."
export PATH="${R_ENV_PATH}/bin:$PATH"
source activate "${R_ENV_PATH}"

echo "R: $(which R)"
echo "R version: $(R --version | head -1)"
echo ""

Rscript "${SRC_DIR}/run_benchmark.R" \
    --region-id ${REGION_ID} \
    --ref-counts "${REF_DIR}/counts.csv" \
    --ref-cell-types "${REF_DIR}/cell_types.csv" \
    --spatial-counts "${TEMP_DIR}/Xenium_region_${REGION_ID}_counts.csv" \
    --spatial-coords "${TEMP_DIR}/Xenium_region_${REGION_ID}_coords.csv" \
    --output-dir "${OUTPUT_DIR}" \
    --prefix "Xenium" \
    --max-cores ${SLURM_CPUS_PER_TASK}

echo ""
echo "=============================================="
echo "Region ${REGION_ID} Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/Xenium_region_${REGION_ID}/"
