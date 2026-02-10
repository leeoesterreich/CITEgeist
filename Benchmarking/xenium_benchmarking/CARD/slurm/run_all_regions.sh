#!/bin/bash
#SBATCH --job-name=card_xenium
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm/slurm_log/xenium_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm/slurm_log/xenium_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Run CARD deconvolution on Xenium pseudo-Visium regions (REFERENCE mode)
# Array job: processes regions 0-4
#
# Prerequisites:
#   1. CARD_env conda environment installed (see Benchmarking/environments/setup_card_env.sh)
#   2. Reference data processed (uses same as RCTD)

set -e

# Get region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

# Directories
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
CARD_DIR="${BASE_DIR}/CARD"
SRC_DIR="${CARD_DIR}/src"
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt"
REF_DIR="${BASE_DIR}/reference_data/GSE156632/processed_protein7/rctd"
OUTPUT_DIR="${CARD_DIR}/output_protein_gt"
TEMP_DIR="${OUTPUT_DIR}/temp_csvs"

echo "=============================================="
echo "CARD Xenium Benchmark - Region ${REGION_ID} (Reference mode)"
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
    echo "Please ensure RCTD reference data is processed first"
    exit 1
fi

# Check if input data exists
if [ ! -f "${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad" ]; then
    echo "ERROR: Input data not found"
    echo "Expected: ${INPUT_DIR}/h5ad_objects/Xenium_region_${REGION_ID}_GEX.h5ad"
    exit 1
fi

# Environment paths
CITEGEIST_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"
CARD_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CARD_env"

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

# Step 2: Run CARD (using R)
echo ""
echo "=============================================="
echo "Step 2: Running CARD (Reference mode)"
echo "=============================================="

# Activate CARD environment
echo "Activating CARD_env from ${CARD_ENV}..."
export PATH="${CARD_ENV}/bin:$PATH"
source activate "${CARD_ENV}"

echo "R: $(which R)"
echo "R version: $(R --version | head -1)"
echo ""

Rscript "${SRC_DIR}/run_benchmark.R" \
    --mode reference \
    --region-id ${REGION_ID} \
    --ref-counts "${REF_DIR}/counts.csv" \
    --ref-cell-types "${REF_DIR}/cell_types.csv" \
    --spatial-counts "${TEMP_DIR}/Xenium_region_${REGION_ID}_counts.csv" \
    --spatial-coords "${TEMP_DIR}/Xenium_region_${REGION_ID}_coords.csv" \
    --output-dir "${OUTPUT_DIR}" \
    --prefix "Xenium"

echo ""
echo "=============================================="
echo "Region ${REGION_ID} Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/Xenium_region_${REGION_ID}/"
