#!/bin/bash
#SBATCH --job-name=card_reffree_xenium
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm/slurm_log/xenium_reffree_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm/slurm_log/xenium_reffree_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Run CARD deconvolution on Xenium pseudo-Visium regions (REFERENCE-FREE mode)
# Array job: processes regions 0-4
#
# Prerequisites:
#   1. CARD_env conda environment installed
#   2. Marker genes generated (run setup_markers.sh first)

set -e

# Get region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

# Directories
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
CARD_DIR="${BASE_DIR}/CARD"
SRC_DIR="${CARD_DIR}/src"
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt"
MARKER_FILE="${CARD_DIR}/marker_genes_protein7.csv"
OUTPUT_DIR="${CARD_DIR}/output_protein_gt_reffree"
TEMP_DIR="${OUTPUT_DIR}/temp_csvs"

echo "=============================================="
echo "CARD Xenium Benchmark - Region ${REGION_ID} (Reference-FREE mode)"
echo "=============================================="
echo "Start time: $(date)"
echo "Array Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo ""
echo "Input directory: ${INPUT_DIR}"
echo "Marker genes: ${MARKER_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Create directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${TEMP_DIR}"

# Check if marker genes exist
if [ ! -f "${MARKER_FILE}" ]; then
    echo "ERROR: Marker genes not found at ${MARKER_FILE}"
    echo "Please run setup_markers.sh first"
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

# Step 2: Run CARD-Free (using R)
echo ""
echo "=============================================="
echo "Step 2: Running CARD (Reference-FREE mode)"
echo "=============================================="

# Activate CARD environment
echo "Activating CARD_env from ${CARD_ENV}..."
export PATH="${CARD_ENV}/bin:$PATH"
source activate "${CARD_ENV}"

echo "R: $(which R)"
echo "R version: $(R --version | head -1)"
echo ""

Rscript "${SRC_DIR}/run_benchmark.R" \
    --mode reference_free \
    --region-id ${REGION_ID} \
    --marker-genes "${MARKER_FILE}" \
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
