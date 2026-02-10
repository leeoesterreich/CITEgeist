#!/bin/bash
#SBATCH --job-name=card_setup_markers
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm/slurm_log/setup_markers_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm/slurm_log/setup_markers_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Generate marker genes for CARD reference-free mode
# Run this once before running run_all_regions_reffree.sh
#
# Uses the reference scRNA-seq data to identify top marker genes per cell type.

set -e

# Directories
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
CARD_DIR="${BASE_DIR}/CARD"
SRC_DIR="${CARD_DIR}/src"
REF_FILE="${BASE_DIR}/reference_data/GSE156632/processed_protein7/rctd/reference.h5ad"
OUTPUT_FILE="${CARD_DIR}/marker_genes_protein7.csv"

echo "=============================================="
echo "CARD Marker Gene Generation"
echo "=============================================="
echo "Start time: $(date)"
echo ""
echo "Reference: ${REF_FILE}"
echo "Output: ${OUTPUT_FILE}"
echo ""

# Check if reference exists
if [ ! -f "${REF_FILE}" ]; then
    echo "ERROR: Reference h5ad not found at ${REF_FILE}"
    exit 1
fi

# Environment
CITEGEIST_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"

# Activate CITEgeist_env (has scanpy)
echo "Activating CITEgeist_env..."
export PATH="${CITEGEIST_ENV}/bin:$PATH"
source activate "${CITEGEIST_ENV}"

# Generate markers
python "${SRC_DIR}/generate_marker_genes.py" \
    --reference "${REF_FILE}" \
    --output "${OUTPUT_FILE}" \
    --cell-type-col "cell_type" \
    --n-markers 50

echo ""
echo "=============================================="
echo "Marker Generation Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Generated marker file:"
ls -la "${OUTPUT_FILE}"
echo ""
echo "Preview:"
head -20 "${OUTPUT_FILE}"
