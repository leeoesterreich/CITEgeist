#!/bin/bash
#SBATCH --job-name=process_ref_p7
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/slurm_log/process_ref_protein7_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/slurm_log/process_ref_protein7_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Process GSE156632 reference data with 7 cell types matching protein GT
# Output: processed_protein7/ with cell2location, tangram, rctd, seurat subdirs

set -e

BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/reference_data/GSE156632"
SRC_DIR="${BASE_DIR}/src"
INPUT_DIR="${BASE_DIR}"
OUTPUT_DIR="${BASE_DIR}/processed_protein7"

mkdir -p "${OUTPUT_DIR}"

echo "=============================================="
echo "Processing GSE156632 Reference (7 Protein-GT Types)"
echo "=============================================="
echo "Start time: $(date)"
echo "Output directory: ${OUTPUT_DIR}"

# Environment
ENV_BASE="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs"
CELL2LOC_ENV="${ENV_BASE}/cell2location_env"
CITEGEIST_ENV="${ENV_BASE}/CITEgeist_env"

if [ -d "${CELL2LOC_ENV}" ]; then
    echo "Activating cell2location_env..."
    export PATH="${CELL2LOC_ENV}/bin:$PATH"
    source activate "${CELL2LOC_ENV}"
elif [ -d "${CITEGEIST_ENV}" ]; then
    echo "Using CITEgeist_env..."
    export PATH="${CITEGEIST_ENV}/bin:$PATH"
    source activate "${CITEGEIST_ENV}"
else
    echo "ERROR: No suitable conda environment found"
    exit 1
fi

echo "Python: $(which python)"
echo ""

python "${SRC_DIR}/process_reference.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --min-genes 200 \
    --max-genes 5000 \
    --max-mt-pct 20.0

echo ""
echo "=============================================="
echo "Processing Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output:"
ls -la "${OUTPUT_DIR}/"
for subdir in cell2location tangram rctd seurat; do
    echo "${subdir}:"
    ls -la "${OUTPUT_DIR}/${subdir}/" 2>/dev/null || echo "  (not created)"
done
