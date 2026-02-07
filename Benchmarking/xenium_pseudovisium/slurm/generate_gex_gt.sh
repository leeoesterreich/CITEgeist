#!/bin/bash
#SBATCH --job-name=gen_gex_gt          # Job name
#SBATCH --output=slurm_log/%x_%j.out   # Standard output
#SBATCH --error=slurm_log/%x_%j.err    # Standard error
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=4              # CPU cores per task
#SBATCH --mem=128G                     # Memory (large for 465K cells)
#SBATCH --time=04:00:00                # Time limit
#SBATCH --cluster=htc                  # Cluster
#SBATCH --partition=htc                # Partition

# ============================================================================
# Generate GEX Ground Truth for Xenium Benchmarking (RNA-based)
# ============================================================================
# Run this BEFORE the main benchmarking pipeline
# Creates per-cell-type gene expression ground truth files
#
# Uses RNA-based cell type classification (recommended) to avoid circular
# logic with protein markers.
#
# Reference:
#   Zhao et al. (2025). "Benchmarking cell type annotation methods for 10x
#   Xenium spatial transcriptomics data." BMC Bioinformatics, 26(1), 25.
#   https://doi.org/10.1186/s12859-025-06044-0
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
OUTPUT_DIR="${PSEUDOVISIUM_DIR}/data_rna_gt"

mkdir -p "${PSEUDOVISIUM_DIR}/slurm/slurm_log"

echo "============================================"
echo "Generate GEX Ground Truth (RNA-based)"
echo "============================================"
echo "Start time: $(date)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Activate environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

# Run the GEX ground truth generation script with RNA-based cell types (default)
python "${PSEUDOVISIUM_DIR}/src/generate_gex_ground_truth.py" \
    --use-rna-types \
    --output-dir "${OUTPUT_DIR}"

echo ""
echo "============================================"
echo "GEX Ground Truth Generation Complete"
echo "============================================"
echo "End time: $(date)"

# List generated files
echo ""
echo "Generated files:"
ls -la "${OUTPUT_DIR}/ground_truth_gex/" 2>/dev/null || echo "No files generated yet"
for region_dir in "${OUTPUT_DIR}/ground_truth_gex/Xenium_region_"*; do
    if [ -d "$region_dir" ]; then
        echo ""
        echo "$(basename $region_dir):"
        ls -la "$region_dir"
    fi
done
