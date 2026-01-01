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
# Generate GEX Ground Truth for Xenium Benchmarking
# ============================================================================
# Run this BEFORE the main benchmarking pipeline
# Creates per-cell-type gene expression ground truth files
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"

mkdir -p "${BENCH_DIR}/slurm/slurm_log"

echo "============================================"
echo "Generate GEX Ground Truth"
echo "============================================"
echo "Start time: $(date)"
echo ""

# Activate environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

# Run the GEX ground truth generation script
python "${BENCH_DIR}/src/generate_gex_ground_truth.py"

echo ""
echo "============================================"
echo "GEX Ground Truth Generation Complete"
echo "============================================"
echo "End time: $(date)"

# List generated files
echo ""
echo "Generated files:"
ls -la "${BENCH_DIR}/data/ground_truth_gex/"
