#!/bin/bash
#SBATCH --job-name=umap_validation
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/slurm/slurm_log/validate_umap_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/slurm/slurm_log/validate_umap_%j.err

# UMAP Validation of Xenium 10-Cell-Type Ground Truth
# This script validates the RNA k-means clustering ground truth by
# generating UMAP visualizations and computing clustering metrics.
#
# Expected runtime: 2-4 hours for 465K cells
# Output: figures/ directory with PNGs and JSON metrics

# Use absolute paths
PROJECT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium"
CONDA_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env"

echo "=============================================="
echo "UMAP VALIDATION OF XENIUM GROUND TRUTH"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start time: $(date)"
echo "=============================================="

# Change to project directory
cd "$PROJECT_DIR"
echo "Working directory: $(pwd)"

# Activate conda environment using full path
source /ihome/alee/alc376/.bashrc
source /ihome/alee/alc376/miniforge3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

# Verify Python and key packages
echo ""
echo "Environment check:"
which python
python --version
python -c "import scanpy; print(f'scanpy: {scanpy.__version__}')"
python -c "import sklearn; print(f'sklearn: {sklearn.__version__}')"
echo ""

# Run validation script with absolute path
echo "Running UMAP validation..."
python "${PROJECT_DIR}/src/validate_ground_truth_umap.py"

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=============================================="
    echo "VALIDATION COMPLETE"
    echo "=============================================="
    echo "Output files in: figures/"
    ls -la figures/
else
    echo ""
    echo "=============================================="
    echo "VALIDATION FAILED"
    echo "=============================================="
    exit 1
fi

echo ""
echo "End time: $(date)"
