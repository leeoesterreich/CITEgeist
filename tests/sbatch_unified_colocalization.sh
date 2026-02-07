#!/bin/bash
#SBATCH --job-name=unified_coloc_test
#SBATCH --output=slurm_log/unified_coloc_%j.out
#SBATCH --error=slurm_log/unified_coloc_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Test unified colocalization on pILC single-cell data
#
# Compares:
# 1. Original approach (bivariate Moran's I) - fails on single-cell
# 2. Unified approach (spatial pattern similarity) - should work
#
# Usage:
#   cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests
#   sbatch sbatch_unified_colocalization.sh

echo "=============================================="
echo "Unified Colocalization Test"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo ""

# Create log directory
mkdir -p slurm_log

# Navigate to tests directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests

# Activate environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Run comparison test
python test_unified_colocalization.py \
    --sample TP08-2202 \
    --max-cells 10000 \
    --n-markers 40 \
    --n-perm 299 \
    --output test_results/unified_colocalization/TP08-2202

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Test Complete"
echo "=============================================="
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

exit $EXIT_CODE
