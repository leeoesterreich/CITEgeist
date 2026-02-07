#!/bin/bash
#SBATCH --job-name=unified_both_res
#SBATCH --output=slurm_log/unified_both_%j.out
#SBATCH --error=slurm_log/unified_both_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Test unified colocalization on BOTH resolutions
#
# Compares scoring methods:
# - weighted (legacy)
# - geometric (sqrt of product)
# - min (bottleneck)
# - harmonic (2ab/(a+b))
#
# On:
# - Single-cell: pILC Visium HD
# - Spot-level: Xenium pseudo-Visium

echo "=============================================="
echo "Unified Colocalization: Both Resolutions"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests
mkdir -p slurm_log

source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Python: $(which python)"
echo ""

python test_unified_both_resolutions.py \
    --max-cells 10000 \
    --max-spots 5000 \
    --output test_results/unified_both_resolutions

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Complete"
echo "=============================================="
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

exit $EXIT_CODE
