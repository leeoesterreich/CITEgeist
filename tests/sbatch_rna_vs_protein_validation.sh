#!/bin/bash
#SBATCH --job-name=rna_protein_validation
#SBATCH --output=slurm_log/rna_protein_validation_%j.out
#SBATCH --error=slurm_log/rna_protein_validation_%j.err
#SBATCH --time=02:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# RNA vs Protein Profile Discovery Validation
#
# This script validates that RNA-based profile discovery can recover
# similar structure to protein-based discovery using real Xenium data.
#
# Usage:
#   sbatch sbatch_rna_vs_protein_validation.sh
#
# Expected runtime: ~30-60 minutes

echo "=============================================="
echo "RNA vs Protein Profile Discovery Validation"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Start time: $(date)"
echo ""

# Navigate to repo root
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Create output directories
mkdir -p tests/slurm_log
mkdir -p tests/test_results/rna_vs_protein_validation

# Activate conda environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Run validation
echo "Running RNA vs Protein validation..."
echo "=============================================="

python tests/test_rna_vs_protein_xenium.py \
    --max-cells 20000 \
    --n-permutations 199 \
    --seed 42

EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Validation Complete"
echo "=============================================="
echo "Exit code: $EXIT_CODE"
echo "End time: $(date)"

if [ $EXIT_CODE -eq 0 ]; then
    echo "STATUS: VALIDATION PASSED"
else
    echo "STATUS: VALIDATION FAILED"
fi

# Show results summary if exists
SUMMARY_FILE="tests/test_results/rna_vs_protein_validation/validation_summary.json"
if [ -f "$SUMMARY_FILE" ]; then
    echo ""
    echo "Results summary:"
    cat "$SUMMARY_FILE"
fi

exit $EXIT_CODE
