#!/bin/bash
#SBATCH --job-name=rna_marker_tests
#SBATCH --output=slurm_log/rna_marker_tests_%j.out
#SBATCH --error=slurm_log/rna_marker_tests_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# RNA Marker Selection Module - Comprehensive Tests
#
# This script runs all unit and integration tests for the RNA marker
# selection module, which enables CITEgeist to work with RNA-only data.
#
# Usage:
#   sbatch sbatch_rna_marker_tests.sh
#
# Expected runtime: ~30-60 minutes depending on cluster load

echo "=============================================="
echo "RNA Marker Selection Tests"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo ""

# Navigate to repo root
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Create slurm_log directory if needed
mkdir -p tests/slurm_log

# Activate conda environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Run tests
echo "Running comprehensive tests..."
echo "=============================================="

python tests/test_rna_marker_selection.py

TEST_EXIT_CODE=$?

echo ""
echo "=============================================="
echo "Test Summary"
echo "=============================================="
echo "Exit code: $TEST_EXIT_CODE"
echo "End time: $(date)"

if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "STATUS: ALL TESTS PASSED"
else
    echo "STATUS: SOME TESTS FAILED"
fi

exit $TEST_EXIT_CODE
