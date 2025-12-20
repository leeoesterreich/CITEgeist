#!/bin/bash
#SBATCH --job-name=coverage_test
#SBATCH --output=benchmarking_logs/coverage_test_%j.log
#SBATCH --error=benchmarking_logs/coverage_test_%j.log
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=8
#SBATCH --partition=HTC

# Module 2b/2c Refactoring Test
# Tests two-support attachment (2b) and coverage-based selection (2c)

# Activate conda environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory: $(pwd)"

# Ensure log directory exists
mkdir -p benchmarking_logs

echo "Starting Coverage-Based Selection Test"
echo "Date: $(date)"
echo "========================================"

python tests/test_coverage_selection.py

echo "========================================"
echo "Test completed: $(date)"
