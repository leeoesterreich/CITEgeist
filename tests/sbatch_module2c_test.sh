#!/bin/bash
#SBATCH --job-name=module2c_test
#SBATCH --output=benchmarking_logs/module2c_test_%j.log
#SBATCH --error=benchmarking_logs/module2c_test_%j.log
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=8
#SBATCH --partition=HTC

# Module 2c: Reconstruction-Based Profile Selection Test
# Tests elbow method for selecting optimal number of profiles

# Activate conda environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory: $(pwd)"

# Ensure log directory exists
mkdir -p benchmarking_logs

echo "Starting Module 2c Profile Selection Test"
echo "Date: $(date)"
echo "========================================"

python tests/test_module2c_profile_selection.py

echo "========================================"
echo "Test completed: $(date)"
