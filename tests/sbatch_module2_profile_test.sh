#!/bin/bash
#SBATCH --job-name=module2_profile_test
#SBATCH --output=benchmarking_logs/module2_profile_test_%j.log
#SBATCH --error=benchmarking_logs/module2_profile_test_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=8
#SBATCH --partition=HTC

# Module 2 Profile Discovery Test
# Tests gap-based lineage splitting and dynamic tree cutting
# on both simulated (9 GT cell types) and real patient data

# Activate conda environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory: $(pwd)"

# Ensure log directory exists
mkdir -p benchmarking_logs

echo "Starting Module 2 Profile Discovery Test"
echo "Date: $(date)"
echo "========================================"

python tests/test_module2_profile_discovery.py

echo "========================================"
echo "Test completed: $(date)"
