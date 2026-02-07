#!/bin/bash
#SBATCH --job-name=perm5000_test
#SBATCH --output=benchmarking_logs/perm5000_test_%j.log
#SBATCH --error=benchmarking_logs/perm5000_test_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=8
#SBATCH --partition=HTC

# Test with 5000 permutations to check timing and FDR

source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory: $(pwd)"

mkdir -p benchmarking_logs

echo "Starting 5000 permutation test"
echo "Date: $(date)"
echo "========================================"

python tests/test_5000_perms.py

echo "========================================"
echo "Test completed: $(date)"
