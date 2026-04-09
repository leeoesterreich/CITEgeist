#!/bin/bash
#SBATCH --job-name=test_joint_real
#SBATCH --output=tests/slurm_log/test_joint_real_%j.out
#SBATCH --error=tests/slurm_log/test_joint_real_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Test discover_joint_programs on real HCC22-088-P1-S1 data

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Activate conda environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Starting test at $(date)"
echo "Python: $(which python)"
echo "Working directory: $(pwd)"

# Run the test with unbuffered output
export PYTHONUNBUFFERED=1
python -u tests/test_joint_on_real_data.py

echo "Test completed at $(date)"
