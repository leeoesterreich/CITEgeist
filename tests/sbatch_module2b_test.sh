#!/bin/bash
#SBATCH --job-name=module2b_test
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests/slurm_log/module2b_test_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests/slurm_log/module2b_test_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=smp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo "Starting Module 2b parameter test at $(date)"

# Activate conda environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run test
python tests/test_module2b_relaxed.py

echo "Finished at $(date)"
