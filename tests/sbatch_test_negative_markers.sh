#!/bin/bash
#SBATCH --job-name=test_neg_markers
#SBATCH --output=slurm_log/test_neg_markers_%j.out
#SBATCH --error=slurm_log/test_neg_markers_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Load required modules
module load gurobi/12.0.3

# Activate conda environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Change to test directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/tests

# Run the test
echo "Running negative marker tests..."
python test_negative_markers.py

echo "Test complete!"
