#!/bin/bash
#SBATCH --job-name=Module4_test
#SBATCH --output=benchmarking_logs/Module4_test_%j.log
#SBATCH --error=benchmarking_logs/Module4_test_%j.log
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc

# Module 4 test - run on a single sample to verify protein-anchored program discovery

# Activate conda environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Load Gurobi module (sets GRB_LICENSE_FILE env var automatically)
module load gurobi/12.0.3
echo "Loaded gurobi module"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Create output directory
mkdir -p output
mkdir -p benchmarking_logs

# Data folder
DATA_FOLDER="/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files"

# Test on HCC22-088-P1-S2 (high-quality biopsy sample)
SAMPLE="HCC22-088-P1-S2/outs"

echo "Testing Module 4 on sample: $SAMPLE"
echo "Using auto-profiles (default)"

# Run with auto-profiles (default now) - unbuffered output for real-time logging
PYTHONUNBUFFERED=1 python examples/compute_sample.py \
    --path "$DATA_FOLDER/$SAMPLE" \
    --output_folder output \
    --radius 400 \
    --min_counts 100 \
    --top-k 3 \
    --discovery-seed 1234 \
    --n-permutations 999 \
    --fdr-threshold 0.05

echo "Test complete"
