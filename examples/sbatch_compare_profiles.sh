#!/bin/bash
#SBATCH --job-name=compare_profiles
#SBATCH --output=slurm_log/compare_profiles_%j.out
#SBATCH --error=slurm_log/compare_profiles_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Phase 2 Step 2: Compare discovered vs curated profiles
# Aggregates Module 1-2 results and generates comparison reports

echo "=========================================="
echo "Profile Comparison: Discovered vs Curated"
echo "Date: $(date)"
echo "=========================================="

# Create log directory
mkdir -p slurm_log

# Change to project directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Use direct python path
PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

# Input: Module 1-2 discovery results
INPUT_DIR="output/module12_discovery"

# Output: Comparison reports
OUTPUT_DIR="output/profile_comparison"
mkdir -p ${OUTPUT_DIR}

echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Run comparison
$PYTHON examples/compare_profiles.py \
    --input-dir ${INPUT_DIR} \
    --output-dir ${OUTPUT_DIR}

echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="
echo ""
echo "Output files:"
ls -la ${OUTPUT_DIR}/
