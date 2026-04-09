#!/bin/bash
#SBATCH --job-name=module12_discovery
#SBATCH --output=slurm_log/module12_discovery_%A_%a.out
#SBATCH --error=slurm_log/module12_discovery_%A_%a.err
#SBATCH --array=0-13
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Module 1-2 Discovery on all 14 patient samples
# Runs profile discovery pipeline to compare with curated profiles

# Create log directory
mkdir -p slurm_log

# Patient samples array (14 total)
SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

# Get current sample from array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

# Change to project directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Use direct python path (conda activation unreliable on compute nodes)
PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

# Output directory
OUTPUT_DIR="output/module12_discovery"
mkdir -p ${OUTPUT_DIR}

# Run Module 1-2 discovery
$PYTHON examples/run_module12_discovery.py \
    --sample ${SAMPLE} \
    --output-dir ${OUTPUT_DIR} \
    --top-k 3 \
    --n-permutations 999 \
    --fdr-threshold 0.05 \
    --variance-target 0.90 \
    --min-marginal-gain 0.005 \
    --seed 42

echo "=========================================="
echo "Completed: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
