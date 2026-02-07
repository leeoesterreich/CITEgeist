#!/bin/bash
#SBATCH --job-name=module4_discovery
#SBATCH --output=slurm_log/module4_discovery_%A_%a.out
#SBATCH --error=slurm_log/module4_discovery_%A_%a.err
#SBATCH --array=0-13
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Phase 3: Module 4 Program Discovery on all 14 patient samples
# Runs both anchored (per cell type) and joint (cross cell type) discovery

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
echo "Module 4 Program Discovery"
echo "Sample: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

# Change to project directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Use direct python path
PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

# Output directory
OUTPUT_DIR="output/module4_discovery"
mkdir -p ${OUTPUT_DIR}

# Run Module 4 discovery (both modes)
$PYTHON examples/run_module4_discovery.py \
    --sample ${SAMPLE} \
    --output-dir ${OUTPUT_DIR} \
    --mode both \
    --k-anchored 5 \
    --k-joint 10 \
    --top-n-genes 50 \
    --seed 42

echo "=========================================="
echo "Completed: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
