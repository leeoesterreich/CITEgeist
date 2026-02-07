#!/bin/bash
#SBATCH --job-name=module3_unified
#SBATCH --output=slurm_log/module3_unified_%A_%a.out
#SBATCH --error=slurm_log/module3_unified_%A_%a.err
#SBATCH --array=0-13
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Module 3 with Unified Profile (10 cell types)
# Runs cell proportion estimation + gene expression deconvolution
# Using unified profile derived from Module 1-2 discovery

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
echo "Module 3 with Unified Profile"
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

# Change to project directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Load Gurobi module for license
module load gurobi/12.0.3

# Use direct python path
PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

# Output directory
OUTPUT_DIR="output/module3_unified"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/checkpoints

# Run Module 3 with unified profile
$PYTHON examples/run_module3_unified.py \
    --sample ${SAMPLE} \
    --output-dir ${OUTPUT_DIR}

echo "=========================================="
echo "Completed: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
