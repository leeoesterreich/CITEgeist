#!/bin/bash
#SBATCH --job-name=module4_unified
#SBATCH --output=slurm_log/module4_unified_%A_%a.out
#SBATCH --error=slurm_log/module4_unified_%A_%a.err
#SBATCH --array=0-13
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Module 4: Anchored + Joint Program Discovery
# Uses Module 3 unified profile outputs (10 cell types)

mkdir -p slurm_log

# All 14 patient samples
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

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Module 4: Program Discovery"
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Load Gurobi (not strictly needed for Module 4 but doesn't hurt)
module load gurobi/12.0.3

PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

OUTPUT_DIR="output/module4_unified"
mkdir -p ${OUTPUT_DIR}

# Run Module 4 with both anchored and joint discovery
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
