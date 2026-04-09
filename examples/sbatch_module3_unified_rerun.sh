#!/bin/bash
#SBATCH --job-name=module3_rerun
#SBATCH --output=slurm_log/module3_rerun_%A_%a.out
#SBATCH --error=slurm_log/module3_rerun_%A_%a.err
#SBATCH --array=0-5
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Re-run Module 3 for 6 failed S2 samples (NaN coordinate fix applied)

mkdir -p slurm_log

# Failed samples (all S2 surgical samples)
SAMPLES=(
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S2"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Module 3 Rerun (NaN fix)"
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Load Gurobi module for license
module load gurobi/12.0.3

PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

OUTPUT_DIR="output/module3_unified"
mkdir -p ${OUTPUT_DIR}
mkdir -p ${OUTPUT_DIR}/checkpoints

$PYTHON examples/run_module3_unified.py \
    --sample ${SAMPLE} \
    --output-dir ${OUTPUT_DIR}

echo "=========================================="
echo "Completed: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
