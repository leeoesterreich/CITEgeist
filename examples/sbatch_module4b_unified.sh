#!/bin/bash
#SBATCH --job-name=module4b_bivariate
#SBATCH --output=slurm_log/module4b_unified_%A_%a.out
#SBATCH --error=slurm_log/module4b_unified_%A_%a.err
#SBATCH --array=0-11
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Module 4b: Bivariate Program Relationship Analysis
# Uses NB-derived Module 4 outputs

mkdir -p slurm_log

# Deduplicated 12 patient samples
SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Module 4b: Bivariate Program Relationships"
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

OUTPUT_DIR="output/module4b_nb"
mkdir -p ${OUTPUT_DIR}

# Run Module 4b
$PYTHON examples/run_module4b_bivariate.py \
    --sample ${SAMPLE} \
    --output-dir ${OUTPUT_DIR} \
    --n-permutations 199 \
    --fdr-threshold 0.05 \
    --min-bivariate-i 0.1 \
    --neighbor-k 8 \
    --seed 42

echo "=========================================="
echo "Completed: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
