#!/bin/bash
#SBATCH --job-name=module12_rerun
#SBATCH --output=slurm_log/module12_rerun_%A_%a.out
#SBATCH --error=slurm_log/module12_rerun_%A_%a.err
#SBATCH --array=0-2
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Rerun the 3 timed-out samples with 4-hour limit

mkdir -p slurm_log

# Only the 3 samples that timed out
SAMPLES=(
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

OUTPUT_DIR="output/module12_discovery"
mkdir -p ${OUTPUT_DIR}

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
