#!/bin/bash
#SBATCH --job-name=eval_morph
#SBATCH --output=logs/eval_morphology_%j.out
#SBATCH --error=logs/eval_morphology_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Morphology Assignment Evaluation
# Evaluates whether morphology features improve cell type assignment

set -e

echo "=========================================="
echo "Morphology Assignment Evaluation"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "=========================================="

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Paths
SC_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution"
GT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium"
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/results/morphology_evaluation"

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Run evaluation
python /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src/evaluate_morphology_assignment.py \
    --sc-dir ${SC_DIR} \
    --gt-dir ${GT_DIR} \
    --regions 0,1,2,4 \
    --output-dir ${OUTPUT_DIR}

echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
