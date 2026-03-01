#!/bin/bash
#SBATCH --job-name=two_stage
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/two_stage_%A_%a.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/two_stage_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load modules
module load gurobi/12.0.3

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Navigate to repo
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run benchmark for this region
REGION=${SLURM_ARRAY_TASK_ID}
OUTPUT_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/two_stage"

echo "=== Job Info ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Region: ${REGION}"
echo "Start time: $(date)"

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage.py \
    --region ${REGION} \
    --output-dir ${OUTPUT_DIR} \
    --device cuda \
    --n-epochs 50

echo "End time: $(date)"
