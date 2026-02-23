#!/bin/bash
#SBATCH --job-name=cellpose_v2
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/cellpose_v2_%A_%a.out
#SBATCH --error=logs/cellpose_v2_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Discrete cellpose benchmark with correct coordinate conversion (v2)
# Fix: use pixel_bounds instead of padding_microns for coordinate transformation

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

REGION_ID=$SLURM_ARRAY_TASK_ID
OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose_v2

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py \
    --region $REGION_ID \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20

echo "Completed region $REGION_ID"
