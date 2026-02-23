#!/bin/bash
#SBATCH --job-name=hybrid_gpu
#SBATCH --output=logs/hybrid_gpu_%A_%a.out
#SBATCH --error=logs/hybrid_gpu_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Hybrid approach with GPU Cellpose:
# 1. Run Cellpose on GPU (much faster - minutes vs hours)
# 2. Run continuous model (CLR preprocessing, QP optimization)
# 3. Discretize using nuclei counts
# Achieves ~94% of continuous model quality with integer cell counts

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_cellpose \
    --use-gpu \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0
