#!/bin/bash
#SBATCH --job-name=mil_assign
#SBATCH --array=0-4
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/logs/mil_assign_%A_%a.out
#SBATCH --error=slurm/logs/mil_assign_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Clean LD_LIBRARY_PATH to avoid GLIBC/CUDA conflicts
export LD_LIBRARY_PATH=""

REGION=${SLURM_ARRAY_TASK_ID}

SIMCLR_CHECKPOINT="output/simclr_ssl/simclr_best.pt"
OUTPUT_DIR="output/mil_assignment"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist

python src/benchmark_mil_assignment.py \
    --region ${REGION} \
    --simclr-checkpoint ${SIMCLR_CHECKPOINT} \
    --output-dir ${OUTPUT_DIR} \
    --lambda-prior 1.0 \
    --n-epochs 100 \
    --device cuda \
    --batch-size 64
