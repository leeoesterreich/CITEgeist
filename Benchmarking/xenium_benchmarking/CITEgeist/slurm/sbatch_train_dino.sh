#!/bin/bash
#SBATCH --job-name=train_dino
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=slurm/logs/train_dino_%j.out
#SBATCH --error=slurm/logs/train_dino_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# train_dino SLURM script

set -e

# Clean LD_LIBRARY_PATH to avoid GLIBC conflicts with system libraries
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v "ondemand-jupyter" | tr '\n' ':' | sed 's/:$//')

echo "============================================"
echo "DINO Training Job"
echo "Started: $(date)"
echo "============================================"

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCHMARK_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
PATCHES_DIR="${BENCHMARK_DIR}/output/vae_masked/patches_combined"
OUTPUT_DIR="${BENCHMARK_DIR}/output/dino_ssl"

cd "${BENCHMARK_DIR}"
mkdir -p slurm/logs

python src/train_dino.py \
    --patches-dir "${PATCHES_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --epochs 300 \
    --batch-size 128 \
    --lr 5e-4 \
    --device cuda \
    --checkpoint-interval 50 \
    --num-workers 8

echo "Completed: $(date)"
