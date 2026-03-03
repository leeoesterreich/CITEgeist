#!/bin/bash
#SBATCH --job-name=eval_ssl
#SBATCH --output=slurm/logs/eval_ssl_%j.out
#SBATCH --error=slurm/logs/eval_ssl_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -e

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
BENCHMARK_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"

cd "${BENCHMARK_DIR}"

# Extract MAE embeddings
python src/extract_ssl_embeddings.py \
    --checkpoint output/mae_ssl/mae_final.pt \
    --model-type mae \
    --patches-dir output/vae_masked/patches_combined \
    --output-dir output/ssl_embeddings

# Extract DINO embeddings
python src/extract_ssl_embeddings.py \
    --checkpoint output/dino_ssl/dino_final.pt \
    --model-type dino \
    --patches-dir output/vae_masked/patches_combined \
    --output-dir output/ssl_embeddings

# Evaluate all
python src/evaluate_ssl_embeddings.py \
    --embeddings-dir output/ssl_embeddings \
    --output-dir output/ssl_evaluation \
    --models mae dino

echo "Evaluation complete"
