#!/bin/bash
#SBATCH --job-name=extract_vit
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/extract_vit_%j.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/extract_vit_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -euo pipefail

echo "=== ViT Feature Extraction ==="
date

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/extract_vit_features.py \
    --patches-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vae_masked/patches_combined \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vit_features \
    --device cuda

echo "=== Done ==="
date
