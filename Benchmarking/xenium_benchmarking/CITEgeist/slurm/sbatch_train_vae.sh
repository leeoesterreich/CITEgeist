#!/bin/bash
#SBATCH --job-name=train_vae
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/train_vae_%j.out
#SBATCH --error=logs/train_vae_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# VAE Training for Nucleus Patch Representation Learning (Stage 1)
# This script trains a Variational Autoencoder on pre-extracted nucleus patches.
# The trained encoder will be used in Stage 2 for prototype learning.

set -e

# Create logs directory if not exists
mkdir -p logs

# Load modules and activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Change to project root
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Define paths (modify as needed)
PATCHES_DIR="${PATCHES_DIR:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_patches}"
OUTPUT_DIR="${OUTPUT_DIR:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_vae}"

echo "========================================"
echo "VAE Training Script"
echo "========================================"
echo "Patches directory: $PATCHES_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Device: cuda"
echo "========================================"

# Run VAE training
python -m CITEgeist.model.train_vae \
    --patches-dir "$PATCHES_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --in-channels 3 \
    --latent-dim 128 \
    --batch-size 64 \
    --epochs 100 \
    --lr 1e-4 \
    --beta 0.5 \
    --device cuda \
    --checkpoint-interval 10 \
    --num-workers 4

echo "========================================"
echo "VAE training completed successfully"
echo "Model saved to: $OUTPUT_DIR/vae_final.pt"
echo "========================================"
