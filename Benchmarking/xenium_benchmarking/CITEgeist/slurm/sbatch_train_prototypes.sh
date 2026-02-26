#!/bin/bash
#SBATCH --job-name=train_proto
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/train_proto_%j.out
#SBATCH --error=logs/train_proto_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Prototype Learning Training (Stage 2)
# This script trains projection heads and prototypes for cell type assignment
# using a pre-trained VAE encoder from Stage 1.

set -e

# Create logs directory if not exists
mkdir -p logs

# Load modules and activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Change to project root
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Define paths (modify as needed)
VAE_CHECKPOINT="${VAE_CHECKPOINT:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_vae/vae_final.pt}"
PATCHES_DIR="${PATCHES_DIR:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_patches}"
PROPORTIONS_FILE="${PROPORTIONS_FILE:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_proportions/proportions.csv}"
OUTPUT_DIR="${OUTPUT_DIR:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_prototypes}"

# Number of cell types (default 7 for Xenium benchmark)
N_TYPES="${N_TYPES:-7}"

echo "========================================"
echo "Prototype Learning Training Script"
echo "========================================"
echo "VAE checkpoint: $VAE_CHECKPOINT"
echo "Patches directory: $PATCHES_DIR"
echo "Proportions file: $PROPORTIONS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Number of cell types: $N_TYPES"
echo "Device: cuda"
echo "========================================"

# Run prototype training
python -m CITEgeist.model.train_prototypes \
    --vae-checkpoint "$VAE_CHECKPOINT" \
    --patches-dir "$PATCHES_DIR" \
    --proportions-file "$PROPORTIONS_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --n-types "$N_TYPES" \
    --latent-dim 128 \
    --projection-dim 32 \
    --epochs 50 \
    --lr 1e-3 \
    --sinkhorn-temp 0.1 \
    --sinkhorn-iters 50 \
    --device cuda \
    --checkpoint-interval 10 \
    --min-cells-per-spot 3

echo "========================================"
echo "Prototype training completed successfully"
echo "Model saved to: $OUTPUT_DIR/prototypes_final.pt"
echo "========================================"
