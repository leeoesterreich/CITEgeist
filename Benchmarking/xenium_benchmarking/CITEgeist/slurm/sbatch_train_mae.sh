#!/bin/bash
#SBATCH --job-name=train_mae
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=logs/train_mae_%j.out
#SBATCH --error=logs/train_mae_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# MAE Training Script for Nucleus Patch Representation Learning
# ============================================================================
# Trains a Masked Autoencoder (MAE) on 2-channel nucleus patches (DAPI + boundary)
# for self-supervised representation learning.
#
# Architecture:
#   - Encoder: ViT-Small (384-dim, 12 layers, 6 heads)
#   - Decoder: Lightweight transformer (192-dim, 4 layers, 3 heads)
#   - Masking: 75% of patches masked during training
#
# Expected runtime: ~8-12 hours for 200 epochs on ~100k patches
# ============================================================================

set -e

# Create logs directory if not exists
mkdir -p logs

# Load modules and activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Change to project root
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Define paths (can be overridden via environment variables)
PATCHES_DIR="${PATCHES_DIR:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output/vae_masked/patches_combined}"
OUTPUT_DIR="${OUTPUT_DIR:-/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output/mae_ssl}"

# Training hyperparameters (can be overridden via environment variables)
EPOCHS="${EPOCHS:-200}"
BATCH_SIZE="${BATCH_SIZE:-256}"
LR="${LR:-1.5e-4}"
WEIGHT_DECAY="${WEIGHT_DECAY:-0.05}"
WARMUP_EPOCHS="${WARMUP_EPOCHS:-10}"
MASK_RATIO="${MASK_RATIO:-0.75}"
CHECKPOINT_INTERVAL="${CHECKPOINT_INTERVAL:-50}"
NUM_WORKERS="${NUM_WORKERS:-8}"

echo "========================================"
echo "MAE Training Script"
echo "========================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo ""
echo "Patches directory: $PATCHES_DIR"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Hyperparameters:"
echo "  Epochs: $EPOCHS"
echo "  Batch size: $BATCH_SIZE"
echo "  Learning rate: $LR"
echo "  Weight decay: $WEIGHT_DECAY"
echo "  Warmup epochs: $WARMUP_EPOCHS"
echo "  Mask ratio: $MASK_RATIO"
echo "  Checkpoint interval: $CHECKPOINT_INTERVAL"
echo "  Num workers: $NUM_WORKERS"
echo "========================================"

# Check if patches directory exists
if [ ! -d "$PATCHES_DIR" ]; then
    echo "ERROR: Patches directory does not exist: $PATCHES_DIR"
    echo "Please run prepare_patches.py first to extract nucleus patches."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run MAE training
python Benchmarking/xenium_benchmarking/CITEgeist/src/train_mae.py \
    --patches-dir "$PATCHES_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --epochs "$EPOCHS" \
    --batch-size "$BATCH_SIZE" \
    --lr "$LR" \
    --weight-decay "$WEIGHT_DECAY" \
    --warmup-epochs "$WARMUP_EPOCHS" \
    --mask-ratio "$MASK_RATIO" \
    --device cuda \
    --checkpoint-interval "$CHECKPOINT_INTERVAL" \
    --num-workers "$NUM_WORKERS"

echo "========================================"
echo "MAE training completed successfully"
echo "========================================"
echo "Model saved to: $OUTPUT_DIR/mae_final.pt"
echo "Encoder weights: $OUTPUT_DIR/mae_final.pt (encoder_state_dict key)"
echo "Training history: $OUTPUT_DIR/training_history.json"
echo "Training config: $OUTPUT_DIR/training_config.json"
echo "========================================"
