#!/bin/bash
#SBATCH --job-name=visiumhd_benchmark
#SBATCH --output=../logs/benchmark_%j.out
#SBATCH --error=../logs/benchmark_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Visium HD H&E Morphology Benchmark
# ==================================
#
# Usage:
#   sbatch sbatch_benchmark.sh
#
# Or with custom sample:
#   SAMPLE=TP12-880 sbatch sbatch_benchmark.sh
#
# To skip certain steps (resume from checkpoint):
#   SKIP_FLAGS="--skip-segmentation --skip-features" sbatch sbatch_benchmark.sh

set -e

# Allow TensorFlow to grow GPU memory (needed by StarDist/csbdeep)
export TF_FORCE_GPU_ALLOW_GROWTH=true

# Load modules
module load gurobi/12.0.3

# Activate environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Set paths
REPO_ROOT=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
PILC_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project
ENACT_DIR=/ihome/alee/alc376/alc376_bgfs/enact-pipeline/output_files

# Sample to process (override via environment variable)
SAMPLE=${SAMPLE:-TP08-2202}

# Skip flags (override via environment variable)
SKIP_FLAGS=${SKIP_FLAGS:-""}

# Map sample to paths
case $SAMPLE in
    TP08-2202)
        H5AD=$PILC_DIR/enact_adatas/TP08-2202_canonical_annotated.h5ad
        WSI=$ENACT_DIR/tp08-2202-pilc/tmap/wsi.tif
        ;;
    TP12-880)
        H5AD=$PILC_DIR/enact_adatas/TP12-880_canonical_annotated.h5ad
        WSI=$ENACT_DIR/tp12-880-pilc/tmap/wsi.tif
        ;;
    TP15-M509)
        H5AD=$PILC_DIR/enact_adatas/TP15-M509_canonical_annotated.h5ad
        WSI=$ENACT_DIR/tp15-m509-pilc/tmap/wsi.tif
        ;;
    *)
        echo "Unknown sample: $SAMPLE"
        echo "Available samples: TP08-2202, TP12-880, TP15-M509"
        exit 1
        ;;
esac

# Verify input files exist
if [ ! -f "$H5AD" ]; then
    echo "ERROR: H5AD file not found: $H5AD"
    exit 1
fi

if [ ! -f "$WSI" ]; then
    echo "ERROR: WSI file not found: $WSI"
    exit 1
fi

OUTPUT=$REPO_ROOT/Benchmarking/visiumhd_benchmarking/results/$SAMPLE

# ViT model settings
# Options: vit_small_patch16_224 (384-dim), vit_base_patch16_224 (768-dim), vit_large_patch16_224 (1024-dim)
VIT_MODEL=${VIT_MODEL:-vit_small_patch16_224}

# UNI weights (optional - requires Hugging Face download)
# VIT_WEIGHTS=$REPO_ROOT/models/uni_weights.pth
VIT_WEIGHTS=${VIT_WEIGHTS:-""}

# Training parameters
N_EPOCHS=${N_EPOCHS:-100}
BATCH_SIZE=${BATCH_SIZE:-32}
LR=${LR:-1e-4}

# Create output and log directories
mkdir -p $OUTPUT
mkdir -p $REPO_ROOT/Benchmarking/visiumhd_benchmarking/logs

# Log configuration
echo "=================================================="
echo "Visium HD H&E Morphology Benchmark"
echo "=================================================="
echo "Sample:      $SAMPLE"
echo "H5AD:        $H5AD"
echo "WSI:         $WSI"
echo "Output:      $OUTPUT"
echo "ViT Model:   $VIT_MODEL"
echo "ViT Weights: ${VIT_WEIGHTS:-"(ImageNet pretrained)"}"
echo "Epochs:      $N_EPOCHS"
echo "Batch Size:  $BATCH_SIZE"
echo "Skip Flags:  ${SKIP_FLAGS:-"(none)"}"
echo "=================================================="
echo ""

# Build command
CMD="python run_benchmark.py \
    --sample $H5AD \
    --wsi $WSI \
    --output $OUTPUT \
    --vit-model $VIT_MODEL \
    --device cuda \
    --n-epochs $N_EPOCHS \
    --batch-size $BATCH_SIZE \
    --lr $LR"

# Add custom weights if specified
if [ -n "$VIT_WEIGHTS" ] && [ -f "$VIT_WEIGHTS" ]; then
    CMD="$CMD --vit-weights $VIT_WEIGHTS"
fi

# Add skip flags
if [ -n "$SKIP_FLAGS" ]; then
    CMD="$CMD $SKIP_FLAGS"
fi

# Run benchmark
cd $REPO_ROOT/Benchmarking/visiumhd_benchmarking/src

echo "Running: $CMD"
echo ""

eval $CMD

echo ""
echo "=================================================="
echo "Benchmark complete for $SAMPLE"
echo "Results saved to: $OUTPUT"
echo "=================================================="
