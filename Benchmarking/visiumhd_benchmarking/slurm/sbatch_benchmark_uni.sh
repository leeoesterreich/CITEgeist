#!/bin/bash
#SBATCH --job-name=visiumhd_uni
#SBATCH --output=../logs/benchmark_uni_%j.out
#SBATCH --error=../logs/benchmark_uni_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Visium HD H&E Benchmark with UNI Foundation Model
# ==================================================
#
# UNI is a pathology-specific ViT-L/16 trained on 100M+ histology patches.
# Expected improvement: +15-25% over ImageNet ViT.
#
# Prerequisites:
#   1. HuggingFace account with institutional email
#   2. Accept license at https://huggingface.co/MahmoodLab/UNI
#   3. Login: huggingface-cli login
#
# Usage:
#   sbatch sbatch_benchmark_uni.sh

set -e

# Clean up LD_LIBRARY_PATH to avoid GLIBC conflicts
export LD_LIBRARY_PATH=""

# Allow TensorFlow to grow GPU memory (needed by StarDist/csbdeep)
export TF_FORCE_GPU_ALLOW_GROWTH=true

# Set cache directories
export HF_HOME=/ihome/alee/alc376/.cache/huggingface
export TORCH_HOME=/ihome/alee/alc376/.cache/torch

# Load modules
module load gurobi/12.0.3

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Set paths
REPO_ROOT=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
PILC_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/pILC_project
ENACT_DIR=/ihome/alee/alc376/alc376_bgfs/enact-pipeline/output_files

# Sample
SAMPLE=${SAMPLE:-TP08-2202}

# UNI model configuration
VIT_MODEL="vit_large_patch16_224"
UNI_WEIGHTS=$REPO_ROOT/UNI_model/pytorch_model.bin

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
        exit 1
        ;;
esac

# Verify UNI weights exist
if [ ! -f "$UNI_WEIGHTS" ]; then
    echo "ERROR: UNI weights not found at: $UNI_WEIGHTS"
    exit 1
fi

# Verify inputs
if [ ! -f "$H5AD" ]; then
    echo "ERROR: H5AD file not found: $H5AD"
    exit 1
fi

if [ ! -f "$WSI" ]; then
    echo "ERROR: WSI file not found: $WSI"
    exit 1
fi

OUTPUT=$REPO_ROOT/Benchmarking/visiumhd_benchmarking/results/${SAMPLE}_UNI

# Training parameters
N_EPOCHS=${N_EPOCHS:-100}
BATCH_SIZE=${BATCH_SIZE:-16}  # Smaller batch for larger model
LR=${LR:-1e-4}

# Create directories
mkdir -p $OUTPUT
mkdir -p $REPO_ROOT/Benchmarking/visiumhd_benchmarking/logs

# Log configuration
echo "=================================================="
echo "Visium HD H&E Benchmark with UNI Foundation Model"
echo "=================================================="
echo "Sample:      $SAMPLE"
echo "H5AD:        $H5AD"
echo "WSI:         $WSI"
echo "Output:      $OUTPUT"
echo "ViT Model:   $VIT_MODEL (UNI pathology foundation)"
echo "UNI Weights: $UNI_WEIGHTS"
echo "Embed Dim:   1024"
echo "Epochs:      $N_EPOCHS"
echo "Batch Size:  $BATCH_SIZE"
echo "=================================================="
echo ""

# Run benchmark
cd $REPO_ROOT/Benchmarking/visiumhd_benchmarking/src

python run_benchmark.py \
    --sample $H5AD \
    --wsi $WSI \
    --output $OUTPUT \
    --vit-model $VIT_MODEL \
    --vit-weights $UNI_WEIGHTS \
    --device cuda \
    --n-epochs $N_EPOCHS \
    --batch-size $BATCH_SIZE \
    --lr $LR

echo ""
echo "=================================================="
echo "UNI Benchmark complete for $SAMPLE"
echo "Results saved to: $OUTPUT"
echo "=================================================="
