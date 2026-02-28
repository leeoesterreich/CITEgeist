#!/bin/bash
#SBATCH --job-name=eval_v2
#SBATCH --output=slurm/logs/eval_v2_%j.out
#SBATCH --error=slurm/logs/eval_v2_%j.err
#SBATCH --time=1:00:00
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Evaluate v2 DirectSoftmaxModel (with global normalization + size features)
set -e

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load cuda/12.1.1

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/vae_sinkhorn_v2"
PATCHES_V2="Benchmarking/xenium_benchmarking/CITEgeist/output/patches_v2"

echo "============================================"
echo "Evaluating v2 Model (global norm + size features)"
echo "============================================"
echo "VAE: ${OUTPUT_DIR}/vae/vae_final.pt"
echo "Prototypes: ${OUTPUT_DIR}/prototypes/prototypes_final.pt"
echo "============================================"

python Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_direct_softmax.py \
    --vae-checkpoint "${OUTPUT_DIR}/vae/vae_final.pt" \
    --model-checkpoint "${OUTPUT_DIR}/prototypes/prototypes_final.pt" \
    --patches-dir "${OUTPUT_DIR}/spot_patches" \
    --proportions-file "${OUTPUT_DIR}/proportions_for_training.csv" \
    --nucleus-features-dir "${PATCHES_V2}" \
    --output-dir "${OUTPUT_DIR}/eval_v2" \
    --regions 0 1 2 3 4 \
    --device cuda

echo ""
echo "============================================"
echo "Evaluation complete!"
echo "============================================"
echo "Results: ${OUTPUT_DIR}/eval_v2/evaluation_results.json"
