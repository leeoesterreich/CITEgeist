#!/bin/bash
#SBATCH --job-name=extract_eval
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm/logs/extract_eval_%j.out
#SBATCH --error=slurm/logs/extract_eval_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist

source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Clean LD_LIBRARY_PATH to avoid GLIBC conflicts
export LD_LIBRARY_PATH=""

echo "=== Step 1: Extract embeddings from masked SimCLR ==="
python src/extract_ssl_embeddings.py \
    --checkpoint output/ssl_models/simclr_masked_l40s/simclr_final.pt \
    --model-type simclr \
    --patches-dir output/patches_v3_masked \
    --output-dir output/ssl_embeddings_masked_l40s

echo "=== Step 2: Evaluate constrained assignment ==="
python src/analyze_distinct_types.py \
    --embeddings output/ssl_embeddings_masked_l40s/simclr_embeddings.npy \
    --max-spots 3000

echo "=== Done ==="
