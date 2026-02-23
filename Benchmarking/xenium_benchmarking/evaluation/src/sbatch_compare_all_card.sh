#!/bin/bash
#SBATCH --job-name=eval_card
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src/eval_card_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src/eval_card_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Run method comparison evaluation including CARD

set -e

echo "=============================================="
echo "Running Method Comparison Evaluation (with CARD)"
echo "=============================================="
echo "Start time: $(date)"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src

python compare_all_methods.py

echo ""
echo "=============================================="
echo "Evaluation Complete!"
echo "=============================================="
echo "End time: $(date)"
