#!/bin/bash
#SBATCH --cluster=gpu
#SBATCH --job-name=pc_mil
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/pc_mil_%j.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/pc_mil_%j.err
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

echo "=== PC-MIL Benchmark ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
date

# Purge default modules to avoid libstdc++ conflicts
module purge 2>/dev/null || true

# Activate environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Ensure conda's libstdc++ is found before system /lib64
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

set -u  # Re-enable nounset after conda activation

# Run benchmark
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_pc_mil.py \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/pc_mil \
    --features-dir Benchmarking/xenium_benchmarking/CITEgeist/output/vit_features \
    --n-epochs 200 \
    --device cuda \
    ${FOLD:+--fold $FOLD}

echo "=== Done ==="
date
