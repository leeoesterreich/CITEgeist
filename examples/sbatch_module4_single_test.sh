#!/bin/bash
#SBATCH --job-name=module4_single
#SBATCH --output=slurm_log/module4_single_%j.out
#SBATCH --error=slurm_log/module4_single_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Test Module 4 on one sample before running full array

echo "=========================================="
echo "Module 4 Test Run"
echo "Sample: HCC22-088-P1-S1"
echo "Date: $(date)"
echo "=========================================="

mkdir -p slurm_log

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

OUTPUT_DIR="output/module4_discovery"
mkdir -p ${OUTPUT_DIR}

$PYTHON examples/run_module4_discovery.py \
    --sample HCC22-088-P1-S1 \
    --output-dir ${OUTPUT_DIR} \
    --mode both \
    --k-anchored 5 \
    --k-joint 10 \
    --top-n-genes 50 \
    --seed 42

echo "=========================================="
echo "Test complete: $(date)"
echo "=========================================="
