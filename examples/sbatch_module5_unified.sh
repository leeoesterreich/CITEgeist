#!/bin/bash
#SBATCH --job-name=module5_integration
#SBATCH --output=slurm_log/module5_integration_%j.out
#SBATCH --error=slurm_log/module5_integration_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Module 5: Cross-Sample Integration
# Integrates all 14 patient samples using Harmony-style batch correction

mkdir -p slurm_log

echo "=========================================="
echo "Module 5: Cross-Sample Integration"
echo "Date: $(date)"
echo "=========================================="

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python

OUTPUT_DIR="output/module5_integration"
mkdir -p ${OUTPUT_DIR}

# Run Module 5 integration
$PYTHON examples/run_module5_integration.py \
    --output-dir ${OUTPUT_DIR} \
    --n-components 30 \
    --n-clusters 50 \
    --theta 2.0 \
    --similarity-threshold 0.7 \
    --seed 42

echo "=========================================="
echo "Module 5 Complete"
echo "Date: $(date)"
echo "=========================================="
