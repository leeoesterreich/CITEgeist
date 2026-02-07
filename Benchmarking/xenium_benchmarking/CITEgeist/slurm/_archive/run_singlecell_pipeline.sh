#!/bin/bash
#SBATCH --job-name=citegeist_singlecell
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/singlecell_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/singlecell_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# =============================================================================
# CITEgeist Single-Cell Pipeline on Xenium Data
# =============================================================================
#
# This script runs the full CITEgeist pipeline (Modules 1-4) on Xenium
# single-cell data to demonstrate resolution independence.
#
# Usage:
#   sbatch run_singlecell_pipeline.sh
#   sbatch --array=0 run_singlecell_pipeline.sh  # Single region for testing
#
# =============================================================================

set -e

# Load environment
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Navigate to project root
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Create log directory
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log

# Get region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

echo "=============================================="
echo "CITEgeist Single-Cell Pipeline"
echo "=============================================="
echo "Date: $(date)"
echo "Region: ${REGION_ID}"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Hostname: $(hostname)"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "=============================================="

# Run full pipeline
python Benchmarking/xenium_benchmarking/CITEgeist/src/run_singlecell_pipeline.py \
    --region ${REGION_ID} \
    --k-programs 5 \
    --seed 42

echo ""
echo "=============================================="
echo "Pipeline complete!"
echo "Date: $(date)"
echo "=============================================="
