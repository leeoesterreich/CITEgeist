#!/bin/bash
#SBATCH --job-name=module5_pseudobulk
#SBATCH --output=slurm_log/module5_pseudobulk_%j.out
#SBATCH --error=slurm_log/module5_pseudobulk_%j.err
#SBATCH --time=4:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Module 5 PyDESeq2 Analysis - Pseudo-bulk Approach
# ==================================================
# Uses sample-level aggregation for faster, more appropriate DE analysis

echo "============================================"
echo "Module 5 PyDESeq2 (Pseudo-bulk)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Started: $(date)"
echo "============================================"

# Load conda
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Set working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Ensure output directories exist
mkdir -p examples/slurm_log
mkdir -p examples/output_module5_pydeseq

# Run analysis
echo ""
echo "Starting pseudo-bulk PyDESeq2 analysis..."
python examples/run_module5_pydeseq_pseudobulk.py

echo ""
echo "============================================"
echo "Completed: $(date)"
echo "============================================"
