#!/bin/bash
#SBATCH --job-name=integration_compare
#SBATCH --output=slurm_log/integration_compare_%j.out
#SBATCH --error=slurm_log/integration_compare_%j.err
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Compare protein-anchored vs RNA-only program discovery

echo "=============================================="
echo "Integration Comparison: CITEgeist vs RNA-only"
echo "=============================================="
echo "Date: $(date)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Hostname: $(hostname)"
echo "=============================================="

# Activate conda environment
module purge
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Navigate to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/src

# Create log directory if needed
mkdir -p ../slurm/slurm_log

# Run comparison for all regions
python compare_protein_rna_integration.py --seed 42

echo "=============================================="
echo "Comparison complete!"
echo "Date: $(date)"
echo "=============================================="
