#!/bin/bash
#SBATCH --job-name=analyze_clusters
#SBATCH --output=slurm_log/analyze_clusters_%j.out
#SBATCH --error=slurm_log/analyze_clusters_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

# Analyze RNA cluster protein profiles for cell type annotation

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium

# Set up environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

echo "Starting cluster profile analysis..."
echo "Job ID: $SLURM_JOB_ID"
echo "Time: $(date)"
echo "Python: $PYTHON"

$PYTHON analyze_cluster_profiles.py

echo "Completed at: $(date)"
