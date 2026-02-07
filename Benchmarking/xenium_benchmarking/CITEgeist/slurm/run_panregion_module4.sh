#!/bin/bash
#SBATCH --job-name=panregion_m4
#SBATCH --output=slurm_log/panregion_m4_%j.out
#SBATCH --error=slurm_log/panregion_m4_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Pan-Region Module 4: Consensus Profile Program Discovery
# Applies unified profiles across all 465K cells

set -e

echo "=========================================="
echo "Pan-Region Module 4: Consensus Profiles"
echo "Started: $(date)"
echo "=========================================="

# Load environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Change to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration/Benchmarking/xenium_benchmarking/CITEgeist/src

echo ""
echo "Running Pan-Region Module 4..."
python run_panregion_module4.py \
    --n-programs 5 \
    --min-cells 1000

echo ""
echo "=========================================="
echo "Pan-Region Module 4 Completed: $(date)"
echo "=========================================="
