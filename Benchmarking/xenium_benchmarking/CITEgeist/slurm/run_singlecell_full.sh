#!/bin/bash
#SBATCH --job-name=singlecell_full
#SBATCH --output=slurm_log/singlecell_full_%j.out
#SBATCH --error=slurm_log/singlecell_full_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell demonstration: Full dataset (465K cells)
# Runs Modules 1-2 then Module 4 sequentially

set -e

echo "=========================================="
echo "Single-Cell Demonstration: Full Dataset"
echo "Started: $(date)"
echo "=========================================="

# Load environment
module load python/ondemand-jupyter-python3.11

# Change to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration/Benchmarking/xenium_benchmarking/CITEgeist/src

# Run Module 1-2
echo ""
echo "Running Modules 1-2..."
python run_singlecell_module12.py --mode full

# Run Module 4
echo ""
echo "Running Module 4..."
python run_singlecell_module4.py --mode full

# Run evaluation
echo ""
echo "Running evaluation..."
python evaluate_singlecell.py --mode full

# Generate figures
echo ""
echo "Generating figures..."
python generate_singlecell_figures.py --mode full

echo ""
echo "=========================================="
echo "Completed: $(date)"
echo "=========================================="
