#!/bin/bash
#SBATCH --job-name=singlecell_Q0_ext
#SBATCH --output=slurm_log/singlecell_Q0_ext_%j.out
#SBATCH --error=slurm_log/singlecell_Q0_ext_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Resubmission of Q0 with extended time limit (12 hours)
# Q0 timed out at 8 hours during Module 2c Profile Selection

set -e

QUADRANT_ID=0

echo "=========================================="
echo "Single-Cell Demonstration: Quadrant ${QUADRANT_ID} (Extended)"
echo "Started: $(date)"
echo "Time limit: 12 hours"
echo "=========================================="

# Load environment
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Change to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration/Benchmarking/xenium_benchmarking/CITEgeist/src

# Run Module 1-2 (will restart from scratch)
echo ""
echo "Running Modules 1-2..."
python run_singlecell_module12.py --mode quadrant --quadrant-id ${QUADRANT_ID}

# Run Module 4
echo ""
echo "Running Module 4..."
python run_singlecell_module4.py --mode quadrant --quadrant-id ${QUADRANT_ID}

# Run evaluation
echo ""
echo "Running evaluation..."
python evaluate_singlecell.py --mode quadrant --quadrant-id ${QUADRANT_ID}

# Generate figures
echo ""
echo "Generating figures..."
python generate_singlecell_figures.py --mode quadrant --quadrant-id ${QUADRANT_ID}

echo ""
echo "=========================================="
echo "Quadrant ${QUADRANT_ID} Completed: $(date)"
echo "=========================================="
