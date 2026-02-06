#!/bin/bash
#SBATCH --job-name=singlecell_quad
#SBATCH --output=slurm_log/singlecell_quad_%A_%a.out
#SBATCH --error=slurm_log/singlecell_quad_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --array=0-3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell demonstration: Quadrants (array job)
# Each array task processes one quadrant (~116K cells)

set -e

QUADRANT_ID=$SLURM_ARRAY_TASK_ID

echo "=========================================="
echo "Single-Cell Demonstration: Quadrant ${QUADRANT_ID}"
echo "Started: $(date)"
echo "=========================================="

# Load environment
module load python/ondemand-jupyter-python3.11

# Change to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/singlecell-demonstration/Benchmarking/xenium_benchmarking/CITEgeist/src

# Run Module 1-2
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
