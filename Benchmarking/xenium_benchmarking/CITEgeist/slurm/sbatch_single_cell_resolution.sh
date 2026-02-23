#!/bin/bash
#SBATCH --job-name=sc_resolution
#SBATCH --output=logs/sc_resolution_%A_%a.out
#SBATCH --error=logs/sc_resolution_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell resolution benchmark using Module 3b nucleus assignment pipeline.
#
# This script tests the new single-cell resolution features:
# 1. Loads existing hybrid results (proportions + deconvolved GEX)
# 2. Runs Cellpose to get full nuclei masks (not just counts)
# 3. Maps nuclei to spots
# 4. Runs run_nucleus_assignment() pipeline (morphology features + soft-label classifier + Hungarian assignment)
# 5. Distributes GEX to individual cells
# 6. Creates single-cell AnnData output
#
# Requires: Existing hybrid_cellpose results for each region

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Check that hybrid results exist
HYBRID_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_cellpose"
REGION_DIR="${HYBRID_DIR}/Xenium_region_${SLURM_ARRAY_TASK_ID}"
if [ ! -d "$REGION_DIR" ]; then
    echo "ERROR: Hybrid results not found at $REGION_DIR"
    echo "Please run sbatch_hybrid_cellpose_gpu.sh first"
    exit 1
fi

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution \
    --use-gpu \
    --spot-diameter-um 55.0
