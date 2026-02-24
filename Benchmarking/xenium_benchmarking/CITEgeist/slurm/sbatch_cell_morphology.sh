#!/bin/bash
#SBATCH --job-name=cell_morph
#SBATCH --output=logs/cell_morph_%A_%a.out
#SBATCH --error=logs/cell_morph_%A_%a.err
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Cell Morphology Benchmark
# Runs Cellpose + watershed segmentation for 5 Xenium regions
# Extracts nuclear and cell-level morphology features

set -e

echo "=== Cell Morphology Benchmark ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo ""

# Activate conda environment
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Set working directory
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
SCRIPT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output/cell_morphology"

cd "${SCRIPT_DIR}"

echo "Working directory: $(pwd)"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Run benchmark
python benchmark_cell_morphology.py \
    --region "${SLURM_ARRAY_TASK_ID}" \
    --output-dir "${OUTPUT_DIR}"

echo ""
echo "=== Job Complete ==="
echo "Exit code: $?"
echo "Date: $(date)"
