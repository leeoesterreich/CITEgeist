#!/bin/bash
#SBATCH --job-name=stage2_morph
#SBATCH --output=logs/stage2_morph_%A_%a.out
#SBATCH --error=logs/stage2_morph_%A_%a.err
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --array=0-4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Two-Stage Morphology Benchmark
# Stage 1: Hybrid proportions (pre-computed)
# Stage 2: Morphology-based single-cell assignment

set -e

# Setup environment
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Paths
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
SCRIPT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output/two_stage_morphology"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "$(dirname $0)/logs"

# Run benchmark for this region
REGION=${SLURM_ARRAY_TASK_ID}

echo "=========================================="
echo "Two-Stage Morphology Benchmark"
echo "Region: ${REGION}"
echo "=========================================="

python "${SCRIPT}" \
    --region "${REGION}" \
    --output-dir "${OUTPUT_DIR}"

echo "Done!"
