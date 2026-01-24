#!/bin/bash
#SBATCH --job-name=citegeist_hierarchical
#SBATCH --output=slurm_log/hierarchical_%A_%a.out
#SBATCH --error=slurm_log/hierarchical_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Hierarchical Profile Discovery Benchmark on Xenium Data
# Tests the new hierarchical NMF-based profile discovery

echo "======================================"
echo "HIERARCHICAL PROFILE DISCOVERY BENCHMARK"
echo "======================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start time: $(date)"
echo ""

# Activate environment (using module system)
module purge
module load gcc/12.2.0
module load gurobi/12.0.3
module load anaconda3/2023.09-0-python_3.11.5
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Set paths
CITEGEIST_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/hierarchical-profiles/Benchmarking/xenium_benchmarking/CITEgeist"
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_granular_gt"
OUTPUT_DIR="${CITEGEIST_DIR}/output_hierarchical"
PYTHON="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${CITEGEIST_DIR}/slurm_log"

# Region ID from array task
REGION_ID=${SLURM_ARRAY_TASK_ID}

echo "Processing region ${REGION_ID}..."
echo "Input: ${INPUT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Run hierarchical benchmark
${PYTHON} "${CITEGEIST_DIR}/src/run_benchmark_hierarchical.py" \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --radius 4.0 \
    --lambda-reg 1.0 \
    --alpha-elastic 0.7 \
    --max-y-change 0.4 \
    --min-counts 25 \
    --improvement-threshold 0.05 \
    --max-depth 5

exit_code=$?

echo ""
echo "======================================"
echo "Benchmark complete"
echo "Exit code: ${exit_code}"
echo "End time: $(date)"
echo "======================================"

exit ${exit_code}
