#!/bin/bash
#SBATCH --job-name=staged_eval
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/staged_eval_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/staged_eval_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Module-by-Module Pipeline Evaluation for CITEgeist
# Runs staged evaluation on all 5 Xenium regions in parallel
#
# Each stage produces diagnostic outputs:
#   Stage 1: Module 1 (Marker Interest Detection)
#   Stage 2a: Module 2a (Colocalization)
#   Stage 2b: Module 2b (Profile Discovery)
#   Stage 2c: Module 2c (Profile Selection)
#   Stage 3b: Profile-to-GT Gap Analysis
#   Stage 4: Oracle vs Auto-Discovery Comparison

echo "=============================================="
echo "STAGED PIPELINE EVALUATION"
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Running on: $(hostname)"
echo "Date: $(date)"
echo "=============================================="

# Activate conda environment using absolute path
source ~/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Load Gurobi module for license
module load gurobi/12.0.3

# Set paths using absolute paths
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
SRC_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src"
INPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_protein_gt"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_staged_evaluation"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Get region ID from array task
REGION_ID=$SLURM_ARRAY_TASK_ID

echo "Processing Region: ${REGION_ID}"
echo "Input Directory: ${INPUT_DIR}"
echo "Output Directory: ${OUTPUT_DIR}"
echo "SRC Directory: ${SRC_DIR}"
echo ""

# Run staged evaluation
cd "${SRC_DIR}"
python evaluate_pipeline_stages.py \
    --region-id ${REGION_ID} \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}"

EXIT_CODE=$?

echo ""
echo "=============================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Staged evaluation completed successfully for Region ${REGION_ID}"
else
    echo "Staged evaluation FAILED for Region ${REGION_ID} with exit code ${EXIT_CODE}"
fi
echo "=============================================="
echo "End time: $(date)"
