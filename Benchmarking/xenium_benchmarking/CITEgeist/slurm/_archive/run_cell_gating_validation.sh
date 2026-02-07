#!/bin/bash
#SBATCH --job-name=cell_gating_validation
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_log/cell_gating_%A_%a.out
#SBATCH --error=slurm_log/cell_gating_%A_%a.err
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Run gating-based cell classification validation on Xenium RCC regions 0-4
# Uses new cell_classification.py (flow-cytometry Boolean gating)

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

if [ $? -ne 0 ]; then
    echo "Error: Failed to activate CITEgeist environment!" >&2
    exit 1
fi

module load gurobi/12.0.3

echo "========================================"
echo "Cell Gating Validation - Region ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID}, Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname), Start: $(date)"
echo "========================================"

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd $REPO_ROOT

python Benchmarking/xenium_benchmarking/CITEgeist/src/run_cell_gating_validation.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --max-cells 10000

echo "========================================"
echo "Completed: $(date)"
echo "========================================"
