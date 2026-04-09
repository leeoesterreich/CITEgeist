#!/bin/bash
#SBATCH --job-name=patient_phase4
#SBATCH --array=0-1
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=output/patient_pipeline/logs/phase4_%A_%a.out
#SBATCH --error=output/patient_pipeline/logs/phase4_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
export HDF5_USE_FILE_LOCKING=FALSE
export TF_FORCE_GPU_ALLOW_GROWTH=true

# Task 0 = baseline, Task 1 = cellularity
if [ "${SLURM_ARRAY_TASK_ID}" -eq 0 ]; then
    VARIANT="baseline"
else
    VARIANT="cellularity"
fi

echo "Running Phase 4 (MIL training + Hungarian assignment) variant=${VARIANT}"

mkdir -p output/patient_pipeline/logs

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist

python examples/run_patient_phase4.py \
    --output-dir output/patient_pipeline/phase4 \
    --mod3-dir "output/patient_pipeline/phase2/${VARIANT}" \
    --features-dir output/patient_pipeline/phase3 \
    --seg-dir output/patient_pipeline/phase1 \
    --variant "${VARIANT}"

echo "Phase 4 complete for variant=${VARIANT}"
