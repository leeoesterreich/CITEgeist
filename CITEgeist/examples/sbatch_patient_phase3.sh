#!/bin/bash
#SBATCH --job-name=patient_phase3
#SBATCH --array=0-11
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=output/patient_pipeline/logs/phase3_%A_%a.out
#SBATCH --error=output/patient_pipeline/logs/phase3_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
export HDF5_USE_FILE_LOCKING=FALSE
export TF_FORCE_GPU_ALLOW_GROWTH=true

SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
echo "Running Phase 3 (ViT feature extraction) for ${SAMPLE}"

# Gate on Phase 1 completion
PHASE1_MARKER="output/patient_pipeline/phase1/${SAMPLE}/.phase1_complete"
if [ ! -f "${PHASE1_MARKER}" ]; then
    echo "Phase 1 not complete for ${SAMPLE}, exiting"
    exit 1
fi

mkdir -p output/patient_pipeline/logs

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist

python examples/run_patient_phase3.py \
    --sample "${SAMPLE}" \
    --output-dir output/patient_pipeline/phase3 \
    --seg-dir output/patient_pipeline/phase1 \
    --data-dir /ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files

echo "Phase 3 complete for ${SAMPLE}"
