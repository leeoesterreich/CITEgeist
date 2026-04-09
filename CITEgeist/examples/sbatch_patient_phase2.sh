#!/bin/bash
#SBATCH --job-name=patient_phase2
#SBATCH --array=0-23
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --output=output/patient_pipeline/logs/phase2_%A_%a.out
#SBATCH --error=output/patient_pipeline/logs/phase2_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

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

# Tasks 0-11 = baseline, 12-23 = cellularity
N_SAMPLES=12
SAMPLE_IDX=$(( SLURM_ARRAY_TASK_ID % N_SAMPLES ))
VARIANT_IDX=$(( SLURM_ARRAY_TASK_ID / N_SAMPLES ))

SAMPLE=${SAMPLES[$SAMPLE_IDX]}

if [ "${VARIANT_IDX}" -eq 0 ]; then
    VARIANT="baseline"
    CELLULARITY_FLAG=""
else
    VARIANT="cellularity"
    CELLULARITY_FLAG="--use-cellularity-prior"
fi

echo "Running Phase 2 (Module 3 deconvolution) for ${SAMPLE} [variant=${VARIANT}]"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist

# Poll for Phase 1 completion (cross-cluster, no SLURM dependency)
PHASE1_MARKER="output/patient_pipeline/phase1/${SAMPLE}/.phase1_complete"
echo "Waiting for Phase 1 marker: ${PHASE1_MARKER}"
for i in $(seq 1 60); do
    if [ -f "${PHASE1_MARKER}" ]; then
        echo "Phase 1 complete for ${SAMPLE}"
        break
    fi
    if [ $i -eq 60 ]; then
        echo "ERROR: Phase 1 not complete for ${SAMPLE} after 60 min, exiting"
        exit 1
    fi
    sleep 60
done

mkdir -p output/patient_pipeline/logs

python examples/run_patient_phase2.py \
    --sample "${SAMPLE}" \
    --output-dir "output/patient_pipeline/phase2/${VARIANT}" \
    --seg-dir "output/patient_pipeline/phase1" \
    --data-dir /ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files \
    ${CELLULARITY_FLAG}

echo "Phase 2 complete for ${SAMPLE} [${VARIANT}]"
