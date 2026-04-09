#!/bin/bash
#SBATCH --job-name=patient_phase5_validate
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --output=output/patient_pipeline/logs/phase5_%j.out
#SBATCH --error=output/patient_pipeline/logs/phase5_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Running Phase 5 (validation)"

mkdir -p output/patient_pipeline/logs

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist-sc-assignment/CITEgeist

python examples/run_patient_phase5_validate.py \
    --phase4-dir output/patient_pipeline/phase4 \
    --output-dir output/patient_pipeline/phase5_validation \
    --data-dir /ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files

echo "Phase 5 validation complete"
