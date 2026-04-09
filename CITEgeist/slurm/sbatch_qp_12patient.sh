#!/bin/bash
#SBATCH --job-name=qp_12patient
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:1
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/slurm/slurm_log/qp_12patient_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/slurm/slurm_log/qp_12patient_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Re-run cuOPT QP for all 12 HCC22-088 samples with the updated Epithelial profile.
# Outputs: output/module3_cuopt_qp/{sample}/{sample}_cell_prop_global_results.csv
#
# Submit via gpu-race-submit (not direct sbatch) to race l40s/a100/a100_nvlink:
#   gpu-race-submit.sh CITEgeist/slurm/sbatch_qp_12patient.sh

set -eo pipefail

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
LOG_DIR="${REPO_ROOT}/CITEgeist/slurm/slurm_log"
DATA_DIR="/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files"
OUTPUT_DIR="${REPO_ROOT}/output/module3_cuopt_qp"
# GPU_RACE_NAMESPACE used for per-partition temp dir namespacing
NUMBA_CACHE_DIR="/tmp/numba_cache_${GPU_RACE_PARTITION:-default}_${SLURM_ARRAY_TASK_ID}"

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

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
echo "=== Array task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE} ==="
echo "Node: $(hostname), GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"

mkdir -p "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}/${SAMPLE}"

# Remove stale completion marker to force re-run with Epithelial profile
rm -f "${OUTPUT_DIR}/${SAMPLE}/.phase2_complete"

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

MPLBACKEND=Agg \
python CITEgeist/examples/run_patient_phase2.py \
    --sample "${SAMPLE}" \
    --output-dir "${OUTPUT_DIR}" \
    --data-dir "${DATA_DIR}"

echo "=== Complete: ${SAMPLE} ==="
