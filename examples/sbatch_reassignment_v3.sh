#!/bin/bash
#SBATCH --job-name=reassign_v3
#SBATCH --output=slurm_log/reassign_v3_%A_%a.out
#SBATCH --error=slurm_log/reassign_v3_%A_%a.err
#SBATCH --array=0-11
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Re-run Stage 5+6 (bayesian assignment + SACE GEX) for all 12 patients.
# Stages 1-2 (segmentation + embedding) reuse existing output from v3.
# Writes to output/morphology_assignment_v3/{sample}/ (overwrites in-place).

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
REPO=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
PYTHON=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python
OUTPUT_DIR=${REPO}/output/morphology_assignment_v3/${SAMPLE}

echo "=========================================="
echo "Bayesian Re-assignment (Stage 5+6)"
echo "Sample: ${SAMPLE}"
echo "Output: ${OUTPUT_DIR}"
echo "Date: $(date)"
echo "=========================================="

mkdir -p ${REPO}/slurm_log
cd ${REPO}

# Check segmentation exists before submitting stage 5
if [ ! -f "${OUTPUT_DIR}/nucleus_spot_mapping.csv" ]; then
    echo "ERROR: nucleus_spot_mapping.csv missing — run stages 1,2 first"
    exit 1
fi

${PYTHON} examples/run_morphology_assignment.py \
    --sample ${SAMPLE} \
    --stages 5,6 \
    --output-dir-override ${OUTPUT_DIR} \
    --device cpu

echo "=========================================="
echo "Done: ${SAMPLE} at $(date)"
echo "=========================================="
