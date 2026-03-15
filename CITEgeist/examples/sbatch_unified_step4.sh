#!/bin/bash
#SBATCH --job-name=unified_step4
#SBATCH --array=0-11
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --output=output/unified_pipeline/logs/step4_%A_%a.out
#SBATCH --error=output/unified_pipeline/logs/step4_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

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

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

MARKER="output/unified_pipeline/${SAMPLE}/.step3_complete"
if [ ! -f "${MARKER}" ]; then
    echo "Step 3 not complete for ${SAMPLE}, exiting"
    exit 1
fi

echo "Running Step 4 (Marker gene validation) for ${SAMPLE}"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist
python examples/run_unified_step4_validate.py --sample ${SAMPLE}
