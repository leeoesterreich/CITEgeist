#!/bin/bash
#SBATCH --job-name=unified_step2
#SBATCH --array=0-11
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=output/unified_pipeline/logs/step2_%A_%a.out
#SBATCH --error=output/unified_pipeline/logs/step2_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
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

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist

# Gate on Step 1 completion
MARKER="output/unified_pipeline/${SAMPLE}/.step1_complete"
if [ ! -f "${MARKER}" ]; then
    echo "Step 1 not complete for ${SAMPLE}, exiting"
    exit 1
fi

echo "Running Step 2 (StarDist + ViT features) for ${SAMPLE}"
python examples/run_unified_step2_features.py --sample ${SAMPLE} --modality he
