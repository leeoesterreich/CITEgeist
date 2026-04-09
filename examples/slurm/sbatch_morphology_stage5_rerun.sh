#!/bin/bash
#SBATCH --job-name=morph_s5
#SBATCH --output=output/morphology_assignment_v3/logs/stage5_%A_%a.out
#SBATCH --error=output/morphology_assignment_v3/logs/stage5_%A_%a.err
#SBATCH --array=0-11
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Stage 5 only re-run — fixes SACE GEX antibody/GEX spot mismatch

export HDF5_USE_FILE_LOCKING=FALSE

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
echo "Stage 5 re-run: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID})"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python examples/run_morphology_assignment.py \
    --sample "${SAMPLE}" \
    --stages 5 \
    --module3-dir output/module3_cuopt_qp \
    --output-dir-override output/morphology_assignment_v3/${SAMPLE} \
    --embeddings-dir output/morphology_assignment_v3/${SAMPLE}

echo "Done: ${SAMPLE}"
