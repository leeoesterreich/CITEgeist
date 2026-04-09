#!/bin/bash
#SBATCH --job-name=morph_sd
#SBATCH --output=output/morphology_assignment_v3/logs/morph_%A_%a.out
#SBATCH --error=output/morphology_assignment_v3/logs/morph_%A_%a.err
#SBATCH --array=0-11
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Full pipeline re-run with StarDist 3-layer defense segmentation.
# Replaces Cellpose (March 11) with production StarDist:
#   Layer 1 (Physics): auto scale from scalefactors (pixel_size = 55µm / spot_diam_fullres)
#   Layer 2 (Model): prob_thresh=0.6 (conservative, removes debris)
#   Layer 3 (Biology): area filter 20-500 µm² (removes sub-nuclear + merged objects)
#
# Runs ALL 6 stages from scratch to output/morphology_assignment_v3/

# Disable HDF5 file locking — NFS doesn't support it and causes OSError
# when multiple array tasks load StarDist model weights simultaneously
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
echo "Processing sample: ${SAMPLE} (task ${SLURM_ARRAY_TASK_ID})"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Create directories
mkdir -p output/morphology_assignment_v3/logs
mkdir -p output/morphology_assignment_v3/${SAMPLE}

# Full pipeline: stages 1-6
# Stage 1: patchwise StarDist on fullres (0.79 µm/px, ~132px patches per spot)
#   Crops 55µm+25µm context around each spot, runs StarDist at native resolution
#   3-layer defense: auto scale=3.2x, prob_thresh=0.6, area filter 20-500 µm²
# Stage 2: centroid-based patch extraction (no global mask needed)
python examples/run_morphology_assignment.py \
    --sample "${SAMPLE}" \
    --stages 1,2,3,4,5,6 \
    --module3-dir output/module3_cuopt_qp \
    --output-dir-override output/morphology_assignment_v3/${SAMPLE} \
    --pooled-checkpoint output/morphology_assignment_v2/pooled_mil_checkpoint.pt \
    --epochs 100 \
    --lr 1e-3

echo "Done: ${SAMPLE}"
