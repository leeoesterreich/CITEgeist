#!/bin/bash
#SBATCH --job-name=gen_cellpose_img
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --array=0-4
#SBATCH --output=logs/gen_cellpose_img_%A_%a.out
#SBATCH --error=logs/gen_cellpose_img_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Generate Cellpose-compatible images for simulation benchmarking
# Usage: sbatch sbatch_generate_images.sh high_seg
#        sbatch sbatch_generate_images.sh mixed

CONDITION=${1:-high_seg}
REPLICATE_ID=$SLURM_ARRAY_TASK_ID

eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR=Benchmarking/simulation_benchmarking/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/generate_cellpose_images.py \
    --replicate-id $REPLICATE_ID \
    --condition $CONDITION \
    --output-dir $OUTPUT_DIR \
    --mode both \
    --image-size 8000

echo "Completed replicate $REPLICATE_ID for condition $CONDITION"
