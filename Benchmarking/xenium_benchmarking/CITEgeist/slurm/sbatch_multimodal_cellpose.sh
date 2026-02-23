#!/bin/bash
#SBATCH --job-name=multimodal_xenium
#SBATCH --output=logs/multimodal_xenium_%A_%a.out
#SBATCH --error=logs/multimodal_xenium_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Multimodal approach:
# 1. Run continuous model (CLR preprocessing, QP optimization) - Pass 1
# 2. Run multimodal EM refinement (anchor gene learning + iterative refinement) - Pass 1.5 + 2
# 3. Discretize using nuclei counts from Cellpose
# 4. Run GEX deconvolution
#
# This tests whether RNA-based EM refinement improves protein-based deconvolution
# for cells with RNA signal but low protein expression

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_multimodal_cellpose.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/multimodal_cellpose \
    --lambda-laplacian 0.1 \
    --spot-diameter-um 55.0 \
    --n-anchors 20 \
    --min-correlation 0.3 \
    --lambda-prior 1.0 \
    --max-iterations 20
