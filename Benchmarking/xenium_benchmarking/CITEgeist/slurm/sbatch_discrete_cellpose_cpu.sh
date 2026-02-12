#!/bin/bash
#SBATCH --job-name=cellpose_discrete_cpu
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/cellpose_discrete_cpu_%A_%a.out
#SBATCH --error=logs/cellpose_discrete_cpu_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

REGION_ID=$SLURM_ARRAY_TASK_ID
OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py \
    --region $REGION_ID \
    --output-dir $OUTPUT_DIR \
    --max-em-iterations 20

echo "Completed region $REGION_ID"
