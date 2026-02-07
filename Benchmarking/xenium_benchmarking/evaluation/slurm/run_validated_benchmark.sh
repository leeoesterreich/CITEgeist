#!/bin/bash
#SBATCH --job-name=validated_benchmark
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/validated_benchmark_%A_%a.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/validated_benchmark_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=smp
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo "Starting validated benchmark for region ${SLURM_ARRAY_TASK_ID} at $(date)"

# Activate conda environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Load Gurobi module with valid license
module load gurobi/12.0.3

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run benchmark for this region
python Benchmarking/xenium_benchmarking/evaluation/src/run_validated_benchmark.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output_dir /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/validated_benchmark \
    --verbose

echo "Finished region ${SLURM_ARRAY_TASK_ID} at $(date)"
