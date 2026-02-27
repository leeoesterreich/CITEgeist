#!/bin/bash
#SBATCH --job-name=det_est
#SBATCH --output=logs/det_est_%A_%a.out
#SBATCH --error=logs/det_est_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Detection + Estimation Benchmark
# Runs two-stage model: GMM detection -> masked IQP estimation
#
# Usage:
#   cd Benchmarking/xenium_benchmarking/CITEgeist/slurm
#   mkdir -p logs
#   sbatch sbatch_detection_estimation.sh

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate /ihome/alee/alc376/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_detection_estimation.py \
    --region-id ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/detection_estimation \
    --max-iter 10 \
    --detection-threshold 0.5

echo "Region ${SLURM_ARRAY_TASK_ID} completed"
