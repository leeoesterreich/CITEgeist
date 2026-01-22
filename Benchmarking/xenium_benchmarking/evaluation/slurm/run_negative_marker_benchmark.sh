#!/bin/bash
#SBATCH --job-name=negmarker_bench
#SBATCH --output=slurm_log/negmarker_bench_%A_%a.out
#SBATCH --error=slurm_log/negmarker_bench_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Load Gurobi
module load gurobi/12.0.3

# Activate conda environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Navigate to script directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/src

# Output directory
OUTPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/negmarker_benchmark"

# Run benchmark for this region
REGION=$SLURM_ARRAY_TASK_ID

echo "Running negative marker benchmark for region $REGION"
echo "Output: $OUTPUT_DIR"

python run_negative_marker_benchmark.py \
    --region $REGION \
    --output_dir $OUTPUT_DIR \
    --lambda_neg 1.0 \
    --residual_threshold 0.10 \
    --fdr_alpha 0.05

echo "Benchmark complete for region $REGION"
