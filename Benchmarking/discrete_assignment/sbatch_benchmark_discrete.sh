#!/bin/bash
#SBATCH --job-name=discrete_bench
#SBATCH --output=discrete_bench_%j.out
#SBATCH --error=discrete_bench_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load Gurobi module (REQUIRED for discrete cell assignment)
module load gurobi/12.0.3

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

# Change to project directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/.worktrees/feat-discrete-cell-assignment

# Run benchmark on Xenium with GEX evaluation
python Benchmarking/discrete_assignment/benchmark_discrete.py \
    --dataset xenium \
    --region 0 \
    --output-dir output/discrete_benchmark \
    --run-gex

echo "Benchmark complete!"
