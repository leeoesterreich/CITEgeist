#!/bin/bash
#SBATCH --job-name=plot_benchmark
#SBATCH --output=slurm_log/plot_benchmark_%j.out
#SBATCH --error=slurm_log/plot_benchmark_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Activate conda environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/figures

python plot_benchmark_overview.py

echo "Done!"
