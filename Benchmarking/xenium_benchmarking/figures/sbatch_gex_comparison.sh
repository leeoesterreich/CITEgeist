#!/bin/bash
#SBATCH --job-name=gex_comparison
#SBATCH --output=slurm_log/gex_comparison_%j.out
#SBATCH --error=slurm_log/gex_comparison_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/figures

python plot_gex_comparison.py

echo "Done!"
