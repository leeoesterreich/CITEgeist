#!/bin/bash
#SBATCH --job-name=singlecell_validation
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/singlecell_validation_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/slurm/slurm_log/singlecell_validation_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --partition=smp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

echo "Starting single-cell validation analysis at $(date)"

# Activate conda environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run analysis
python Benchmarking/xenium_benchmarking/evaluation/src/singlecell_validation.py \
    --output_dir /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/evaluation/singlecell_validation \
    --verbose

echo "Finished at $(date)"
