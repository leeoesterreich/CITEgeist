#!/bin/bash
#SBATCH --job-name=spatial_panels
#SBATCH --output=slurm_log/spatial_panels_%j.out
#SBATCH --error=slurm_log/spatial_panels_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Generate spatial visualization panels for manuscript figures

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures

# Activate environment
source /ihome/alee/alc376/.bashrc
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Create log directory
mkdir -p slurm_log

# Run the spatial panel generation
python generate_spatial_panels.py

echo "Done!"
