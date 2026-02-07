#!/bin/bash
#SBATCH --job-name=script10
#SBATCH --output=slurm_log/script10_%j.out
#SBATCH --error=slurm_log/script10_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# MDK Saturation Pipeline
# Runs scripts 02-07 (skipping 01 which needs vignette4 outputs)

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine/mdk_saturation_pipeline

# Load conda environment
module load gcc/8.2.0
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Running Script 10: Unsupervised HSP90B1 Identification"
echo "================================"
date

# Run pipeline starting from script 02
python scripts/10_unsupervised_hsp90b1_identification.py

echo ""
echo "Pipeline finished"
date
