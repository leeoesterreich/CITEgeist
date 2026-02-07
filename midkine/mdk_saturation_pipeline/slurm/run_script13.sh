#!/bin/bash
#SBATCH --job-name=script13
#SBATCH --output=slurm_log/script13_%j.out
#SBATCH --error=slurm_log/script13_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# MDK Saturation Pipeline - Script 13: HSP90B1 Mechanism Discrimination

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine/mdk_saturation_pipeline

# Load conda environment
module load gcc/8.2.0
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Running Script 13: HSP90B1 Mechanism Discrimination"
echo "================================"
date

python scripts/13_discriminate_hsp90b1_mechanism.py

echo ""
echo "Script 13 finished"
date
