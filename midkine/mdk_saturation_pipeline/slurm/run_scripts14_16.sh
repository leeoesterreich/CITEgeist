#!/bin/bash
#SBATCH --job-name=fig14_16
#SBATCH --output=slurm_log/fig14_16_%j.out
#SBATCH --error=slurm_log/fig14_16_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --cluster=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/midkine/mdk_saturation_pipeline

module load gcc/8.2.0
source /ix1/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/etc/profile.d/conda.sh
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Running Scripts 14-16: Attractor State Evidence Figures"
echo "======================================================="
date

python scripts/14_cistrome_divergence.py && \
python scripts/15_proteostasis_setpoint.py && \
python scripts/16_functional_consequence.py

echo ""
echo "All figure scripts finished"
date
