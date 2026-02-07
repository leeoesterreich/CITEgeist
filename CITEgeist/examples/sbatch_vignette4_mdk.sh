#!/bin/bash
#SBATCH --job-name=vignette4_mdk
#SBATCH --output=slurm_log/vignette4_mdk_%j.out
#SBATCH --error=slurm_log/vignette4_mdk_%j.err
#SBATCH --time=4:00:00
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Vignette 4: ESR1-D538G MDK Analysis
# Runs full spatial analysis + RNA-seq validation

echo "=========================================="
echo "Vignette 4: ESR1-D538G MDK Analysis"
echo "Start time: $(date)"
echo "Node: $(hostname)"
echo "=========================================="

# Environment setup
source ~/.bashrc
eval "$(conda shell.bash hook)"
module load gurobi/12.0.3
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Create log directory
mkdir -p slurm_log

# Change to examples directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/CITEgeist/examples

# Run the analysis
python run_vignette4_mdk_analysis.py

echo "=========================================="
echo "End time: $(date)"
echo "=========================================="
