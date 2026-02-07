#!/bin/bash
#SBATCH --job-name=gex_lam_sweep
#SBATCH --output=slurm_log/gex_lambda_sweep_%j.out
#SBATCH --error=slurm_log/gex_lambda_sweep_%j.err
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e
module load gurobi/12.0.3

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

export PYTHONPATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist:$PYTHONPATH"

echo "=============================================="
echo "GEX Lambda Sweep - Region 0"
echo "=============================================="
echo "Start time: $(date)"

python Benchmarking/xenium_benchmarking/CITEgeist/src/sweep_gex_lambda.py \
    --region-id 0 \
    --lambdas "0.01,0.1,1.0,5.0,10.0,50.0,100.0"

echo ""
echo "End time: $(date)"
