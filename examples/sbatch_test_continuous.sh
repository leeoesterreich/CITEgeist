#!/bin/bash
#SBATCH --job-name=test_continuous
#SBATCH --output=slurm_log/test_continuous_%j.out
#SBATCH --error=slurm_log/test_continuous_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/examples
mkdir -p ../output/module12_discovery_continuous

# Limit joblib workers to match allocated CPUs
export CITEGEIST_N_JOBS=16

/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin/python run_module12_discovery.py \
    --sample HCC22-088-P1-S1 \
    --output-dir ../output/module12_discovery_continuous \
    --top-k 3

echo "Done"
