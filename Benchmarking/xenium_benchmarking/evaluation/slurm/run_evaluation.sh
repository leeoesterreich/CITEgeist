#!/bin/bash
#SBATCH --job-name=eval_discovery
#SBATCH --output=slurm_log/%x.out
#SBATCH --error=slurm_log/%x.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
cd "${REPO_ROOT}"

python Benchmarking/xenium_benchmarking/evaluation/src/evaluate_discovery_methods.py \
    --input-dir Benchmarking/xenium_benchmarking/evaluation/results/discovery_comparison

echo "Evaluation complete"
