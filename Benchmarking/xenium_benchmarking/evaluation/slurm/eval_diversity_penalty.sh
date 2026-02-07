#!/bin/bash
#SBATCH --job-name=eval_div_penalty
#SBATCH --output=slurm_log/eval_div_penalty_%j.out
#SBATCH --error=slurm_log/eval_div_penalty_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

python "${REPO_ROOT}/Benchmarking/xenium_benchmarking/evaluation/src/evaluate_gex.py" \
    --gt-dir "${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_granular_gt" \
    --pred-dir "${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_achievable_7_no_unknown" \
    --n-regions 1 \
    --prefix "Xenium" \
    --output "${REPO_ROOT}/Benchmarking/xenium_benchmarking/evaluation/results_gex_diversity_penalty.json"

echo "Done!"
