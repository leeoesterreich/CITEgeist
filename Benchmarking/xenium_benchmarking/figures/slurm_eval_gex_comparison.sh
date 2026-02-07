#!/bin/bash
#SBATCH --job-name=eval_gex_comp
#SBATCH --output=slurm_log/eval_gex_comp_%j.out
#SBATCH --error=slurm_log/eval_gex_comp_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

python "${REPO_ROOT}/Benchmarking/xenium_benchmarking/figures/plot_gex_comparison.py"

echo "Done!"
