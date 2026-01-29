#!/bin/bash
#SBATCH --job-name=test_beta_fixes
#SBATCH --output=slurm_log/test_beta_fixes_%j.out
#SBATCH --error=slurm_log/test_beta_fixes_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

echo "Running per-marker beta bug fix tests..."
echo "Python: $(which python)"
echo "Gurobi license: ${GRB_LICENSE_FILE:-default}"
echo ""

python -m pytest tests/test_per_marker_beta_fixes.py -v --tb=short 2>&1

echo ""
echo "Exit code: $?"
