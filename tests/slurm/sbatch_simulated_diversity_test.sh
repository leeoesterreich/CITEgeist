#!/bin/bash
#SBATCH --job-name=sim_div_test
#SBATCH --output=slurm_log/sim_div_test_%j.out
#SBATCH --error=slurm_log/sim_div_test_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Simulated benchmark: test quadratic diversity penalty (lambda_reg_gex=0.01)
# Runs 1 replicate (Wu_rep_0) with current Xenium-validated settings
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

echo "Running simulated benchmark with diversity penalty (lambda_reg_gex=0.01)..."
python tests/test_citegeist_simulated.py \
    --radius 4.0 \
    --lambda_reg 1.0 \
    --alpha_elastic 0.7 \
    --max_y_change 0.05 \
    --input_folder replicates/high_seg \
    --output_folder test_results/diversity_penalty_test \
    --sample_prefix Wu_rep \
    --lambda-laplacian 0.0 \
    --lambda-reg-gex 0.01

echo "Done!"
