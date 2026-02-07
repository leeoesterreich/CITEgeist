#!/bin/bash
#SBATCH --job-name=sim_baseline
#SBATCH --output=slurm_log/sim_baseline_%j.out
#SBATCH --error=slurm_log/sim_baseline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Simulated benchmark baseline: NO diversity penalty (lambda_reg_gex=0.0)
# Same settings but with original winner-take-all GEX objective
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

echo "Running simulated benchmark WITHOUT diversity penalty (baseline)..."
python tests/test_citegeist_simulated.py \
    --radius 4.0 \
    --lambda_reg 1.0 \
    --alpha_elastic 0.7 \
    --max_y_change 0.05 \
    --input_folder replicates/high_seg \
    --output_folder test_results/baseline_no_diversity \
    --sample_prefix Wu_rep \
    --lambda-laplacian 0.0 \
    --lambda-reg-gex 0.0

echo "Done!"
