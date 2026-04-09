#!/bin/bash
#SBATCH --job-name=sim_div_v4
#SBATCH --output=slurm_log/sim_div_v4_%j.out
#SBATCH --error=slurm_log/sim_div_v4_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Simulated benchmark v4: full fix of lambda_reg_gex wiring
# Fixed in: gurobi_impl.py, citegeist_model.py
# Sequential runs with checkpoint clearing between them
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

# Clear any existing shared checkpoints for Wu_rep
rm -f checkpoints/Wu_rep_*_gene_expression_*.npz 2>/dev/null

echo "Running simulated benchmark WITH diversity penalty (lambda_reg_gex=0.01)..."
python tests/test_citegeist_simulated.py \
    --radius 4.0 \
    --lambda_reg 1.0 \
    --alpha_elastic 0.7 \
    --max_y_change 0.05 \
    --input_folder replicates/high_seg \
    --output_folder test_results/diversity_penalty_v4 \
    --sample_prefix Wu_rep \
    --lambda-laplacian 0.0 \
    --lambda-reg-gex 0.01

echo "=== Diversity penalty done ==="

# Clear checkpoints for clean baseline run
rm -f checkpoints/Wu_rep_*_gene_expression_*.npz 2>/dev/null

echo "Running simulated benchmark WITHOUT diversity penalty (baseline)..."
python tests/test_citegeist_simulated.py \
    --radius 4.0 \
    --lambda_reg 1.0 \
    --alpha_elastic 0.7 \
    --max_y_change 0.05 \
    --input_folder replicates/high_seg \
    --output_folder test_results/baseline_v4 \
    --sample_prefix Wu_rep \
    --lambda-laplacian 0.0 \
    --lambda-reg-gex 0.0

echo "=== Baseline done ==="
echo "Both runs complete!"
