#!/bin/bash
#SBATCH --job-name=sim_unknown_cmp
#SBATCH --output=slurm_log/%x_%j.out
#SBATCH --error=slurm_log/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Unknown Fraction Comparison: Simulated Benchmark
# ============================================================================
# Runs sequentially: first WITH unknown, then WITHOUT unknown
# (avoids checkpoint contamination between parallel runs)
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
INPUT_FOLDER="${REPO_ROOT}/replicates/high_seg"
SLURM_LOG_DIR="${REPO_ROOT}/tests/slurm_log"
mkdir -p "${SLURM_LOG_DIR}"

# Environment Setup
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env
module load gurobi/12.0.3

cd "${REPO_ROOT}"

# ============================================================================
# Run 1: WITH Unknown (default behavior)
# ============================================================================
echo "=========================================="
echo "CONDITION: WITH Unknown fraction"
echo "=========================================="

OUTPUT_WITH="${REPO_ROOT}/test_results/unknown_comparison/with_unknown"
mkdir -p "${OUTPUT_WITH}"

# Clean checkpoints to avoid contamination
rm -f checkpoints/Wu_rep_*_gene_expression_*.npz

python tests/test_citegeist_simulated.py \
    --radius 4 \
    --lambda_reg 1 \
    --alpha_elastic 0.7 \
    --max_y_change 0.4 \
    --input_folder "${INPUT_FOLDER}" \
    --output_folder "${OUTPUT_WITH}" \
    --per-marker-beta \
    --beta-min 0.1 \
    --beta-max 2.0

echo "WITH Unknown completed (exit=$?)"

# ============================================================================
# Run 2: WITHOUT Unknown
# ============================================================================
echo "=========================================="
echo "CONDITION: WITHOUT Unknown fraction"
echo "=========================================="

OUTPUT_WITHOUT="${REPO_ROOT}/test_results/unknown_comparison/no_unknown"
mkdir -p "${OUTPUT_WITHOUT}"

# Clean checkpoints to avoid contamination
rm -f checkpoints/Wu_rep_*_gene_expression_*.npz

python tests/test_citegeist_simulated.py \
    --radius 4 \
    --lambda_reg 1 \
    --alpha_elastic 0.7 \
    --max_y_change 0.4 \
    --input_folder "${INPUT_FOLDER}" \
    --output_folder "${OUTPUT_WITHOUT}" \
    --per-marker-beta \
    --beta-min 0.1 \
    --beta-max 2.0 \
    --no-unknown

echo "WITHOUT Unknown completed (exit=$?)"

echo "=========================================="
echo "Both conditions complete. Compare results in:"
echo "  WITH:    ${OUTPUT_WITH}"
echo "  WITHOUT: ${OUTPUT_WITHOUT}"
echo "=========================================="
