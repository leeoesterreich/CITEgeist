#!/bin/bash
#SBATCH --job-name=eval_vim_epi
#SBATCH --output=slurm_log/%x_%j.out
#SBATCH --error=slurm_log/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# ============================================================================
# Evaluate VIM-Epithelial experiment vs baseline
# ============================================================================
# Run AFTER all 5 benchmark regions complete.
# Compares baseline achievable-7 with VIM-Epithelial variant.
# ============================================================================

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
CITEGEIST_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
GT_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium/data_granular_gt"
SLURM_LOG_DIR="${CITEGEIST_DIR}/slurm/slurm_log"

# Try both possible baseline locations
BASELINE_DIR="${CITEGEIST_DIR}/output/manual"
if [ ! -d "${BASELINE_DIR}" ] || [ -z "$(ls -A ${BASELINE_DIR} 2>/dev/null)" ]; then
    BASELINE_DIR="${CITEGEIST_DIR}/output_achievable_7"
fi

EXPERIMENT_DIR="${CITEGEIST_DIR}/output/vim_epithelial"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/evaluation/results_vim_epithelial"

mkdir -p "${SLURM_LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"

echo "=============================================="
echo "VIM-Epithelial Experiment Evaluation"
echo "=============================================="
echo "Baseline: ${BASELINE_DIR}"
echo "Experiment: ${EXPERIMENT_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Environment setup
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

python "${CITEGEIST_DIR}/src/evaluate_vim_epithelial.py" \
    --baseline-dir "${BASELINE_DIR}" \
    --experiment-dir "${EXPERIMENT_DIR}" \
    --gt-dir "${GT_DIR}" \
    --n-regions 5 \
    --output-dir "${OUTPUT_DIR}"

echo ""
echo "=============================================="
echo "Evaluation Complete"
echo "=============================================="
echo "Results: ${OUTPUT_DIR}/"
ls -la "${OUTPUT_DIR}/"
