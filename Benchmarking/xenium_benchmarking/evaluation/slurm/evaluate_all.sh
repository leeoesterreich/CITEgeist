#!/bin/bash
#SBATCH --job-name=xenium_eval            # Job name
#SBATCH --output=slurm_log/%x_%j.out      # Standard output
#SBATCH --error=slurm_log/%x_%j.err       # Standard error
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=2                 # CPU cores per task
#SBATCH --mem=16G                         # Memory per node
#SBATCH --time=01:00:00                   # Time limit
#SBATCH --cluster=htc                     # Cluster
#SBATCH --partition=htc                   # Partition

# ============================================================================
# Evaluate Xenium Benchmark Results
# ============================================================================
# Run this script AFTER all regions have been processed by xenium_benchmark.sh
# ============================================================================

# Directories
REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
PSEUDOVISIUM_DIR="${REPO_ROOT}/Benchmarking/xenium_pseudovisium"
BENCH_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking"
CITEGEIST_DIR="${BENCH_DIR}/CITEgeist"
EVALUATION_DIR="${BENCH_DIR}/evaluation"
INPUT_DIR="${PSEUDOVISIUM_DIR}/data"
OUTPUT_DIR="${CITEGEIST_DIR}/output"
METRICS_DIR="${BENCH_DIR}/metrics"

# Create metrics directory
mkdir -p "${METRICS_DIR}"

# Activate environment
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO_ROOT}"

echo "============================================"
echo "Evaluating Xenium Benchmark Results"
echo "============================================"

# Run evaluation for all regions
python -c "
import sys
import logging
logging.basicConfig(level=logging.INFO)

sys.path.insert(0, 'Benchmarking/xenium_benchmarking/evaluation/src')
from evaluate_metrics import evaluate_all_regions

results = evaluate_all_regions(
    input_dir='${INPUT_DIR}',
    output_dir='${METRICS_DIR}',
    predictions_dir='${OUTPUT_DIR}',
    n_regions=5,
    prefix='Xenium',
)

print('\n' + '='*60)
print('FINAL RESULTS')
print('='*60)
print(f'Regions: {results[\"n_regions\"]}')
print(f'Total spots: {results[\"total_spots\"]}')
print(f'Mean JSD: {results[\"JSD_median_mean\"]:.4f} ± {results[\"JSD_median_std\"]:.4f}')
print(f'Mean RMSE: {results[\"RMSE_global_mean\"]:.4f} ± {results[\"RMSE_global_std\"]:.4f}')
print(f'Mean MAE: {results[\"MAE_global_mean\"]:.4f} ± {results[\"MAE_global_std\"]:.4f}')
print(f'Mean Pearson r: {results[\"Pearson_r_mean\"]:.4f} ± {results[\"Pearson_r_std\"]:.4f}')
print('='*60)
"

echo ""
echo "Metrics saved to: ${METRICS_DIR}"
echo "============================================"
