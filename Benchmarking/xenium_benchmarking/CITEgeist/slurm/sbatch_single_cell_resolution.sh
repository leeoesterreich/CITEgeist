#!/bin/bash
#SBATCH --job-name=sc_resolution
#SBATCH --output=logs/sc_resolution_%A_%a.out
#SBATCH --error=logs/sc_resolution_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s

# Single-cell resolution benchmark using Module 3b nucleus assignment pipeline.
#
# This script tests the new single-cell resolution features:
# 1. Loads existing hybrid results (proportions + deconvolved GEX)
# 2. Runs Cellpose to get full nuclei masks (not just counts)
# 3. Maps nuclei to spots
# 4. Runs run_nucleus_assignment() pipeline (morphology features + soft-label classifier + Hungarian assignment)
# 5. Distributes GEX to individual cells
# 6. Creates single-cell AnnData output
#
# Requires: Existing hybrid_cellpose results for each region

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

cd ${CITEGEIST_REPO_ROOT}

# Check that hybrid results exist
HYBRID_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/hybrid_cellpose"
REGION_DIR="${HYBRID_DIR}/Xenium_region_${SLURM_ARRAY_TASK_ID}"
if [ ! -d "$REGION_DIR" ]; then
    echo "ERROR: Hybrid results not found at $REGION_DIR"
    echo "Please run sbatch_hybrid_cellpose_gpu.sh first"
    exit 1
fi

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py \
    --region ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution \
    --use-gpu \
    --spot-diameter-um 55.0
