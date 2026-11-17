#!/bin/bash
#SBATCH --job-name=stage2_morph
#SBATCH --output=logs/stage2_morph_%A_%a.out
#SBATCH --error=logs/stage2_morph_%A_%a.err
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --array=0-4

# Two-Stage Morphology Benchmark
# Stage 1: Hybrid proportions (pre-computed)
# Stage 2: Morphology-based single-cell assignment

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
set -e

# Setup environment
eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

# Paths
REPO_ROOT="${CITEGEIST_REPO_ROOT}"
SCRIPT="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage_morphology.py"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output/two_stage_morphology"

# Create output directory
mkdir -p "${OUTPUT_DIR}"
SLURM_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/slurm"
mkdir -p "${SLURM_DIR}/logs"

# Run benchmark for this region
REGION=${SLURM_ARRAY_TASK_ID}

echo "=========================================="
echo "Two-Stage Morphology Benchmark"
echo "Region: ${REGION}"
echo "=========================================="

python "${SCRIPT}" \
    --region "${REGION}" \
    --output-dir "${OUTPUT_DIR}"

echo "Done!"
