#!/bin/bash
#SBATCH --job-name=sace_xenium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --gres=gpu:1
#SBATCH --array=0-4

# GPU_RACE_READONLY_OK=1
# CITEgeist-SACE Xenium GEX benchmark (own QP proportions, 7 cell types)
#
# Submit via gpu-race-submit.sh (DO NOT add --partition or --cluster here):
#   cd Benchmarking/xenium_benchmarking/CITEgeist/slurm
#   gpu-race-submit.sh sbatch_sace_xenium.sh
#
# To run a single region manually (overrides array):
#   SACE_REGION=2 gpu-race-submit.sh sbatch_sace_xenium.sh

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
set -eo pipefail

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

export OMP_NUM_THREADS=1

# Allow SACE_REGION env var to override array task ID (for single-region testing)
REGION=${SACE_REGION:-${SLURM_ARRAY_TASK_ID}}

REPO=${CITEGEIST_REPO_ROOT}
OUTPUT_DIR="${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/results_sace_singler_7type"

echo "=== CITEgeist-SACE Xenium benchmark ==="
echo "Region:    ${REGION}"
echo "Date:      $(date)"
echo "Host:      $(hostname)"
echo "GPU:       $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'none')"
echo "OutputDir: ${OUTPUT_DIR}"

mkdir -p "$(dirname "${OUTPUT_DIR}")"

python "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_sace_xenium.py" \
    --region "${REGION}" \
    --output_dir "${OUTPUT_DIR}"

echo "=== Done. Output at ${OUTPUT_DIR}/Xenium_region_${REGION}/layers/ ==="
echo "Date: $(date)"
