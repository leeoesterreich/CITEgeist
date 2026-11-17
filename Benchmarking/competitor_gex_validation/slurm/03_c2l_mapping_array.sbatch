#!/bin/bash
#SBATCH --job-name=c2l_map
#SBATCH --output=Benchmarking/competitor_gex_validation/logs/03_c2l_mapping_%A_%a.out
#SBATCH --error=Benchmarking/competitor_gex_validation/logs/03_c2l_mapping_%A_%a.err
#SBATCH --array=0-11
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --gres=gpu:1

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../_paths.sh"
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

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
cd "${REPO_ROOT}"
export MPLCONFIGDIR=/tmp/mpl_${SLURM_JOB_ID}
export HDF5_USE_FILE_LOCKING=FALSE

module load cuda/12.6

# Read sample name from samples.txt (0-indexed SLURM_ARRAY_TASK_ID → 1-indexed sed)
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" Benchmarking/competitor_gex_validation/slurm/samples.txt)

echo "=== C2L Spatial Mapping ==="
echo "Job ID: ${SLURM_JOB_ID}, Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Sample: ${SAMPLE}"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "race namespace=${GPU_RACE_NAMESPACE:-none} partition=${GPU_RACE_PARTITION:-none} head=$(git rev-parse HEAD)"
echo "Start: $(date)"

"$(cg_env_python cell2location_env)" -u Benchmarking/competitor_gex_validation/run_c2l_mapping.py --sample "${SAMPLE}"

echo "=== Done: $(date) ==="
