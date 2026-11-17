#!/bin/bash
#SBATCH --job-name=qp_sim_clean
#SBATCH --array=0-9
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --gres=gpu:1
# GPU_RACE_READONLY_OK=1

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

# Map SLURM_ARRAY_TASK_ID (0-9) to (condition, replicate_id)
# IDs 0-4: high_seg reps 0-4
# IDs 5-9: mixed reps 0-4
if [ "${SLURM_ARRAY_TASK_ID}" -lt 5 ]; then
    CONDITION="high_seg"
    REP_ID="${SLURM_ARRAY_TASK_ID}"
else
    CONDITION="mixed"
    REP_ID=$((SLURM_ARRAY_TASK_ID - 5))
fi

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
OUTPUT_DIR="${REPO_ROOT}/Benchmarking/simulation_benchmarking/CITEgeist/results_qp"
ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/cuopt_test"

mkdir -p "${OUTPUT_DIR}"

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate "${ENV}"

export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

echo "=== QP Clean Sim: condition=${CONDITION} rep=${REP_ID} ==="
echo "Start: $(date)"
echo "Node: $(hostname)"
nvidia-smi --query-gpu=name --format=csv,noheader || true

cd "${REPO_ROOT}"

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_simulation_unified.py \
    --replicate-id "${REP_ID}" \
    --condition "${CONDITION}" \
    --mode dapi \
    --output-dir "${OUTPUT_DIR}"

echo "=== Done: condition=${CONDITION} rep=${REP_ID} exit=$? $(date) ==="
