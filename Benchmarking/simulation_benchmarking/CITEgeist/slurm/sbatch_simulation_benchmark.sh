#!/bin/bash
#SBATCH --job-name=sim_bench
#SBATCH --partition=htc
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-29

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

module load gurobi/12.0.3

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

TASK_ID=${SLURM_ARRAY_TASK_ID}

# Decode array index to (condition, mode, replicate)
# 0-4:  high_seg/dapi/rep0-4
# 5-9:  high_seg/h_and_e/rep0-4
# 10-14: high_seg/gt_ablation/rep0-4
# 15-19: mixed/dapi/rep0-4
# 20-24: mixed/h_and_e/rep0-4
# 25-29: mixed/gt_ablation/rep0-4
CONDITIONS=(high_seg high_seg high_seg mixed mixed mixed)
MODES=(dapi h_and_e gt_ablation dapi h_and_e gt_ablation)

BLOCK=$((TASK_ID / 5))
REP=$((TASK_ID % 5))
CONDITION=${CONDITIONS[$BLOCK]}
MODE=${MODES[$BLOCK]}

echo "=== Task ${TASK_ID}: condition=${CONDITION} mode=${MODE} replicate=${REP} ==="

BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CITEgeist"
SCRIPT_DIR="${BASE_DIR}/src"
RESULTS_DIR="${BASE_DIR}/results"

python -u "${SCRIPT_DIR}/benchmark_simulation_unified.py" \
    --replicate-id "${REP}" \
    --condition "${CONDITION}" \
    --mode "${MODE}" \
    --output-dir "${RESULTS_DIR}"

echo "=== Done ==="
