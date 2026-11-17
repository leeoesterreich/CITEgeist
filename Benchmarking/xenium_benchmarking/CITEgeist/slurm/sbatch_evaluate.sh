#!/bin/bash
#SBATCH --job-name=citegeist_evaluate
#SBATCH --cluster=htc
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=slurm/logs/evaluate_%j.out
#SBATCH --error=slurm/logs/evaluate_%j.err

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

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env
set -u  # Re-enable nounset after conda activation

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
EVAL_DIR="${SCRIPT_DIR}/../evaluation/src"
OUTPUT_6TYPE="${SCRIPT_DIR}/output_singler_6type"
OUTPUT_7TYPE="${SCRIPT_DIR}/output_singler_7type"

# Poll for completion (cross-cluster dependency doesn't work)
echo "Waiting for benchmark jobs to complete..."
for region_id in 0 1 2 3 4; do
    marker_6="${OUTPUT_6TYPE}/Xenium_region_${region_id}/run_stats.json"
    marker_7="${OUTPUT_7TYPE}/Xenium_region_${region_id}/run_stats.json"
    attempts=0
    while [[ ! -f "${marker_6}" || ! -f "${marker_7}" ]]; do
        if [[ ${attempts} -ge 120 ]]; then
            echo "ERROR: Timeout waiting for region ${region_id} (waited 2 hours)"
            exit 1
        fi
        sleep 60
        attempts=$((attempts + 1))
    done
    echo "Region ${region_id} complete"
done

echo "All regions complete. Running evaluation at $(date)"

echo "=== 6-type evaluation ==="
python "${EVAL_DIR}/compare_all_methods.py" --gt-type singler_6type

echo "=== 7-type evaluation ==="
python "${EVAL_DIR}/compare_all_methods.py" --gt-type singler_7type

echo "Evaluation complete at $(date)"
