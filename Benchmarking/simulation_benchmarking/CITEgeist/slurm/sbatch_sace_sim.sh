#!/bin/bash
#SBATCH --job-name=sace_sim_gex
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --gres=gpu:1
#SBATCH --array=0-4

# GPU_RACE_READONLY_OK=1
# Usage: sbatch sbatch_sace_sim.sh <condition>
# Example: sbatch sbatch_sace_sim.sh high_seg

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

CONDITION="${SACE_CONDITION:-${1:?Usage: SACE_CONDITION=high_seg gpu-race-submit.sh sbatch_sace_sim.sh}}"

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

module load gurobi/12.0.3

SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CITEgeist/src"

echo "$(date): Starting SACE benchmark: condition=$CONDITION rep=$SLURM_ARRAY_TASK_ID"

python "$SCRIPT_DIR/benchmark_sace_simulation.py" \
    --replicate "$SLURM_ARRAY_TASK_ID" \
    --condition "$CONDITION"

echo "$(date): SACE benchmark complete: condition=$CONDITION rep=$SLURM_ARRAY_TASK_ID"
