#!/bin/bash
#SBATCH --job-name=nb_regression
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --cluster=htc

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

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"

source activate "$ENV"
module load gurobi/12.0.3

export OMP_NUM_THREADS=1

cd "$REPO_ROOT"

CONDITION="${1:-high_seg}"
REPS="${@:2}"
if [ -z "$REPS" ]; then
    REPS="0 1 2 3 4"
fi

echo "=== NB Regression Benchmark ==="
echo "Condition: $CONDITION"
echo "Reps: $REPS"
echo "Date: $(date)"
echo ""

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_nb_regression.py \
    --condition "$CONDITION" \
    --reps $REPS

echo "Exit code: $?"
