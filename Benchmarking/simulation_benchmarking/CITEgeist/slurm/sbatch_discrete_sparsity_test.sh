#!/bin/bash
#SBATCH --job-name=sparsity_test
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/sparsity_test_%A_%a.out
#SBATCH --error=logs/sparsity_test_%A_%a.err

# Run discrete cell assignment benchmark with sparsity penalty
# Tests lambda_sparse=0.1 to validate Phase 0 sparsity fix
# Usage: sbatch sbatch_discrete_sparsity_test.sh high_seg dapi
#        sbatch sbatch_discrete_sparsity_test.sh mixed dapi

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
CONDITION=${1:-high_seg}
MODE=${2:-dapi}
LAMBDA_SPARSE=${3:-0.1}
REPLICATE_ID=$SLURM_ARRAY_TASK_ID

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

cd ${CITEGEIST_REPO_ROOT}

OUTPUT_DIR=Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_sparsity_test

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id $REPLICATE_ID \
    --condition $CONDITION \
    --mode $MODE \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20 \
    --lambda-sparse $LAMBDA_SPARSE

echo "Completed replicate $REPLICATE_ID, condition $CONDITION, mode $MODE, lambda_sparse=$LAMBDA_SPARSE"
