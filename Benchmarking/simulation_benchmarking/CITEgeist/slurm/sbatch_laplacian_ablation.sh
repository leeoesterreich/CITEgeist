#!/bin/bash
#SBATCH --job-name=laplacian_ablation
#SBATCH --output=logs/laplacian_ablation_%A_%a.out
#SBATCH --error=logs/laplacian_ablation_%A_%a.err
#SBATCH --array=0-9
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc

# Array tasks 0-4: lambda_laplacian=0.0 (NO spatial smoothing)
# Array tasks 5-9: lambda_laplacian=0.1 (WITH spatial smoothing)

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

if [ ${SLURM_ARRAY_TASK_ID} -lt 5 ]; then
    REPLICATE_ID=${SLURM_ARRAY_TASK_ID}
    LAMBDA=0.0
else
    REPLICATE_ID=$((SLURM_ARRAY_TASK_ID - 5))
    LAMBDA=0.1
fi

echo "Running replicate ${REPLICATE_ID} with lambda_laplacian=${LAMBDA}"

python Benchmarking/simulation_benchmarking/CITEgeist/src/test_laplacian_ablation.py \
    --replicate-id ${REPLICATE_ID} \
    --condition mixed \
    --lambda-laplacian ${LAMBDA}
