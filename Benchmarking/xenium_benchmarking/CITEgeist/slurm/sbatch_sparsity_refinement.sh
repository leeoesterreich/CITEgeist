#!/bin/bash
#SBATCH --job-name=sparsity_refine
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/results_sparsity_refinement/slurm_%j.log
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/results_sparsity_refinement/slurm_%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1

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

set -u

cd ${CITEGEIST_REPO_ROOT}

# Ensure output directory exists for SLURM logs
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/results_sparsity_refinement

# Parse env vars (defaults: all regions, all original configs)
REGIONS="${REGIONS:-0 1 2 3 4}"
CONFIGS="${CONFIGS:-baseline gex_only refine_only both}"

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_sparsity_refinement.py \
    --regions $REGIONS \
    --configs $CONFIGS

echo "DONE: $(date)"
