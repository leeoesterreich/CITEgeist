#!/bin/bash
#SBATCH --job-name=sc_res_r3
#SBATCH --output=logs/sc_resolution_r3_%j.out
#SBATCH --error=logs/sc_resolution_r3_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s

# Single-cell resolution benchmark for region 3 with increased memory (64GB)
# Region 3 has a larger image that requires more RAM for Cellpose segmentation

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

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py \
    --region 3 \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution \
    --use-gpu \
    --spot-diameter-um 55.0
