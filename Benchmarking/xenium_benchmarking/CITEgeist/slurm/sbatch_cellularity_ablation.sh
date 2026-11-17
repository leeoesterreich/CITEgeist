#!/bin/bash
# GPU_RACE_READONLY_OK=1
#SBATCH --job-name=cell_abl
#SBATCH --output=output/cellularity_ablation/slurm_%A_%a.out
#SBATCH --error=output/cellularity_ablation/slurm_%A_%a.err
#SBATCH --array=0-24%3
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
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

module purge
eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

cd ${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist

# Decode array task: 5 regions x 5 conditions
REGION=$((SLURM_ARRAY_TASK_ID / 5))
CONDITION=$((SLURM_ARRAY_TASK_ID % 5))

echo "Task ${SLURM_ARRAY_TASK_ID}: region=${REGION}, condition=${CONDITION}"
echo "$(date) conda done, starting python..."

python -u src/benchmark_cellularity_ablation.py \
    --region ${REGION} \
    --condition ${CONDITION} \
    --output-dir output/cellularity_ablation
