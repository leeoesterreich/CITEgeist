#!/bin/bash
#SBATCH --job-name=train_dino
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --output=slurm/logs/train_dino_%j.out
#SBATCH --error=slurm/logs/train_dino_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00

# train_dino SLURM script

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
set -e

# Clean LD_LIBRARY_PATH to avoid GLIBC conflicts with system libraries
export LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v "ondemand-jupyter" | tr '\n' ':' | sed 's/:$//')

echo "============================================"
echo "DINO Training Job"
echo "Started: $(date)"
echo "============================================"

eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
BENCHMARK_DIR="${REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist"
PATCHES_DIR="${BENCHMARK_DIR}/output/vae_masked/patches_combined"
OUTPUT_DIR="${BENCHMARK_DIR}/output/dino_ssl"

cd "${BENCHMARK_DIR}"
mkdir -p slurm/logs

python src/train_dino.py \
    --patches-dir "${PATCHES_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --epochs 300 \
    --batch-size 128 \
    --lr 5e-4 \
    --device cuda \
    --checkpoint-interval 50 \
    --num-workers 8

echo "Completed: $(date)"
