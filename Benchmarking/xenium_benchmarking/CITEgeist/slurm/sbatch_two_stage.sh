#!/bin/bash
#SBATCH --job-name=two_stage
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/two_stage_%A_%a.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs/two_stage_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# Load modules
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

# Activate environment
eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

# Navigate to repo
cd ${CITEGEIST_REPO_ROOT}

# Run benchmark for this region
REGION=${SLURM_ARRAY_TASK_ID}
OUTPUT_DIR="Benchmarking/xenium_benchmarking/CITEgeist/output/two_stage"

echo "=== Job Info ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Region: ${REGION}"
echo "Start time: $(date)"

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_two_stage.py \
    --region ${REGION} \
    --output-dir ${OUTPUT_DIR} \
    --device cuda \
    --n-epochs 50

echo "End time: $(date)"
