#!/bin/bash
#SBATCH --job-name=sace_benchmark
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --output=slurm-%j.out

# SACE per-cell GEX benchmark on Xenium pseudo-Visium
# Usage: sbatch sbatch_sace_benchmark.sh <region_id>
#   e.g. for region in 0 1 2 3 4; do sbatch sbatch_sace_benchmark.sh $region; done

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

REGION=${1:-0}

module load gurobi/12.0.3

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

export OMP_NUM_THREADS=1

REPO=${CITEGEIST_REPO_ROOT}
OUTPUT_DIR="${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/results_sace"

echo "Running SACE benchmark region=${REGION}"
python "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/src/run_sace_benchmark.py" \
    --region "${REGION}" \
    --output_dir "${OUTPUT_DIR}" \
    --max_iter 10 \
    --n_0 10.0

echo "Done. Results in ${OUTPUT_DIR}/sace_region_${REGION}.json"
