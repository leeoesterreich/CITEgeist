#!/bin/bash
#SBATCH --job-name=gen_sim_images
#SBATCH --partition=htc
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-19
#SBATCH --output=slurm/logs/gen_images_%A_%a.out
#SBATCH --error=slurm/logs/gen_images_%A_%a.err

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
eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

TASK_ID=${SLURM_ARRAY_TASK_ID}

# 20 configs: 2 conditions x 2 modes x 5 replicates
# 0-4:  high_seg/dapi/rep0-4
# 5-9:  high_seg/h_and_e/rep0-4
# 10-14: mixed/dapi/rep0-4
# 15-19: mixed/h_and_e/rep0-4
CONDITIONS=(high_seg high_seg mixed mixed)
MODES=(dapi h_and_e dapi h_and_e)

BLOCK=$((TASK_ID / 5))
REP=$((TASK_ID % 5))
CONDITION=${CONDITIONS[$BLOCK]}
MODE=${MODES[$BLOCK]}

echo "=== Task ${TASK_ID}: Generating ${CONDITION}/${MODE}/rep${REP} ==="

BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CITEgeist"
SCRIPT_DIR="${BASE_DIR}/src"
OUTPUT_DIR="${BASE_DIR}"

python "${SCRIPT_DIR}/generate_simulation_images.py" \
    --replicate-id "${REP}" \
    --condition "${CONDITION}" \
    --mode "${MODE}" \
    --output-dir "${OUTPUT_DIR}" \
    --seed $((42 + REP))

echo "=== Done: ${CONDITION}/${MODE}/rep${REP} ==="
