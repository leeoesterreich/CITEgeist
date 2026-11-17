#!/bin/bash
#SBATCH --job-name=tangram_refsens
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --array=0-19

# Array layout (20 jobs):
#   0-4:   high_seg, full ref
#   5-9:   high_seg, downsampled ref
#   10-14: mixed, full ref
#   15-19: mixed, downsampled ref

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../_paths.sh"
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

cd ${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/refsens

module load cuda/12.6.3
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/tangram_env

TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ ${TASK_ID} -lt 5 ]; then
    REP=${TASK_ID}
    CONDITION="high_seg"
    REF_COND="full"
elif [ ${TASK_ID} -lt 10 ]; then
    REP=$((TASK_ID - 5))
    CONDITION="high_seg"
    REF_COND="downsampled"
elif [ ${TASK_ID} -lt 15 ]; then
    REP=$((TASK_ID - 10))
    CONDITION="mixed"
    REF_COND="full"
else
    REP=$((TASK_ID - 15))
    CONDITION="mixed"
    REF_COND="downsampled"
fi

OUTPUT_DIR="./results/Tangram/${CONDITION}/${REF_COND}/tangram_map_${REP}"
mkdir -p "${OUTPUT_DIR}"

START_TIME=$(date +%s)
echo "[$(date)] Tangram refsens: rep=${REP} condition=${CONDITION} ref=${REF_COND}"

python -u run_tangram_refsens.py \
    --replicate-id ${REP} \
    --condition ${CONDITION} \
    --ref-condition ${REF_COND} \
    --ref-dir ./reference \
    --output-dir "${OUTPUT_DIR}"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "[$(date)] Done in ${RUNTIME}s"
