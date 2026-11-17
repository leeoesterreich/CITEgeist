#!/bin/bash
#SBATCH --job-name=gs_sace
#SBATCH --cluster=htc
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# REGION_ID, OFFSET, PROPS_DIR passed via --export
# COMBO_ID = OFFSET + SLURM_ARRAY_TASK_ID

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
COMBO_ID=$((OFFSET + SLURM_ARRAY_TASK_ID))
SCRIPT_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/gridsearch"

# Redirect output (can't use ${REGION_ID} in #SBATCH --output)
mkdir -p "${SCRIPT_DIR}/results/sace/logs"
exec > "${SCRIPT_DIR}/results/sace/logs/combo_${COMBO_ID}_region_${REGION_ID}_${SLURM_JOB_ID}.out" 2>&1

module load python/ondemand-jupyter-python3.8 2>/dev/null
source /software/rhel9/manual/install/python/ondemand-jupyter-python3.8-2021.05/etc/profile.d/conda.sh
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

echo "combo_id=${COMBO_ID} region_id=${REGION_ID} array_task=${SLURM_ARRAY_TASK_ID} offset=${OFFSET}"

python -u "${SCRIPT_DIR}/run_sace_gridsearch_v3.py" \
    --combo-id "${COMBO_ID}" \
    --region-id "${REGION_ID}" \
    --props-dir "${PROPS_DIR}"
