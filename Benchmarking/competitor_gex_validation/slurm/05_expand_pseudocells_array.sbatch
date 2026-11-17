#!/bin/bash
#SBATCH --job-name=expand_pseudo
#SBATCH --output=Benchmarking/competitor_gex_validation/logs/05_expand_pseudocells_%A_%a.out
#SBATCH --error=Benchmarking/competitor_gex_validation/logs/05_expand_pseudocells_%A_%a.err
#SBATCH --partition=htc
#SBATCH --cluster=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --array=0-23

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
set -eo pipefail

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
cd "${REPO_ROOT}"

source activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

# ---- Map array index to method + sample ----
SAMPLES=(
    "HCC22-088-P1-S1" "HCC22-088-P1-S2"
    "HCC22-088-P2-S1" "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A" "HCC22-088-P3-S2"
    "HCC22-088-P4-S1" "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1" "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1" "HCC22-088-P6-S2_D"
)
N_SAMPLES=${#SAMPLES[@]}

TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ ${TASK_ID} -lt ${N_SAMPLES} ]; then
    METHOD="cell2location"
    SAMPLE_IDX=${TASK_ID}
else
    METHOD="tangram"
    SAMPLE_IDX=$((TASK_ID - N_SAMPLES))
fi

SAMPLE=${SAMPLES[${SAMPLE_IDX}]}

echo "=== Expand pseudo-cells ==="
echo "Job ID: ${SLURM_JOB_ID}, Array Task: ${TASK_ID}"
echo "Method: ${METHOD}, Sample: ${SAMPLE}"
echo "Node: $(hostname)"
echo "Start: $(date)"

python -u Benchmarking/competitor_gex_validation/expand_pseudocells.py \
    --method "${METHOD}" \
    --sample "${SAMPLE}"

echo "=== Done: $(date) ==="
