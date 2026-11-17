#!/bin/bash
#SBATCH --job-name=sace_abl_var
#SBATCH --cluster=htc
#SBATCH --array=0-29
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%A_%a.out

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
REPO_ROOT="${CITEGEIST_REPO_ROOT}"
CONDA_ENV="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env"

VARIANTS=(marker_i1_lr prop_i1_cont prop_i1_lr)
VARIANT_IDX=$(( SLURM_ARRAY_TASK_ID / 10 ))
REP_WITHIN=$(( SLURM_ARRAY_TASK_ID % 10 ))

VARIANT=${VARIANTS[$VARIANT_IDX]}

if [ $REP_WITHIN -lt 5 ]; then
    CONDITION="high_seg"
    REPLICATE=$REP_WITHIN
else
    CONDITION="mixed"
    REPLICATE=$(( REP_WITHIN - 5 ))
fi

source activate "${CONDA_ENV}" 2>/dev/null || conda activate "${CONDA_ENV}"
cd "${REPO_ROOT}/CITEgeist"

python -u ../Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_sace_simulation.py \
    --variant "$VARIANT" \
    --replicate "$REPLICATE" \
    --condition "$CONDITION"
