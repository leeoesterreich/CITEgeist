#!/bin/bash
#SBATCH --job-name=regen_card
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-9
#SBATCH --time=04:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
# A10 Step 3: rerun card with the ORIGINAL per-condition script and settings; outputs only under regen/results.
# Array task t: t<5 -> high_seg rep t ; t>=5 -> mixed rep t-5.
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
B=${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking
if [ "$SLURM_ARRAY_TASK_ID" -lt 5 ]; then COND=high_seg; REP=$SLURM_ARRAY_TASK_ID; else COND=mixed; REP=$((SLURM_ARRAY_TASK_ID-5)); fi
OUT=$B/regen/results/card/$COND/Wu_rep_$REP
mkdir -p $OUT
eval "$(conda shell.bash hook 2>/dev/null)"
T0=$(date +%s)
export PATH=${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CARD_env/bin:$PATH
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CARD_env
test -s $B/CARD/reference_csv/reference_counts.csv
cd $B/CARD/$COND
Rscript --vanilla CARD_pipeline_${COND}.R --replicates $REP --output_dir $OUT --mode reference
echo "RUNTIME_SEC $(( $(date +%s) - T0 )) method=CARD cond=$COND rep=$REP"
