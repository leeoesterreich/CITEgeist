#!/bin/bash
#SBATCH --job-name=d1_card
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-9
#SBATCH --time=12:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
# A10 diagnostic D1: full-gene ER+ reference; fixed script copies, only the reference changes.
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../../../_paths.sh"
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
D=${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/regen/diag/d1
ENVS=${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}
eval "$(conda shell.bash hook 2>/dev/null)"
if [ "$SLURM_ARRAY_TASK_ID" -lt 5 ]; then COND=high_seg; REP=$SLURM_ARRAY_TASK_ID; else COND=mixed; REP=$((SLURM_ARRAY_TASK_ID-5)); fi
OUT=$D/results/card/$COND/Wu_rep_$REP; mkdir -p $OUT; T0=$(date +%s)
export PATH=$ENVS/CARD_env/bin:$PATH
conda activate $ENVS/CARD_env
test -s $D/reference/reference_csv/reference_counts.csv
cd $D/src
Rscript --vanilla CARD_fullgene_${COND}.R --replicates $REP --output_dir $OUT --mode reference
echo "RUNTIME_SEC $(( $(date +%s) - T0 )) method=CARD cond=$COND rep=$REP"
