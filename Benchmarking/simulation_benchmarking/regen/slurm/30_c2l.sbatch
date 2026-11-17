#!/bin/bash
#SBATCH --job-name=regen_c2l
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
# A10 Step 2 (C2L): reconstruct from saved model.pt (default) or retrain (C2L_RETRAIN=1, pre-authorised fallback).
# Submit: gpu-race-submit.sh regen/slurm/30_c2l.sbatch --array=0-9   (array task t<5 high_seg rep t, else mixed rep t-5)
# Race-safe: writes under results/c2l/.ns_$GPU_RACE_NAMESPACE, then claims the final dir with an atomic mv.
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
TID=${SLURM_ARRAY_TASK_ID:?}
if [ "$TID" -lt 5 ]; then COND=high_seg; REP=$TID; else COND=mixed; REP=$((TID-5)); fi
NS=${GPU_RACE_NAMESPACE:-${GPU_RACE_PARTITION:-default}}
echo "race namespace/partition: $NS  host: $(hostname)  gpu: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader)"
ROOT=$B/regen/results/c2l${C2L_RETRAIN:+_retrain}
TMPOUT=$ROOT/.ns_$NS/$COND/Wu_rep_$REP
FINAL=$ROOT/$COND/Wu_rep_$REP
mkdir -p $(dirname $TMPOUT) $(dirname $FINAL); rm -rf $TMPOUT
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/cell2location_env
cd ${CITEGEIST_REPO_ROOT}
T0=$(date +%s)
python -u $B/regen/src/c2l_reconstruct.py $COND $REP $TMPOUT --seed 0 ${C2L_RETRAIN:+--retrain}
echo "RUNTIME_SEC $(($(date +%s)-T0)) method=C2L cond=$COND rep=$REP ns=$NS"
if [ -e $FINAL ]; then echo "final exists, keeping first winner"; else mv $TMPOUT $FINAL; fi
