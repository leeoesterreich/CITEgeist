#!/bin/bash
#SBATCH --job-name=v2_nuclei
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/v2_nuclei_%A_%a.log
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
# v2 breast chain step 2 (stage B2): nuclei counts for EXACTLY the p0 spots, stamped. One specimen per array task. Submit via gpu-race-submit.sh
# (.commands/submit-v2-breast-chain.sh nuclei). Env: CODE_DIR, SPECIMENS (space-separated, indexed by the array task), P0_DIR, NUCLEI_DIR (NEW <root>/v2b2).
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
set -euo pipefail
umask 022   # outputs readable by the lab
cd "${CODE_DIR:?CODE_DIR required}"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
SP_ARR=($SPECIMENS)
SPECIMEN=${SP_ARR[$SLURM_ARRAY_TASK_ID]}
# race namespace line (write-safety gate); the output dir is per specimen, refused if non-empty
echo "race namespace=${GPU_RACE_NAMESPACE:-none} partition=${GPU_RACE_PARTITION:-none} specimen=$SPECIMEN out=$NUCLEI_DIR/$SPECIMEN host=$(hostname) start=$(date -Is)"
nvidia-smi -L || true
"$PY" -u Benchmarking/midkine_remediation/src/v2_inputs.py nuclei --specimen "$SPECIMEN" --p0-dir "${P0_DIR:?P0_DIR required}" --nuclei-dir "${NUCLEI_DIR:?NUCLEI_DIR required}"
echo "done $SPECIMEN end=$(date -Is)"
