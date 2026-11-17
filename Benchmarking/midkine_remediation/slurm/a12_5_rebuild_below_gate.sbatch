#!/bin/bash
#SBATCH --job-name=a12_5_rebuild_below_gate
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_5b_%x_%A_%a.log
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
# A12.5b BELOW-GATE per-cell rebuild (one specimen per array task) via the scoped bypass
# wrapper src/rebuild_below_gate.py. Use ONLY for the two specimens gated OUT of the normal
# A12.5 rebuild because their fullres H&E is >1.0 um/px. CPU job on htc (no GPU, no race).
# Env: SPECIMENS (space-separated, indexed by the array task), RUN_ID, OUT_ROOT, DATA_ROOT.
# The bypass wrapper does not edit rebuild_single_cell.py or protocol.json; it shims the
# gate checks in this process only. Output carries the current protocol's hash (protocol.json is read, not edited).
# Env: MODULE3_DIR (stamped B1 dir) and NUCLEI_DIR (stamped B2 dir) are REQUIRED; B2 covers all 12 specimens, so nothing here reads a v40 r5 path.
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
cd "${CODE_DIR:-${CITEGEIST_REPO_ROOT}}"
MAIN=${CITEGEIST_REPO_ROOT}   # data/outputs stay in the main tree; CODE_DIR (if set) is a clean detached worktree
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
SP_ARR=($SPECIMENS)
SPECIMEN=${SP_ARR[$SLURM_ARRAY_TASK_ID]}
echo "B3 CPU (no segmentation: nuclei from the stamped B2 artifact) specimen=$SPECIMEN out=$OUT_ROOT/$SPECIMEN host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/rebuild_below_gate.py \
  --specimen "$SPECIMEN" --data-root "$DATA_ROOT" --module3-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}" --nuclei-dir "${NUCLEI_DIR:?NUCLEI_DIR (the B2 artifact dir) is required}" \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" --protocol Benchmarking/midkine_remediation/protocol.json
echo "done $SPECIMEN end=$(date -Is)"
