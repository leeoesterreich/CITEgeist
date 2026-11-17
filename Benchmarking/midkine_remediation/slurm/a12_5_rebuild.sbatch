#!/bin/bash
#SBATCH --job-name=a12_5_rebuild
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_5_%x_%A_%a.log
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
# v2 B3 Hungarian per-cell rebuild (one specimen per array task). CPU job on htc (no GPU, no race): nuclei come from the stamped B2 artifact; see .commands/submit-a12-rebuild.sh. Memory 48G: max RSS of the v2b3 GPU jobs (which also ran StarDist) was 31.4G (P1-S2, 88k cells).
# Env: MODULE3_DIR (stamped B1 dir) and NUCLEI_DIR (stamped B2 dir) are REQUIRED, V2_NAMESPACE_ROOT (main checkout), CODE_DIR (worktree to run from; default main tree), SPECIMENS (space-separated, indexed by the array task), RUN_ID, OUT_ROOT (= output/morphology_assignment_v4_hungarian/<protocol_sha12>/<RUN_ID>), DATA_ROOT.
# The rebuild writes <OUT_ROOT>/<specimen>/ via temp dir + atomic rename and refuses a non-empty target.
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
"$PY" -u Benchmarking/midkine_remediation/src/rebuild_single_cell.py \
  --specimen "$SPECIMEN" --data-root "$DATA_ROOT" --module3-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}" --nuclei-dir "${NUCLEI_DIR:?NUCLEI_DIR (the B2 artifact dir) is required}" \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" --protocol Benchmarking/midkine_remediation/protocol.json
echo "done $SPECIMEN end=$(date -Is)"
