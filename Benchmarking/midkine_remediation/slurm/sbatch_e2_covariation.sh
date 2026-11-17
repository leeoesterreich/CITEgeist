#!/bin/bash
#SBATCH --job-name=v2_b21_e2
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/b21_e2_%j.log
# v2 B21 (run id from v42_run_ids.yaml, currently v2b21): keratin presence and the estrogen programme (check 7), ESR1 regulon activity vs MDK (check 8). Descriptive only.
# Promoted from the exploratory checks 7 and 8; same math and seeds. Writes a NEW dir results/<protocol_sha12>/<RUN_ID> (refuses to overwrite).
# Required env (no defaults; set by .commands/submit-v2-e2.sh via chain_preflight.py): COMP_RUN (the B5 run dir), COMMOT_ROOT (the B6 run dir), HURDLE_RUN (the B8 run dir),
# OUT_ROOT, RUN_ID, V2_NAMESPACE_ROOT.
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
cd "${CODE_DIR:-${CITEGEIST_REPO_ROOT}}"
MAIN=${CITEGEIST_REPO_ROOT}
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
PY="$(cg_env_python CITEgeist_env)"
: "${COMP_RUN:?COMP_RUN (the B5 run dir) is required}" "${COMMOT_ROOT:?COMMOT_ROOT (the B6 run dir) is required}" "${HURDLE_RUN:?HURDLE_RUN (the B8 run dir) is required}" "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B21 id) is required}"
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
[ ! -e "$OUT_ROOT" ] || { echo "refusing: $OUT_ROOT already exists (a run dir is never overwritten)" >&2; exit 1; }
"$PY" -u Benchmarking/midkine_remediation/e2_covariation.py \
  --compartments-run "$COMP_RUN" --commot-root "$COMMOT_ROOT" --hurdle-run "$HURDLE_RUN" --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --protocol Benchmarking/midkine_remediation/protocol.json
echo "done $RUN_ID end=$(date -Is)"
