#!/bin/bash
#SBATCH --job-name=v2_b18_celltype_map
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/b18_celltype_map_%j.log
# v2 B18 (run id from v42_run_ids.yaml, currently v2b18): which cell types hold the MDK-NCL sender and receiver COMMOT score. Descriptive only.
# Promoted from the exploratory check 3; same math. Writes a NEW dir results/<protocol_sha12>/<RUN_ID> (refuses to overwrite).
# Required env (no defaults; set by .commands/submit-v2-celltype-map.sh via chain_preflight.py): COMMOT_ROOT (the B6 run dir), HURDLE_RUN (the B8 run dir),
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
PY="$(cg_env_python CITEgeist_env)"
: "${COMMOT_ROOT:?COMMOT_ROOT (the B6 run dir) is required}" "${HURDLE_RUN:?HURDLE_RUN (the B8 run dir) is required}" "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B18 id) is required}"
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
[ ! -e "$OUT_ROOT" ] || { echo "refusing: $OUT_ROOT already exists (a run dir is never overwritten)" >&2; exit 1; }
"$PY" -u Benchmarking/midkine_remediation/celltype_map.py \
  --commot-run "$COMMOT_ROOT" --hurdle-run "$HURDLE_RUN" --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --protocol Benchmarking/midkine_remediation/protocol.json
echo "done $RUN_ID end=$(date -Is)"
