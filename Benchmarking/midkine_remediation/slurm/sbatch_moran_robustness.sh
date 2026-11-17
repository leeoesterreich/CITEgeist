#!/bin/bash
#SBATCH --job-name=v2_b22_moran_robustness
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/b22_moran_robustness_%j.log
# v2 B22 / B221 (run id from v42_run_ids.yaml: v2b22, rerun v2b221): robustness of the MDK-secretory bivariate Moran's I (checks 9, 10, 11, 13, 14, 15 in ONE job, about 40 min of compute; B221 adds check 17, about 3 min).
# Promoted from the exploratory checks; same math and seeds. Writes a NEW dir results/<protocol_sha12>/<RUN_ID> (refuses to overwrite).
# Required env (no defaults; set by .commands/submit-v2-moran.sh via chain_preflight.py): COMP_RUN (the B5 run dir), MORANS_RUN (B9), CONTROLS_RUN (B17), KERATIN_RUN (B12),
# E2_RUN (B21), OUT_ROOT, RUN_ID, V2_NAMESPACE_ROOT. MORAN_STAGE (B22 default, B221 for the rerun; set by the submit script) picks the stage; B221 is gated arm by arm against the registered B22 run.
# 16 threads, as the exploratory jobs the stage reproduces (the checks 9-15 jobs ran with OMP/OPENBLAS/MKL = 16).
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
: "${COMP_RUN:?COMP_RUN (the B5 run dir) is required}" "${MORANS_RUN:?MORANS_RUN (the B9 run dir) is required}" "${CONTROLS_RUN:?CONTROLS_RUN (the B17 run dir) is required}"
: "${KERATIN_RUN:?KERATIN_RUN (the B12 run dir) is required}" "${E2_RUN:?E2_RUN (the B21 run dir) is required}" "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B22 id) is required}"
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
[ ! -e "$OUT_ROOT" ] || { echo "refusing: $OUT_ROOT already exists (a run dir is never overwritten)" >&2; exit 1; }
"$PY" -u Benchmarking/midkine_remediation/moran_robustness.py \
  --compartments-run "$COMP_RUN" --morans-run "$MORANS_RUN" --controls-run "$CONTROLS_RUN" --keratin-run "$KERATIN_RUN" --e2-run "$E2_RUN" \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" --stage "${MORAN_STAGE:-B22}" --protocol Benchmarking/midkine_remediation/protocol.json
echo "done $RUN_ID end=$(date -Is)"
