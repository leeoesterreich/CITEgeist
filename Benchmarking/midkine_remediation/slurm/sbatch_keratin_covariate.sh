#!/bin/bash
#SBATCH --job-name=v2_b23_keratin_cov
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/b23_keratin_cov_%j.log
# v2 B23 (run id from v42_run_ids.yaml, currently v2b23): which single covariate absorbs the basal-keratin vs secretory association (check 1). Descriptive only.
# Promoted from the exploratory check 1; same math and seeds as B12 (v2b121), which it is gated against (every specimen, fail closed). Writes a NEW dir results/<protocol_sha12>/<RUN_ID> (refuses to overwrite).
# Required env (no defaults; set by .commands/submit-v2-keratin-cov.sh via chain_preflight.py): COMP_RUN (the B5 run dir), KERATIN_RUN (the B12 run dir), OUT_ROOT, RUN_ID, V2_NAMESPACE_ROOT.
# EXPLORATORY_RESULTS (optional, defaults to the check 1 run of 2026-11-02 in the main tree): compared with the stage, never gated.
# 16 threads, as the B12 run (canonical_eigh makes the surrogates independent of the thread count; the pin is kept anyway).
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
: "${COMP_RUN:?COMP_RUN (the B5 run dir) is required}" "${KERATIN_RUN:?KERATIN_RUN (the B12 run dir) is required}" "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B23 id) is required}"
EXPLORATORY_RESULTS=${EXPLORATORY_RESULTS:-$MAIN/Benchmarking/midkine_remediation/results_exploratory_v42/check1_covariate_absorption_20261102T231500/check1_results.json}
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
[ ! -e "$OUT_ROOT" ] || { echo "refusing: $OUT_ROOT already exists (a run dir is never overwritten)" >&2; exit 1; }
"$PY" -u Benchmarking/midkine_remediation/keratin_covariate.py \
  --compartments-run "$COMP_RUN" --keratin-run "$KERATIN_RUN" --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --protocol Benchmarking/midkine_remediation/protocol.json --exploratory-results "$EXPLORATORY_RESULTS"
echo "done $RUN_ID end=$(date -Is)"
