#!/bin/bash
#SBATCH --job-name=v2_b24_prod_coverage
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/b24_prod_coverage_%j.log
# v2 B24 (run id from v42_run_ids.yaml, currently v2b241): coverage of the nominal 95% GMR interval through the production descriptive.gmr_interval with rna="total"
# (protocol deviation 2026-11-06). Same four cases, 500 datasets each (seeds 1000..1499), 2000 resamples; the four cases run one after the other (about 10 minutes with 8 processes
# in the per-cell series: 60 s, 60 s, 250 s, 204 s). Writes a NEW dir results/<protocol_sha12>/<RUN_ID> (refuses to overwrite) with the stamped prod_coverage_total.json, a stamped input of B8.
# Required env (no defaults; set by .commands/submit-v2-coverage.sh via chain_preflight.py): OUT_ROOT, RUN_ID, V2_NAMESPACE_ROOT.
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY="$(cg_env_python CITEgeist_env)"
: "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B24 id) is required}"
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
[ ! -e "$OUT_ROOT" ] || { echo "refusing: $OUT_ROOT already exists (a run dir is never overwritten)" >&2; exit 1; }
"$PY" -u Benchmarking/midkine_remediation/src/prod_coverage.py --out-root "$OUT_ROOT" --run-id "$RUN_ID" --procs "${SLURM_CPUS_PER_TASK}" \
  --protocol Benchmarking/midkine_remediation/protocol.json
echo "done $RUN_ID end=$(date -Is)"
