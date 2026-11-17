#!/bin/bash
#SBATCH --job-name=a12_10_morans
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_10_%j.log
# A12.10 real-data bivariate Moran's I: Definition C (12 specimens), Definition A on v4, the three-method
# comparator shared object, the superseded v3 recompute, aggregate/BH/G8/patient-sign-flip/trajectory, and
# one stamped namespaced morans_results.json. CITEgeist_env (esda/libpysal needed for the v3 recompute).
# Env: COMPARTMENTS_RUN, C2L_ROOT, TANGRAM_ROOT, OUT_ROOT (.../<protocol_sha12>/<RUN_ID>), RUN_ID, V4_ROOT (REQUIRED: the B3 run dir <...>/morphology_assignment_v4_hungarian/<sha12>/<B3 id>);
# V3_ROOT must stay unset in v42 (output/morphology_assignment_v3 is a v4 stand-in, so there is no v3 recompute).
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
MAIN=${CITEGEIST_REPO_ROOT}   # data/outputs stay in the main tree; CODE_DIR (if set) is a clean detached worktree
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=16
PY="$(cg_env_python CITEgeist_env)"
EXTRA=(); [ -n "${V3_ROOT:-}" ] && EXTRA+=(--v3-root "$V3_ROOT")
EXTRA+=(--v4-root "${V4_ROOT:?V4_ROOT (the B3 run dir) is required}")
"$PY" -u Benchmarking/midkine_remediation/src/a12_10_real.py \
  --compartments-run "$COMPARTMENTS_RUN" --c2l-root "$C2L_ROOT" --tangram-root "$TANGRAM_ROOT" \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --protocol Benchmarking/midkine_remediation/protocol.json "${EXTRA[@]}"
echo "done A12.10 end=$(date -Is)"
