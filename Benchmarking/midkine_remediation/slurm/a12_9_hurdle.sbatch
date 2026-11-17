#!/bin/bash
#SBATCH --job-name=a12_9_hurdle
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_9_%j.log
# A12.9 real-data hurdle: H1 for all 12 specimens + A01-A07 sweep + H2/H3/H4 + negative control + stamped
# commot_compartment_hurdle.json (CITEgeist_env; the COMMOT-injection positive control runs separately in envs/COMMOT).
# Env: CODE_DIR (checkout to run from; default main tree), COMPARTMENTS_RUN, COMMOT_ROOT, COVERAGE_RUN (the stamped B24 run: measured synthetic coverage under rna="total"),
# OUT_ROOT (.../<protocol_sha12>/<RUN_ID>), RUN_ID; optional POSITIVE_CONTROL, SPECIMEN.
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=8
PY="$(cg_env_python CITEgeist_env)"
EXTRA=(); if [ -n "${POSITIVE_CONTROL:-}" ]; then EXTRA=(--positive-control "$POSITIVE_CONTROL"); fi
if [ -n "${SPECIMEN:-}" ]; then EXTRA+=(--specimen "$SPECIMEN"); fi
if [ "${FINALIZE:-0}" = 1 ]; then EXTRA+=(--finalize); fi
"$PY" -u Benchmarking/midkine_remediation/src/a12_9_real.py \
  --compartments-run "$COMPARTMENTS_RUN" --commot-root "$COMMOT_ROOT" --coverage-run "${COVERAGE_RUN:?COVERAGE_RUN (the stamped B24 run dir) is required}" --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --protocol Benchmarking/midkine_remediation/protocol.json "${EXTRA[@]}"
echo "done A12.9 end=$(date -Is)"
