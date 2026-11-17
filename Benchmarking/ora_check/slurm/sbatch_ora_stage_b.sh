#!/bin/bash
#SBATCH --job-name=v2_b19_ora
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --output=Benchmarking/ora_check/slurm/logs/b19_ora_%x_%A_%a.log
# v2 B19 (run id from v42_run_ids.yaml, currently v2b19): corrected ORA Test B on the rebuilt breast spots, plus the Xenium calibration. Post-hoc, descriptive.
# Promoted from run_ora_rebuilt.py (exploratory run results_rebuilt_20261103_1001); same math and seeds. Writes a NEW dir runs/<protocol_sha12>/<RUN_ID> (the submit wrapper refuses an existing one).
# STEP=task    array task 0..16 (0..11 breast, 12..16 Xenium regions); 4 CPU, 48 GB, 1.8-3.6 min each in the exploratory run (time limit 40 min)
# STEP=summary no array; checks the 17 task markers, writes ora_test_b_results.json and its stamp (run it with --dependency=afterok:<array job id>)
# Required env (no defaults; set by .commands/submit-v2-ora.sh via chain_preflight.py): MODULE3_DIR (B1 dir with p0/), NUCLEI_DIR (B2 dir), OUT_ROOT, RUN_ID, V2_NAMESPACE_ROOT.
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 PYTHONDONTWRITEBYTECODE=1
PY="$(cg_env_python CITEgeist_env)"
STEP=${STEP:?STEP (task or summary) is required}
: "${MODULE3_DIR:?MODULE3_DIR (the B1 dir) is required}" "${NUCLEI_DIR:?NUCLEI_DIR (the B2 dir) is required}" "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B19 id) is required}"
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
echo "step=$STEP task=${SLURM_ARRAY_TASK_ID:-none} host=$(hostname) code=$(git --no-optional-locks rev-parse HEAD) start=$(date -Is)"
ARGS=(--run-id "$RUN_ID" --out-root "$OUT_ROOT" --module3-dir "$MODULE3_DIR" --nuclei-dir "$NUCLEI_DIR" --protocol Benchmarking/midkine_remediation/protocol.json)
if [ "$STEP" = task ]; then
  "$PY" -u Benchmarking/ora_check/ora_stage_b.py task "${SLURM_ARRAY_TASK_ID:?array task id required}" "${ARGS[@]}"
else
  "$PY" -u Benchmarking/ora_check/ora_stage_b.py summary "${ARGS[@]}"
fi
echo "done $RUN_ID end=$(date -Is)"
