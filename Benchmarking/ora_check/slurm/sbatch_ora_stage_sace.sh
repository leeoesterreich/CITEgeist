#!/bin/bash
#SBATCH --job-name=v2_b20_ora_sace
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:30:00
#SBATCH --output=Benchmarking/ora_check/slurm/logs/b20_ora_sace_%x_%A_%a.log
# v2 B20 (run id from v42_run_ids.yaml, currently v2b20): ORA on the SACE layers of the rebuilt breast spots (Part 1 protocol, Part 2 Test A) plus a Xenium calibration.
# Post-hoc, descriptive. Promoted from run_ora_sace_rebuilt.py (exploratory run results_sace_rebuilt_20261103_1115); same math and seeds. Writes a NEW dir runs/<protocol_sha12>/<RUN_ID>
# (the submit wrapper refuses an existing one). The run dir holds the SACE layers (about 0.5 GB).
# STEP=task    array task 0..16 (0..11 breast, 12..16 Xenium regions); 8 CPU, 96 GB; breast 5-23 min, Xenium 2-8 min with the regenerated SACE arms (time limit 150 min; with 17 tasks running at once the slowest breast task took 49 min, the exploratory one 12)
# STEP=summary no array; checks the 17 task markers, writes ora_sace_results.json and its stamp (run it with --dependency=afterok:<array job id>)
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
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 PYTHONDONTWRITEBYTECODE=1
PY="$(cg_env_python CITEgeist_env)"
STEP=${STEP:?STEP (task or summary) is required}
: "${MODULE3_DIR:?MODULE3_DIR (the B1 dir) is required}" "${NUCLEI_DIR:?NUCLEI_DIR (the B2 dir) is required}" "${OUT_ROOT:?OUT_ROOT is required}" "${RUN_ID:?RUN_ID (the B20 id) is required}"
[ ! -x "$MAIN/.commands/_refuse_protected_outdir.sh" ] || "$MAIN/.commands/_refuse_protected_outdir.sh" "$OUT_ROOT"
echo "step=$STEP task=${SLURM_ARRAY_TASK_ID:-none} host=$(hostname) code=$(git --no-optional-locks rev-parse HEAD) start=$(date -Is)"
ARGS=(--run-id "$RUN_ID" --out-root "$OUT_ROOT" --module3-dir "$MODULE3_DIR" --nuclei-dir "$NUCLEI_DIR" --protocol Benchmarking/midkine_remediation/protocol.json)
if [ "$STEP" = task ]; then
  "$PY" -u Benchmarking/ora_check/ora_stage_sace.py task "${SLURM_ARRAY_TASK_ID:?array task id required}" "${ARGS[@]}"
else
  "$PY" -u Benchmarking/ora_check/ora_stage_sace.py summary "${ARGS[@]}"
fi
echo "done $RUN_ID end=$(date -Is)"
