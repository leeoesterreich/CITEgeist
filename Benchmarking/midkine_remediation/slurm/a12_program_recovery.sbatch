#!/bin/bash
#SBATCH --job-name=a12_prog_recovery
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=08:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_prog_recovery_%j.log
# A12 program recovery (Figure 5D): legacy recovery rule on Cell2Location / Tangram programs (recomputed from their
# unchanged GEX layers) and on the A12.11 CITEgeist Module 4 programs, 12 specimens. CPU only (sklearn NMF + Module 4).
# Env: RUN_ID, OUT_ROOT (absolute, .../results/<protocol_sha12>/<RUN_ID>), CODE_DIR (checkout to run the code from;
# default the main tree), DATA_DIR (tree holding output/; default the main tree).
# Required (no defaults; the v40 defaults were output/module3_7type, r1, n1): MODULE3_DIR (stamped B1 dir), LAYERS_RUN_ID (B4 id), NMF_RUN (the B10 run dir).
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
CODE_DIR=${CODE_DIR:-${CITEGEIST_REPO_ROOT}}
DATA_DIR=${DATA_DIR:-${CITEGEIST_REPO_ROOT}}
cd "$CODE_DIR"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=8
PY="$(cg_env_python CITEgeist_env)"
echo "run_id=$RUN_ID out=$OUT_ROOT code=$CODE_DIR head=$(git rev-parse HEAD) host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/a12_program_recovery.py \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --layers-root "$DATA_DIR/output/sace_layers_module3_7type" --layers-run-id "${LAYERS_RUN_ID:?LAYERS_RUN_ID (the B4 id) is required}" \
  --module3-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}" \
  --c2l-root "$DATA_DIR/output/competitor_gex/cell2location" \
  --tangram-root "$DATA_DIR/output/competitor_gex/tangram" \
  --legacy-programs-root "$DATA_DIR/output/competitor_gex" \
  --a12-11-counts "${NMF_RUN:?NMF_RUN (the B10 run dir) is required}/nmf_counts.json" \
  --protocol Benchmarking/midkine_remediation/protocol.json
echo "done end=$(date -Is)"
