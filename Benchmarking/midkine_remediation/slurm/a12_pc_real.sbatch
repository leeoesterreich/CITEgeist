#!/bin/bash
#SBATCH --job-name=a12_pc_real
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/pcreal_%j.log
# Positive-control COMMOT rerun on the REAL headline specimen (A12.9 step; not submitted in A12.4-A12.8). envs/COMMOT.
# Env: CODE_DIR (checkout to run from; default main tree), LAYERS_ROOT, COMP_RUN, OUT_ROOT (.../<protocol_sha12>/<RUN_ID>), RUN_ID.
# Required env MODULE3_DIR: the stamped B1 proportions dir (output/module3_rna_refined_7type_<ts>/<B1 id>); there is no default (the old default output/module3_7type is v40, and
# marker-share dirs are v41). Set by .commands via chain_preflight.py; V2_NAMESPACE_ROOT is the main checkout (stamped paths are relative to it).
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
M3=(--qp-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}")
"$(cg_env_python COMMOT)" -u Benchmarking/midkine_remediation/src/positive_control.py \
  --specimen HCC22-088-P4-S2_1i_rep --layers-root "$LAYERS_ROOT" --compartments-run "$COMP_RUN" --out-root "$OUT_ROOT" --run-id "$RUN_ID" \
  --protocol Benchmarking/midkine_remediation/protocol.json "${M3[@]}"
