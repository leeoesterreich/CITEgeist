#!/bin/bash
#SBATCH --job-name=a12_8_commot
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --array=0-11
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --output=Benchmarking/midkine_remediation/slurm/logs/a12_8_%A_%a.log
# A12.8 COMMOT run matrix, one specimen per array task, envs/COMMOT (one `import commot`, 100-220 s on NFS).
# Env: LAYERS_ROOT, COMP_RUN (A12.7 dir), OUT_ROOT (.../<protocol_sha12>/<RUN_ID>), RUN_ID, DATA_ROOT, optional VERIFY_SUBSET=1.
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
MAIN=${CITEGEIST_REPO_ROOT}   # data/outputs stay in the main tree; CODE_DIR (if set) is a clean detached worktree
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID OMP_NUM_THREADS=8
PY="$(cg_env_python COMMOT)"
SPECS=(HCC22-088-P1-S1 HCC22-088-P1-S2 HCC22-088-P2-S1 HCC22-088-P2-S2 HCC22-088-P3-S1_A HCC22-088-P3-S2 HCC22-088-P4-S1 HCC22-088-P4-S2_1i_rep HCC22-088-P5-S1 HCC22-088-P5-S2_F_rep HCC22-088-P6-S1 HCC22-088-P6-S2_D)
SPECIMEN=${SPECS[$SLURM_ARRAY_TASK_ID]}
M3=(--qp-dir "${MODULE3_DIR:?MODULE3_DIR (the B1 output dir) is required}")
EXTRA=(); if [ "${VERIFY_SUBSET:-0}" = 1 ]; then EXTRA=(--verify-subset); fi
if [ "${RESUME:-0}" = 1 ]; then EXTRA+=(--resume); fi
echo "specimen=$SPECIMEN out=$OUT_ROOT/$SPECIMEN host=$(hostname) start=$(date -Is)"
"$PY" -u Benchmarking/midkine_remediation/src/commot_runs.py --specimen "$SPECIMEN" --layers-root "$LAYERS_ROOT" --compartments-run "$COMP_RUN" \
  --out-root "$OUT_ROOT" --run-id "$RUN_ID" --data-root "$DATA_ROOT" --protocol Benchmarking/midkine_remediation/protocol.json "${M3[@]}" "${EXTRA[@]}"
echo "done $SPECIMEN end=$(date -Is)"
