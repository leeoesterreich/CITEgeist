#!/bin/bash
#SBATCH --job-name=v2loro
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --exclude=htc-1024-n2,htc-n37
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# v42 stage X11: the leave-one-region-out tuning audit of frozen v2 (v2_loro.py). Env: SHA (the commit that holds this code).
# Submitted by .commands/submit-v2-loro.sh. CPU only (72 settings x 5 regions x 2 type sets of closed-form EM fits). Runs from the main checkout after
# proving every tracked file of the code tree equals commit SHA; v2_loro.py proves it again, computes, gates against the stored r8 row / X2 / X1
# (a miss writes gate_failure.json and exits non-zero with nothing stamped) and stamps loro_audit.json last.
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
: "${SHA:?set SHA to the commit that holds this code}"
PY="$(cg_env_python CITEgeist_env)"
REPO=${CITEGEIST_REPO_ROOT}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 HDF5_USE_FILE_LOCKING=FALSE
export MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID PYTHONPATH=$REPO:$REPO/Benchmarking
cd "$REPO"
RF=Benchmarking/xenium_benchmarking/rna_refinement
"$PY" -u $RF/v2_official.py check-code "$SHA"
echo "stage=X11 sha=$SHA host=$(hostname) start=$(date -Is)"
"$PY" -u $RF/v2_loro.py run --stamp
echo "stage=X11 done $(date -Is)"
