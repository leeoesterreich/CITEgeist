#!/bin/bash
#SBATCH --job-name=spots_baseline
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
# SPOTS baseline on the Xenium RCC benchmark (6- and 7-type). Env: CODE_DIR (clean detached worktree for a clean
# provenance stamp), OUT_DIR (absolute results dir). Inputs are read relative to CODE_DIR's REPO via spots_baseline.py,
# so DATA_DIR must hold the untracked h5ad inputs: the script resolves them under CODE_DIR, hence the symlink made by
# the submit wrapper.
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
cd "$CODE_DIR"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$SLURM_JOB_ID NUMBA_CACHE_DIR=/tmp/numba_$SLURM_JOB_ID
PY="$(cg_env_python CITEgeist_env)"
echo "code=$CODE_DIR head=$(git rev-parse HEAD) out=$OUT_DIR host=$(hostname) start=$(date -Is)"
for N in 6 7; do "$PY" -u Benchmarking/xenium_benchmarking/spots_baseline/spots_baseline.py "$N" "$OUT_DIR/x.json"; done
echo "done end=$(date -Is)"
