#!/bin/bash
#SBATCH --job-name=spots_sens
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# SPOTS (reading a) under the CITEgeist protein-sensitivity perturbations. Env: CODE_DIR (clean detached worktree), OUT_DIR.
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
"$PY" -u Benchmarking/protein_sensitivity/spots_sensitivity.py run "$OUT_DIR"
"$PY" -u Benchmarking/protein_sensitivity/spots_sensitivity.py aggregate "$OUT_DIR"
echo "done end=$(date -Is)"
