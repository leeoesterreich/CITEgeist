#!/usr/bin/env bash
#SBATCH --job-name=m1_recount
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --exclude=htc-n28,htc-n50,htc-n31
#SBATCH --array=0-13
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
# Module 1 recount, one specimen folder per array task. Output dir must be NEW (a rerun writes a new dir, never over module1_recount_v42).
# Usage: OUT=Benchmarking/dictionary_discovery/module1_recount_<tag> sbatch Benchmarking/dictionary_discovery/run_module1_recount.sbatch
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../_paths.sh"
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
cd ${CITEGEIST_REPO_ROOT}
: "${OUT:?set OUT to a new output directory}"
mkdir -p "$OUT"
export HDF5_USE_FILE_LOCKING=FALSE MPLCONFIGDIR=/tmp/mpl_$$ OMP_NUM_THREADS=4
"$(cg_env_python CITEgeist_env)" -u Benchmarking/dictionary_discovery/run_module1_recount.py "$OUT"
