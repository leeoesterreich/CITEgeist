#!/usr/bin/env bash
#SBATCH --job-name=comp_refs_full
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=192G
#SBATCH --time=12:00:00
# Build the full_size (all 54,850 reference cells) export through the arm pipeline, then compare it to processed/
# (semantic, full-content equivalence). The gate runs in the same job; a content difference fails the job (STOP).
_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
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
export MPLCONFIGDIR=/tmp HDF5_USE_FILE_LOCKING=FALSE
export TMPDIR="${SLURM_TMPDIR:-/tmp}"   # node-local scratch, not the inherited /vast TMPDIR
REPO=${CITEGEIST_REPO_ROOT}
PY="$(cg_env_python CITEgeist_env)"
cd "$REPO/Benchmarking/xenium_benchmarking"
$PY -u -m gt_rebuild.build_composition_reference --endpoint full_size
$PY -u -m gt_rebuild.build_composition_reference --endpoint full_size --verify-only
