#!/bin/bash
#SBATCH --job-name=xen_native_split
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --time=03:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
# Xenium: Tangram native layers + '<Method> proportional split' control layers, then re-run the GEX evaluator.
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
PY="$(cg_env_python CITEgeist_env)"
cd ${CITEGEIST_REPO_ROOT}
$PY -u Benchmarking/xenium_benchmarking/evaluation/src/build_native_and_split_layers.py
echo BUILD_DONE
