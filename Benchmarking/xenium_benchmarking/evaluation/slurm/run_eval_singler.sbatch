#!/bin/bash
#SBATCH --job-name=eval_singler
#SBATCH --output=eval_singler_%j.log
#SBATCH --error=eval_singler_%j.log
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc
#SBATCH --partition=htc

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

module purge

PYTHON="$(cg_env_python CITEgeist_env)"
SRC=${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/evaluation/src

echo "=== Evaluating 6-type (SingleR GT) ==="
$PYTHON -u "${SRC}/compare_all_methods.py" --gt-type singler_6type

echo ""
echo "=== Evaluating 7-type (SingleR GT) ==="
$PYTHON -u "${SRC}/compare_all_methods.py" --gt-type singler_7type

echo ""
echo "=== Done ==="
