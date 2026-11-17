#!/bin/bash
#SBATCH --job-name=m35_val
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --output=slurm_log/m35_val_%j.out
#SBATCH --error=slurm_log/m35_val_%j.err

# Usage: REGION=0 gpu-race-submit.sh sbatch_validate_m35.sh
# Requires GPU because CitegeistModel imports cuOPT at module level.
# GPU_RACE_READONLY_OK=1

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

REGION="${REGION:?Set REGION=0..4 before submitting}"

REPO_ROOT="${CITEGEIST_REPO_ROOT}"
SCRIPT_DIR="$REPO_ROOT/Benchmarking/xenium_benchmarking/CITEgeist/src"
PYTHON="$(cg_env_python CITEgeist_env)"

cd "$REPO_ROOT"

echo "=== M3.5 Validation: Region ${REGION} ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"

$PYTHON -u "$SCRIPT_DIR/run_module3_5_validation.py" --region "${REGION}"

echo "Completed: $(date)"

# After all regions complete, run aggregation manually:
# $PYTHON -u "$SCRIPT_DIR/run_module3_5_validation.py" --aggregate
