#!/bin/bash
# Resubmit only the reference-based CARD jobs after MuSiC is installed
#
# Usage: ./resubmit_reference_jobs.sh [MUSIC_JOB_ID]

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
set -e

MUSIC_JOB_ID=${1:-""}
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CARD/slurm"

cd "${BASE_DIR}"

echo "=============================================="
echo "CARD Reference-Based Jobs Resubmission"
echo "=============================================="

if [ -n "$MUSIC_JOB_ID" ]; then
    DEP="--dependency=afterok:${MUSIC_JOB_ID}"
    echo "Waiting for MuSiC install job: ${MUSIC_JOB_ID}"
else
    DEP=""
fi

# Resubmit Xenium reference mode
echo ""
echo "Submitting Xenium reference mode..."
REF_JOB=$(sbatch ${DEP} run_all_regions.sh | awk '{print $4}')
echo "  Xenium reference job: ${REF_JOB}"

echo ""
echo "=============================================="
echo "Reference jobs resubmitted!"
echo "=============================================="
echo "Job: ${REF_JOB}"
