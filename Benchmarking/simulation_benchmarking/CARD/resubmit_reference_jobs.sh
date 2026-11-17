#!/bin/bash
# Resubmit only the reference-based CARD simulation jobs after MuSiC is installed
#
# Usage: ./resubmit_reference_jobs.sh [MUSIC_JOB_ID]

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
set -e

MUSIC_JOB_ID=${1:-""}
BASE_DIR="${CITEGEIST_REPO_ROOT}/Benchmarking/simulation_benchmarking/CARD"

echo "=============================================="
echo "CARD Simulation Reference Jobs Resubmission"
echo "=============================================="

if [ -n "$MUSIC_JOB_ID" ]; then
    DEP="--dependency=afterok:${MUSIC_JOB_ID}"
    echo "Waiting for MuSiC install job: ${MUSIC_JOB_ID}"
else
    DEP=""
fi

# high_seg reference
echo ""
echo "Submitting high_seg reference mode..."
cd "${BASE_DIR}/high_seg"
HS_JOB=$(sbatch ${DEP} CARD_pipeline_high_seg.sh | awk '{print $4}')
echo "  high_seg reference job: ${HS_JOB}"

# mixed reference
echo ""
echo "Submitting mixed reference mode..."
cd "${BASE_DIR}/mixed"
MX_JOB=$(sbatch ${DEP} CARD_pipeline_mixed.sh | awk '{print $4}')
echo "  mixed reference job: ${MX_JOB}"

echo ""
echo "=============================================="
echo "Reference jobs resubmitted!"
echo "=============================================="
echo "  high_seg: ${HS_JOB}"
echo "  mixed:    ${MX_JOB}"
