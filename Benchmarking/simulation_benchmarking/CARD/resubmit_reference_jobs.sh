#!/bin/bash
# Resubmit only the reference-based CARD simulation jobs after MuSiC is installed
#
# Usage: ./resubmit_reference_jobs.sh [MUSIC_JOB_ID]

set -e

MUSIC_JOB_ID=${1:-""}
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CARD"

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
