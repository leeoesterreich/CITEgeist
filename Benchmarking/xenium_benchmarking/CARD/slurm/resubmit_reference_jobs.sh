#!/bin/bash
# Resubmit only the reference-based CARD jobs after MuSiC is installed
#
# Usage: ./resubmit_reference_jobs.sh [MUSIC_JOB_ID]

set -e

MUSIC_JOB_ID=${1:-""}
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm"

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
