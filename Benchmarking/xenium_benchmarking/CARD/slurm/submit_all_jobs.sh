#!/bin/bash
# Submit all CARD benchmark jobs with dependencies
#
# Usage: ./submit_all_jobs.sh [ENV_JOB_ID]
#   If ENV_JOB_ID is provided, jobs wait for environment setup to complete

set -e

ENV_JOB_ID=${1:-""}
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CARD/slurm"

cd "${BASE_DIR}"

echo "=============================================="
echo "CARD Xenium Benchmark Job Submission"
echo "=============================================="

# Build dependency string
if [ -n "$ENV_JOB_ID" ]; then
    DEP_ENV="--dependency=afterok:${ENV_JOB_ID}"
    echo "Waiting for environment setup job: ${ENV_JOB_ID}"
else
    DEP_ENV=""
    echo "No environment dependency specified"
fi

# Step 1: Generate marker genes (for reference-free mode)
echo ""
echo "Submitting marker genes generation..."
if [ -n "$DEP_ENV" ]; then
    MARKER_JOB=$(sbatch ${DEP_ENV} setup_markers.sh | awk '{print $4}')
else
    MARKER_JOB=$(sbatch setup_markers.sh | awk '{print $4}')
fi
echo "  Marker genes job: ${MARKER_JOB}"

# Step 2: Reference mode (can run after env setup, doesn't need markers)
echo ""
echo "Submitting reference mode benchmark..."
if [ -n "$DEP_ENV" ]; then
    REF_JOB=$(sbatch ${DEP_ENV} run_all_regions.sh | awk '{print $4}')
else
    REF_JOB=$(sbatch run_all_regions.sh | awk '{print $4}')
fi
echo "  Reference mode job: ${REF_JOB}"

# Step 3: Reference-free mode (needs markers)
echo ""
echo "Submitting reference-free mode benchmark..."
REFFREE_JOB=$(sbatch --dependency=afterok:${MARKER_JOB} run_all_regions_reffree.sh | awk '{print $4}')
echo "  Reference-free mode job: ${REFFREE_JOB}"

echo ""
echo "=============================================="
echo "All jobs submitted!"
echo "=============================================="
echo ""
echo "Job summary:"
echo "  Environment setup: ${ENV_JOB_ID:-'(not specified)'}"
echo "  Marker genes:      ${MARKER_JOB}"
echo "  Reference mode:    ${REF_JOB}"
echo "  Reference-free:    ${REFFREE_JOB}"
echo ""
echo "Monitor with: squeue -u \$USER"
