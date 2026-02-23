#!/bin/bash
# Submit all CARD simulation benchmark jobs with dependencies
#
# Usage: ./submit_all_jobs.sh [ENV_JOB_ID]

set -e

ENV_JOB_ID=${1:-""}
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CARD"

echo "=============================================="
echo "CARD Simulation Benchmark Job Submission"
echo "=============================================="

# Build dependency string
if [ -n "$ENV_JOB_ID" ]; then
    DEP_ENV="--dependency=afterok:${ENV_JOB_ID}"
    echo "Waiting for environment setup job: ${ENV_JOB_ID}"
else
    DEP_ENV=""
    echo "No environment dependency specified"
fi

# Step 0: Convert reference to CSV (for simulation benchmarks)
echo ""
echo "Submitting reference CSV conversion..."
REF_CONV_JOB=$(sbatch --job-name=card_ref_csv \
    --output="${BASE_DIR}/ref_convert_%j.out" \
    --error="${BASE_DIR}/ref_convert_%j.err" \
    --time=00:30:00 --mem=64G --cpus-per-task=2 \
    --cluster=htc --partition=htc \
    --mail-type=FAIL --mail-user=alc376@pitt.edu \
    --wrap="source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env && cd ${BASE_DIR} && python convert_reference_to_csv.py" \
    | awk '{print $4}')
echo "  Reference conversion job: ${REF_CONV_JOB}"

# Step 1: Generate marker genes for simulation (for reference-free mode)
echo ""
echo "Submitting marker genes generation for simulation..."
MARKER_JOB=$(sbatch ${DEP_ENV} --job-name=card_sim_markers \
    --output="${BASE_DIR}/marker_gen_%j.out" \
    --error="${BASE_DIR}/marker_gen_%j.err" \
    --time=01:00:00 --mem=64G --cpus-per-task=4 \
    --cluster=htc --partition=htc \
    --mail-type=FAIL --mail-user=alc376@pitt.edu \
    --wrap="source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env && cd ${BASE_DIR} && python generate_marker_genes_simulation.py" \
    | awk '{print $4}')
echo "  Marker genes job: ${MARKER_JOB}"

# Step 2: high_seg reference mode (needs env + reference CSV)
echo ""
echo "Submitting high_seg reference mode..."
cd "${BASE_DIR}/high_seg"
mkdir -p logs
# Depends on both env setup and reference conversion
if [ -n "$DEP_ENV" ]; then
    HS_REF_JOB=$(sbatch ${DEP_ENV},afterok:${REF_CONV_JOB} CARD_pipeline_high_seg.sh | awk '{print $4}')
else
    HS_REF_JOB=$(sbatch --dependency=afterok:${REF_CONV_JOB} CARD_pipeline_high_seg.sh | awk '{print $4}')
fi
echo "  high_seg reference job: ${HS_REF_JOB}"

# Step 3: high_seg reference-free mode (needs markers)
echo ""
echo "Submitting high_seg reference-free mode..."
HS_REFFREE_JOB=$(sbatch --dependency=afterok:${MARKER_JOB} CARD_pipeline_high_seg_reffree.sh | awk '{print $4}')
echo "  high_seg ref-free job: ${HS_REFFREE_JOB}"

# Step 4: mixed reference mode (needs env + reference CSV)
echo ""
echo "Submitting mixed reference mode..."
cd "${BASE_DIR}/mixed"
mkdir -p logs
if [ -n "$DEP_ENV" ]; then
    MX_REF_JOB=$(sbatch ${DEP_ENV},afterok:${REF_CONV_JOB} CARD_pipeline_mixed.sh | awk '{print $4}')
else
    MX_REF_JOB=$(sbatch --dependency=afterok:${REF_CONV_JOB} CARD_pipeline_mixed.sh | awk '{print $4}')
fi
echo "  mixed reference job: ${MX_REF_JOB}"

# Step 5: mixed reference-free mode (needs markers)
echo ""
echo "Submitting mixed reference-free mode..."
MX_REFFREE_JOB=$(sbatch --dependency=afterok:${MARKER_JOB} CARD_pipeline_mixed_reffree.sh | awk '{print $4}')
echo "  mixed ref-free job: ${MX_REFFREE_JOB}"

echo ""
echo "=============================================="
echo "All simulation jobs submitted!"
echo "=============================================="
echo ""
echo "Job summary:"
echo "  Environment setup:     ${ENV_JOB_ID:-'(not specified)'}"
echo "  Reference CSV:         ${REF_CONV_JOB}"
echo "  Marker genes:          ${MARKER_JOB}"
echo "  high_seg reference:    ${HS_REF_JOB}"
echo "  high_seg ref-free:     ${HS_REFFREE_JOB}"
echo "  mixed reference:       ${MX_REF_JOB}"
echo "  mixed ref-free:        ${MX_REFFREE_JOB}"
echo ""
echo "Monitor with: squeue -u \$USER"
