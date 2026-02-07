#!/bin/bash
#SBATCH --job-name=permarker_beta_compare
#SBATCH --output=benchmarking_logs/permarker_beta_compare_%A_%a.log
#SBATCH --error=benchmarking_logs/permarker_beta_compare_%A_%a.log
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --partition=HTC
#SBATCH --array=0-1  # 0 = per-marker beta (new), 1 = per-celltype beta (legacy)

# Comparison test: Per-marker beta vs Per-celltype beta
# This script runs both approaches on simulated data to compare performance

# Activate conda environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Load Gurobi module
module load gurobi/12.0.3
echo "Loaded gurobi module"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Fixed parameters for comparison (use best settings from previous benchmarks)
RADIUS=4
LAMBDA_REG=1
ALPHA_ELASTIC=0.7
MAX_Y_CHANGE=0.4
TEST_SET="high_seg"
PROFILE_MODE="manual"  # Use manual profiles for controlled comparison

# Per-marker beta parameters
BETA_MIN=0.1
BETA_MAX=2.0

# Determine which beta mode to use based on array task ID
if [ "$SLURM_ARRAY_TASK_ID" == "0" ]; then
    BETA_MODE="per_marker"
    PER_MARKER_FLAG="--per-marker-beta"
    echo "=== RUNNING WITH PER-MARKER BETA (NEW APPROACH) ==="
else
    BETA_MODE="per_celltype"
    PER_MARKER_FLAG="--no-per-marker-beta"
    echo "=== RUNNING WITH PER-CELLTYPE BETA (LEGACY APPROACH) ==="
fi

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/permarker_beta_comparison/${BETA_MODE}/"

mkdir -p "$OUTPUT_FOLDER"
mkdir -p benchmarking_logs

echo "Running comparison test with parameters:"
echo "  - beta_mode=$BETA_MODE"
echo "  - radius=$RADIUS"
echo "  - lambda_reg=$LAMBDA_REG"
echo "  - alpha_elastic=$ALPHA_ELASTIC"
echo "  - max_y_change=$MAX_Y_CHANGE"
echo "  - test_set=$TEST_SET"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Build and execute command
CMD="python tests/test_citegeist_simulated.py \
    --radius $RADIUS \
    --lambda_reg $LAMBDA_REG \
    --alpha_elastic $ALPHA_ELASTIC \
    --max_y_change $MAX_Y_CHANGE \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER \
    --beta-min $BETA_MIN \
    --beta-max $BETA_MAX \
    $PER_MARKER_FLAG"

echo "Executing: $CMD"
eval $CMD

echo "Job completed at $(date)"
echo "Results saved to: $OUTPUT_FOLDER"
