#!/bin/bash
#SBATCH --job-name=CITEgeist_array
#SBATCH --output=benchmarking_logs/CITEgeist_array_%A_%a.log
#SBATCH --error=benchmarking_logs/CITEgeist_array_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --partition=HTC
#SBATCH --array=0-23 # 2 profile modes (manual, auto), 2 test sets, 2 alpha, 1 lambda_reg, 3 max_y_change = 24 combinations

# Activate conda environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Load Gurobi module
module load gurobi/12.0.3
echo "Loaded gurobi module"

# Change to working directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Define parameter arrays
lambda_reg=(1)
alpha_elastic=(0.7 0.9)
max_y_change=(0.2 0.4 0.8)
TEST_SETS=("mixed" "high_seg")
PROFILE_MODES=("manual" "auto")  # manual or auto-discovery (spatial colocalization pipeline)

# Auto-discovery parameters (only used when profile_mode="auto")
# Uses the new spatial colocalization pipeline (Modules 1 + 2a + 2b + 2c)
TOP_K=3                    # Mutual top-k for profile pairing (default: 3, tested value)
DISCOVERY_SEED=1234

# Spatial colocalization pipeline parameters
MORANS_K=8                 # Number of neighbors for spatial statistics (Moran's I, bivariate Moran's I)
SMOOTH_K=6                 # Number of neighbors for spatial smoothing before Moran's I calculations
                           # Smoothing reduces noise and improves detection in mixed tissue contexts
N_PERMUTATIONS=999         # Permutations for significance testing (999 for stable p-values)
FDR_THRESHOLD=0.05         # FDR threshold for significant markers/pairs (bivariate Moran's I p-value)
GMM_SNR_THRESHOLD=1.0      # Minimum GMM SNR for marker to be interesting
# NOTE: Using adaptive thresholds (learned from data) for kurtosis and Moran's I
# The adaptive approach works correctly on RAW data (before CLR transformation)
# IMPORTANT: Bivariate Moran's I is now the default p-value source for FDR correction (more robust for mixed data)

# Module 2c: Spatial variance-based profile selection
VARIANCE_TARGET=0.90       # Target fraction of spatial variance to explain
MIN_MARGINAL_GAIN=0.005    # Minimum marginal variance gain (0.5%) - lower = more profiles
MAX_PROFILES=15            # Maximum number of profiles to select

# Laplacian smoothing parameters (for spatial coherence in proportion optimization)
LAMBDA_LAPLACIAN=0.1       # Laplacian smoothing weight (0 to disable)
LAPLACIAN_K=8              # Number of neighbors for Laplacian graph

# Calculate indices based on array task ID
# Order: profile_mode -> test_set -> alpha_elastic -> max_y_change
# Total: 2 * 2 * 2 * 3 = 24 combinations
n_max_y=3
n_alpha=2
n_test_sets=2

profile_mode_index=$((SLURM_ARRAY_TASK_ID / (n_test_sets * n_alpha * n_max_y)))
test_set_index=$(((SLURM_ARRAY_TASK_ID / (n_alpha * n_max_y)) % n_test_sets))
alpha_elastic_index=$(((SLURM_ARRAY_TASK_ID / n_max_y) % n_alpha))
max_y_change_index=$((SLURM_ARRAY_TASK_ID % n_max_y))

# Get parameter values
PROFILE_MODE=${PROFILE_MODES[$profile_mode_index]}
TEST_SET=${TEST_SETS[$test_set_index]}
lambda_reg=${lambda_reg[0]}
alpha_elastic=${alpha_elastic[$alpha_elastic_index]}
max_y_change=${max_y_change[$max_y_change_index]}

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/${TEST_SET}/${PROFILE_MODE}_profiles/"

mkdir -p "$OUTPUT_FOLDER"
mkdir -p benchmarking_logs

echo "Running with parameters:"
echo "  - profile_mode=$PROFILE_MODE"
echo "  - lambda_reg=$lambda_reg"
echo "  - alpha_elastic=$alpha_elastic"
echo "  - max_y_change=$max_y_change"
echo "  - test_set=$TEST_SET"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Build command with conditional auto-profiles flag
CMD="python tests/test_citegeist_simulated.py \
    --radius 4 \
    --lambda_reg $lambda_reg \
    --alpha_elastic $alpha_elastic \
    --max_y_change $max_y_change \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER"

# Add auto-profiles flags if in auto mode (spatial colocalization pipeline)
if [ "$PROFILE_MODE" == "auto" ]; then
    CMD="$CMD --auto-profiles"
    CMD="$CMD --top-k $TOP_K"
    CMD="$CMD --discovery-seed $DISCOVERY_SEED"
    # Spatial colocalization pipeline parameters (Module 1 + 2a)
    CMD="$CMD --morans-k $MORANS_K"
    CMD="$CMD --smooth-k $SMOOTH_K"
    CMD="$CMD --n-permutations $N_PERMUTATIONS"
    CMD="$CMD --fdr-threshold $FDR_THRESHOLD"
    CMD="$CMD --gmm-snr-threshold $GMM_SNR_THRESHOLD"
    # Module 2c: Spatial variance-based profile selection
    CMD="$CMD --variance-target $VARIANCE_TARGET"
    CMD="$CMD --min-marginal-gain $MIN_MARGINAL_GAIN"
    # Laplacian smoothing parameters
    CMD="$CMD --lambda-laplacian $LAMBDA_LAPLACIAN"
    CMD="$CMD --laplacian-k $LAPLACIAN_K"
    echo "  - top_k=$TOP_K"
    echo "  - discovery_seed=$DISCOVERY_SEED"
    echo "  - morans_k=$MORANS_K"
    echo "  - smooth_k=$SMOOTH_K"
    echo "  - n_permutations=$N_PERMUTATIONS"
    echo "  - fdr_threshold=$FDR_THRESHOLD"
    echo "  - gmm_snr_threshold=$GMM_SNR_THRESHOLD"
    echo "  - variance_target=$VARIANCE_TARGET"
    echo "  - min_marginal_gain=$MIN_MARGINAL_GAIN"
    echo "  - lambda_laplacian=$LAMBDA_LAPLACIAN"
    echo "  - laplacian_k=$LAPLACIAN_K"
fi

echo "Executing: $CMD"
eval $CMD

echo "Job completed at $(date)"
