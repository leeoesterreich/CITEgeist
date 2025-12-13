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
#SBATCH --array=0-23 # 2 profile modes, 2 test sets, 2 alpha, 1 lambda_reg, 3 max_y_change = 24 combinations

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
PROFILE_MODES=("manual" "auto")  # New: manual vs auto-discovery

# Auto-discovery parameters (only used when profile_mode="auto")
MAX_PROFILE_SIZE=2
DISCOVERY_SEED=1234

# Robust auto-discovery parameters (new defaults for mixed/dysplastic environments)
MORANS_K="3,5,8,12"        # Multi-scale spatial filtering
MORANS_I_THRESHOLD="0.1"   # Stricter threshold to filter nonspecific markers (was -0.1)
MODEL_SELECTION="cv"        # CV-based model selection (prevents premature stopping)
CV_FOLDS=5
MIN_PROFILES=2
# Note: hierarchical mode disabled by default, enable for tissues with shared markers (e.g., T-cell subtypes)

# NEW: Reconstruction-based discovery parameters
SELECTION_METHOD="miqp_hierarchical"  # Options: reconstruction, permutation, miqp, miqp_hierarchical
MIN_RECON_IMPROVEMENT=0.05         # Minimum reconstruction improvement to include profile
ABUNDANCE_ADAPTIVE="--abundance-adaptive"  # Enable abundance-adaptive marker classification
UBIQUITOUS_CV_THRESHOLD=0.5        # CV threshold for ubiquitous classification
RARE_PRESENCE_THRESHOLD=0.15       # Presence threshold for rare classification (high_frac < this)
REDUNDANCY_THRESHOLD=0.7           # Correlation threshold for redundancy check (LOWERED from 0.9)

# NEW: MIQP-specific parameters (only used when SELECTION_METHOD="miqp" or "miqp_hierarchical")
MIQP_LAMBDA_SPATIAL=0.1            # Spatial penalty weight (balanced priority)
MIQP_LAMBDA_COMPLEXITY=0.15        # Profile count penalty (INCREASED from 0.05 - strongly penalize many profiles)
MIQP_TIME_LIMIT=300.0              # Solver time limit in seconds
MIQP_GAP=0.01                      # Acceptable MIP optimality gap (1%)

# NEW: Hierarchical MIQP parameters (only used when SELECTION_METHOD="miqp_hierarchical")
MIQP_LAMBDA_OVERLAP=0.8            # Same-level competition penalty (INCREASED from 0.7)
MIQP_LAMBDA_ORPHAN=0.3             # Orphan penalty (child without parent) (INCREASED)
MIQP_LAMBDA_SPARSITY=0.6           # Spot-level sparsity penalty (INCREASED from 0.5)
MIQP_SPARSITY_AGGREGATION="mean"   # How to compute spot-profile fit: mean, min, geometric

# NEW: Spatial Mixture Model (SMM) background correction parameters
USE_SMM="--use-smm"                # Enable SMM background correction
SMM_K_NEIGHBORS=6                  # KNN neighbors for spatial graph
SMM_SNR_THRESHOLD=2.0              # Minimum SNR to pass filter (INCREASED - aggressive filtering)
SMM_MIN_SIGNAL_FRACTION=0.05       # Min fraction of spots with signal
SMM_MAX_SIGNAL_FRACTION=0.95       # Max fraction of spots with signal
SMM_BETA_INIT=1.0                  # Initial spatial regularization
SMM_MAX_ITER=50                    # Max EM iterations

# NEW: Spatial smoothing parameters (applied AFTER background correction)
# Pipeline: GMM(raw) -> Soft-scale correction (X * P(signal)) -> Smoothing
SMM_SMOOTHING_SIGMA=1.5            # Gaussian smoothing bandwidth
SMM_SMOOTHING_K=15                 # Number of neighbors for smoothing

# NEW: Visualization parameters
VISUALIZE_MARKERS="--visualize-markers"  # Generate spatial diagnostic plots

# NEW: Laplacian smoothing parameters (for spatial coherence in proportion optimization)
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

# Add auto-profiles flags if in auto mode
if [ "$PROFILE_MODE" == "auto" ]; then
    CMD="$CMD --auto-profiles --max-profile-size $MAX_PROFILE_SIZE --discovery-seed $DISCOVERY_SEED"
    CMD="$CMD --morans-k $MORANS_K --morans-i-threshold $MORANS_I_THRESHOLD"
    CMD="$CMD --model-selection $MODEL_SELECTION --cv-folds $CV_FOLDS --min-profiles $MIN_PROFILES"
    # NEW: Reconstruction-based discovery parameters
    CMD="$CMD --selection-method $SELECTION_METHOD"
    CMD="$CMD --min-reconstruction-improvement $MIN_RECON_IMPROVEMENT"
    CMD="$CMD $ABUNDANCE_ADAPTIVE"
    CMD="$CMD --ubiquitous-cv-threshold $UBIQUITOUS_CV_THRESHOLD"
    CMD="$CMD --rare-presence-threshold $RARE_PRESENCE_THRESHOLD"
    CMD="$CMD --redundancy-threshold $REDUNDANCY_THRESHOLD"
    # MIQP-specific parameters (only applied when selection-method=miqp or miqp_hierarchical)
    CMD="$CMD --miqp-lambda-spatial $MIQP_LAMBDA_SPATIAL"
    CMD="$CMD --miqp-lambda-complexity $MIQP_LAMBDA_COMPLEXITY"
    CMD="$CMD --miqp-time-limit $MIQP_TIME_LIMIT"
    CMD="$CMD --miqp-gap $MIQP_GAP"
    # Hierarchical MIQP parameters (only applied when selection-method=miqp_hierarchical)
    CMD="$CMD --miqp-lambda-overlap $MIQP_LAMBDA_OVERLAP"
    CMD="$CMD --miqp-lambda-orphan $MIQP_LAMBDA_ORPHAN"
    CMD="$CMD --miqp-lambda-sparsity $MIQP_LAMBDA_SPARSITY"
    CMD="$CMD --miqp-sparsity-aggregation $MIQP_SPARSITY_AGGREGATION"
    # SMM background correction parameters
    CMD="$CMD $USE_SMM"
    CMD="$CMD --smm-k-neighbors $SMM_K_NEIGHBORS"
    CMD="$CMD --smm-snr-threshold $SMM_SNR_THRESHOLD"
    CMD="$CMD --smm-min-signal-fraction $SMM_MIN_SIGNAL_FRACTION"
    CMD="$CMD --smm-max-signal-fraction $SMM_MAX_SIGNAL_FRACTION"
    CMD="$CMD --smm-beta-init $SMM_BETA_INIT"
    CMD="$CMD --smm-max-iter $SMM_MAX_ITER"
    # Spatial smoothing parameters
    CMD="$CMD --smm-smoothing-sigma $SMM_SMOOTHING_SIGMA"
    CMD="$CMD --smm-smoothing-k $SMM_SMOOTHING_K"
    # Visualization
    CMD="$CMD $VISUALIZE_MARKERS"
    # Laplacian smoothing parameters
    CMD="$CMD --lambda-laplacian $LAMBDA_LAPLACIAN"
    CMD="$CMD --laplacian-k $LAPLACIAN_K"
    echo "  - max_profile_size=$MAX_PROFILE_SIZE"
    echo "  - discovery_seed=$DISCOVERY_SEED"
    echo "  - morans_k=$MORANS_K"
    echo "  - morans_i_threshold=$MORANS_I_THRESHOLD"
    echo "  - model_selection=$MODEL_SELECTION"
    echo "  - cv_folds=$CV_FOLDS"
    echo "  - min_profiles=$MIN_PROFILES"
    echo "  - selection_method=$SELECTION_METHOD"
    echo "  - min_recon_improvement=$MIN_RECON_IMPROVEMENT"
    echo "  - abundance_adaptive=true"
    echo "  - miqp_lambda_spatial=$MIQP_LAMBDA_SPATIAL"
    echo "  - miqp_lambda_complexity=$MIQP_LAMBDA_COMPLEXITY"
    echo "  - miqp_time_limit=$MIQP_TIME_LIMIT"
    echo "  - miqp_gap=$MIQP_GAP"
    echo "  - miqp_lambda_overlap=$MIQP_LAMBDA_OVERLAP"
    echo "  - miqp_lambda_orphan=$MIQP_LAMBDA_ORPHAN"
    echo "  - miqp_lambda_sparsity=$MIQP_LAMBDA_SPARSITY"
    echo "  - miqp_sparsity_aggregation=$MIQP_SPARSITY_AGGREGATION"
    echo "  - use_smm=true"
    echo "  - smm_snr_threshold=$SMM_SNR_THRESHOLD"
    echo "  - smm_k_neighbors=$SMM_K_NEIGHBORS"
    echo "  - smm_smoothing_sigma=$SMM_SMOOTHING_SIGMA"
    echo "  - smm_smoothing_k=$SMM_SMOOTHING_K"
    echo "  - visualize_markers=true"
    echo "  - lambda_laplacian=$LAMBDA_LAPLACIAN"
    echo "  - laplacian_k=$LAPLACIAN_K"
fi

echo "Executing: $CMD"
eval $CMD

echo "Job completed at $(date)"
