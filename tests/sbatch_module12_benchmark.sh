#!/bin/bash
#SBATCH --job-name=M12_bench
#SBATCH --output=benchmarking_logs/M12_bench_%A_%a.log
#SBATCH --error=benchmarking_logs/M12_bench_%A_%a.log
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --partition=HTC
#SBATCH --array=0-11  # 2 test sets × 2 min_marginal_gain × 3 max_y_change = 12 combinations

# Module 1-2 Benchmark with looser profile selection
# Tests impact of min_marginal_gain on profile discovery

# Activate conda environment
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"

# Load Gurobi module
module load gurobi/12.0.3
echo "Loaded gurobi module"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Fixed parameters
RADIUS=4.0
LAMBDA_REG=1.0
ALPHA_ELASTIC=0.7

# Module 1-2 parameters
TOP_K=3
DISCOVERY_SEED=1234
MORANS_K=8
SMOOTH_K=6
N_PERMUTATIONS=999
FDR_THRESHOLD=0.05
GMM_SNR_THRESHOLD=1.0

# Module 2c: Test different min_marginal_gain values
VARIANCE_TARGET=0.90
# Key change: Test 0.005 (0.5%) vs 0.001 (0.1%) vs 0.0001 (0.01%)
MIN_MARGINAL_GAINS=(0.001 0.0001)  # Much lower = include more profiles
MAX_PROFILES=15

# Laplacian and beta parameters
LAMBDA_LAPLACIAN=0.1
LAPLACIAN_K=8
PER_MARKER_BETA=true
BETA_MIN=0.1
BETA_MAX=2.0

# Test sets and max_y_change
TEST_SETS=("mixed" "high_seg")
MAX_Y_CHANGES=(0.2 0.4 0.8)

# Calculate indices
n_max_y=3
n_marginal=2
n_test_sets=2

test_set_index=$((SLURM_ARRAY_TASK_ID / (n_marginal * n_max_y)))
marginal_index=$(((SLURM_ARRAY_TASK_ID / n_max_y) % n_marginal))
max_y_index=$((SLURM_ARRAY_TASK_ID % n_max_y))

# Get parameter values
TEST_SET=${TEST_SETS[$test_set_index]}
MIN_MARGINAL_GAIN=${MIN_MARGINAL_GAINS[$marginal_index]}
MAX_Y_CHANGE=${MAX_Y_CHANGES[$max_y_index]}

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/${TEST_SET}/auto_profiles_v2/"

mkdir -p "$OUTPUT_FOLDER"
mkdir -p benchmarking_logs

echo "Running Module 1-2 benchmark with parameters:"
echo "  - test_set=$TEST_SET"
echo "  - min_marginal_gain=$MIN_MARGINAL_GAIN"
echo "  - max_y_change=$MAX_Y_CHANGE"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Build command
CMD="python tests/test_citegeist_simulated.py \
    --radius $RADIUS \
    --lambda_reg $LAMBDA_REG \
    --alpha_elastic $ALPHA_ELASTIC \
    --max_y_change $MAX_Y_CHANGE \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER \
    --auto-profiles \
    --top-k $TOP_K \
    --discovery-seed $DISCOVERY_SEED \
    --morans-k $MORANS_K \
    --smooth-k $SMOOTH_K \
    --n-permutations $N_PERMUTATIONS \
    --fdr-threshold $FDR_THRESHOLD \
    --gmm-snr-threshold $GMM_SNR_THRESHOLD \
    --variance-target $VARIANCE_TARGET \
    --min-marginal-gain $MIN_MARGINAL_GAIN \
    --lambda-laplacian $LAMBDA_LAPLACIAN \
    --laplacian-k $LAPLACIAN_K \
    --per-marker-beta \
    --beta-min $BETA_MIN \
    --beta-max $BETA_MAX"

echo "Executing: $CMD"
eval $CMD

echo "Job completed at $(date)"
