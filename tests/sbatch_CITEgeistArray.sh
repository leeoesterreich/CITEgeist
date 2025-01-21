#!/bin/bash
#SBATCH --job-name=CITEgeist_array
#SBATCH --output=benchmarking_logs/CITEgeist_array_%A_%a.log
#SBATCH --error=benchmarking_logs/CITEgeist_array_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --mem=64G
##SBATCH --mail-type=ALL
##SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=16
#SBATCH --partition=HTC
#SBATCH --array=0-29  # 3 radii * 5 lambda_prior_weights * 2 test sets = 30 combinations

# Activate conda environment
source /bgfs/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/singlecell/
echo "Activated conda environment"

# Load Gurobi module
module load gurobi/11.0.2
echo "Loaded gurobi module"

# Change to working directory
cd /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Define parameter arrays
radii=(4 6 8)
lambda_prior_weights=(1.0 10.0 50.0 100.0 500.0)
TEST_SETS=("mixed" "high_seg")

# Calculate indices for each parameter based on array task ID
# We have 3 radii * 5 lambda_prior_weights * 2 test sets = 30 total combinations
# First, get the test set index (0 or 1)
test_set_index=$((SLURM_ARRAY_TASK_ID % 2))

# Then, get the combined radius and lambda index
combined_index=$((SLURM_ARRAY_TASK_ID / 2))

# Finally, separate radius and lambda indices
radius_index=$((combined_index / 5))  # Integer division by number of lambda values
lambda_index=$((combined_index % 5))  # Modulo by number of lambda values

# Get the actual parameter values
radius=${radii[$radius_index]}
lambda_prior_weight=${lambda_prior_weights[$lambda_index]}
TEST_SET=${TEST_SETS[$test_set_index]}

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/${TEST_SET}/"

mkdir -p "$OUTPUT_FOLDER"

echo "Running with parameters:"
echo "  - radius=$radius"
echo "  - lambda_prior_weight=$lambda_prior_weight"
echo "  - test_set=$TEST_SET"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Run the Python script with these parameters
python tests/test_citegeist_simulated.py \
    --radius $radius \
    --lambda_prior_weight $lambda_prior_weight \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER \
    --skip_pass2
echo "Job completed at $(date)"