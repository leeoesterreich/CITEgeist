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
#SBATCH --array=0-29  # 3 radii * 5 weight pairs * 2 test sets = 30 combinations

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
# Define global/local weight pairs (global_weight,local_weight)
weight_pairs=("0.2,0.8" "0.4,0.6" "0.5,0.5" "0.6,0.4" "0.8,0.2")
TEST_SETS=("mixed" "high_seg")

# Calculate indices for each parameter based on array task ID
test_set_index=$((SLURM_ARRAY_TASK_ID % 2))
combined_index=$((SLURM_ARRAY_TASK_ID / 2))
radius_index=$((combined_index / 5))
weight_index=$((combined_index % 5))

# Get the actual parameter values
radius=${radii[$radius_index]}
weight_pair=${weight_pairs[$weight_index]}
TEST_SET=${TEST_SETS[$test_set_index]}

# Split weight pair into global and local weights
IFS=',' read -r global_weight local_weight <<< "${weight_pair}"

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/${TEST_SET}/"

mkdir -p "$OUTPUT_FOLDER"

echo "Running with parameters:"
echo "  - radius=$radius"
echo "  - global_weight=$global_weight"
echo "  - local_weight=$local_weight"
echo "  - test_set=$TEST_SET"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Run the Python script with these parameters
python tests/test_citegeist_simulated.py \
    --radius $radius \
    --global_enrichment_weight $global_weight \
    --local_enrichment_weight $local_weight \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER \
    --skip_pass2
echo "Job completed at $(date)"