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
#SBATCH --array=0-39  # 4 lambda_reg * 5 alpha_elastic * 2 test sets = 40 combinations

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
lambda_reg_values=(0.0001 0.001 0.01 0.1)
alpha_elastic_values=(0.0 0.2 0.5 0.8 1.0)  # 0=L2, 1=L1
TEST_SETS=("mixed" "high_seg")

# Calculate indices for each parameter based on array task ID
test_set_index=$((SLURM_ARRAY_TASK_ID % 2))
combined_index=$((SLURM_ARRAY_TASK_ID / 2))
lambda_index=$((combined_index / 5))
alpha_index=$((combined_index % 5))

# Get the actual parameter values
lambda_reg=${lambda_reg_values[$lambda_index]}
alpha_elastic=${alpha_elastic_values[$alpha_index]}
TEST_SET=${TEST_SETS[$test_set_index]}

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/${TEST_SET}/"

mkdir -p "$OUTPUT_FOLDER"

echo "Running with parameters:"
echo "  - lambda_reg=$lambda_reg"
echo "  - alpha_elastic=$alpha_elastic"
echo "  - test_set=$TEST_SET"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Run the Python script with these parameters
python tests/test_citegeist_simulated.py \
    --radius 4 \
    --lambda_reg $lambda_reg \
    --alpha_elastic $alpha_elastic \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER \
    --skip_pass2
echo "Job completed at $(date)"