#!/bin/bash
#SBATCH --job-name=CITEgeist_array
#SBATCH --output=benchmarking_logs/CITEgeist_array_%A_%a.log
#SBATCH --error=benchmarking_logs/CITEgeist_array_%A_%a.log
#SBATCH --time=72:00:00
##SBATCH --mail-type=ALL
##SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --partition=HTC
#SBATCH --array=0-11 # 2 test sets, 2 alpha, 1 lambda_reg, 3 max_y_change = 12 combinations

# Activate conda environment
source /bgfs/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/singlecell/
echo "Activated conda environment"

# Load Gurobi module
module load gurobi/11.0.2
echo "Loaded gurobi module"

# Change to working directory
cd /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Define parameter arrays (fixed values)
lambda_reg=(1)
alpha_elastic=(0.7 0.9)
max_y_change=(0.2 0.4 0.8)
TEST_SETS=("mixed" "high_seg")

# Calculate indices for the test set, lambda_reg, alpha_elastic, and max_y_change based on array task ID
test_set_index=$((SLURM_ARRAY_TASK_ID / (1 * 2 * 3)))
lambda_reg_index=$(((SLURM_ARRAY_TASK_ID / (2 * 3)) % 1))
alpha_elastic_index=$(((SLURM_ARRAY_TASK_ID / 3) % 2))
max_y_change_index=$((SLURM_ARRAY_TASK_ID % 3))
# Get the test set, lambda_reg, and alpha_elastic values
TEST_SET=${TEST_SETS[$test_set_index]}
lambda_reg=${lambda_reg[$lambda_reg_index]}
alpha_elastic=${alpha_elastic[$alpha_elastic_index]}
max_y_change=${max_y_change[$max_y_change_index]}

INPUT_FOLDER="replicates/${TEST_SET}/"
OUTPUT_FOLDER="test_results/${TEST_SET}/"



mkdir -p "$OUTPUT_FOLDER"

echo "Running with parameters:"
echo "  - lambda_reg=$lambda_reg"
echo "  - alpha_elastic=$alpha_elastic"
echo "  - max_y_change=$max_y_change"
echo "  - test_set=$TEST_SET"
echo "Input folder: $INPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Run the Python script with these parameters
/bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/singlecell/bin/python tests/test_citegeist_simulated.py \
    --radius 4 \
    --lambda_reg $lambda_reg \
    --alpha_elastic $alpha_elastic \
    --max_y_change $max_y_change \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER 
    
    
    

echo "Job completed at $(date)"