#!/bin/bash
#SBATCH --job-name=CITEgeist_array
#SBATCH --output=logs/CITEgeist_array_%A_%a.log
#SBATCH --error=logs/CITEgeist_array_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=16
#SBATCH --array=0-3  # 3 radii = 3 combinations (0-44)

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
radii=(4 6)


# Calculate indices for each parameter based on array task ID
radius_idx=$(( SLURM_ARRAY_TASK_ID / (${#alphas[@]} * ${#gex_lambdas[@]}) ))

# Get parameter values
radius=${radii[$radius_idx]}


TEST_SET="mixed"

INPUT_FOLDER="replicates/${TEST_SET}"/

OUTPUT_FOLDER="test_results/${TEST_SET}"/

mkdir -p "$OUTPUT_FOLDER"


echo "Running with parameters: radius=$radius, input_folder=$INPUT_FOLDER, output_folder=$OUTPUT_FOLDER"

# Run the Python script with these parameters
python CITEgeistBenchmarkGridSearchArray.py \
    --radius $radius \
    --input_folder $INPUT_FOLDER \
    --output_folder $OUTPUT_FOLDER


echo "Job completed at $(date)"