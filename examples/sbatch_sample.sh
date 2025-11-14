#!/bin/bash
#SBATCH --job-name=CITEgeist_sample
#SBATCH --output=benchmarking_logs/CITEgeist_sample_%A_%a.log
#SBATCH --error=benchmarking_logs/CITEgeist_sample_%A_%a.log
#SBATCH --time=72:00:00
#SBATCH --mem=64G
##SBATCH --mail-type=ALL
##SBATCH --mail-user=alc376@pitt.edu
#SBATCH --cluster=htc
#SBATCH --cpus-per-task=16
#SBATCH --partition=HTC
#SBATCH --array=1-14

# Activate conda environment
source /bgfs/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/
echo "Activated conda environment"


# CHANGE THIS FOR REPRODUCIBILITY ON YOUR OWN COMPUTER
DATA_FOLDER="/bgfs/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files"


# Load Gurobi module
module load gurobi/11.0.2
echo "Loaded gurobi module"

# Change to working directory
cd /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"

# Get the path for this array job
path_to_visium_folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" examples/sample_paths.txt)


echo "Processing sample: $path_to_visium_folder"

python examples/compute_sample.py --path "$DATA_FOLDER/$path_to_visium_folder" --output_folder output --radius 400 --min_counts 100
