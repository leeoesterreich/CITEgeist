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
#SBATCH --array=1-14%3

# Activate conda environment
conda activate CITEgeist_env
echo "Activated conda environment"

# CHANGE THESE VALUES FOR REPRODUCIBILITY ON YOUR OWN COMPUTER
CITEGEIST_FOLDER="/path/to/local/install/CITEgeist/CITEgeist"
GUROBI_LICENSE_FILE="/path/to/your/gurobi.lic"

# Change to working directory
cd $CITEGEIST_FOLDER
echo "Changed to working directory"

# Load Gurobi module (lmod)
module load gurobi/11.0.2
echo "Loaded gurobi module"

# Get the path for this array job
DATA_FOLDER="./data/GEO_data"
path_to_visium_folder=$(sed -n "${SLURM_ARRAY_TASK_ID}p" data/sample_paths.txt)


echo "Processing sample: $path_to_visium_folder"

python examples/compute_sample.py --path "$DATA_FOLDER/$path_to_visium_folder" --output_folder output --radius 400 --min_counts 100 --gurobi_license "$GUROBI_LICENSE_FILE"
