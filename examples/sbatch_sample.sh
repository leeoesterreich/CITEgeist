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

# Activate conda environment
source /bgfs/alee/LO_LAB/Personal/Alexander_Chang/miniconda3/bin/activate /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/singlecell/
echo "Activated conda environment"

# Load Gurobi module
module load gurobi/11.0.2
echo "Loaded gurobi module"

# Change to working directory
cd /bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
echo "Changed to working directory"


path_to_visium_folder ="/bgfs/alee/LO_LAB/General/Lab_Data/20240510_Neil_SpatialSequencing_AgePatients/SpaceRanger_OUT/MWS23-10650/outs"

/bgfs/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/singlecell/bin/python examples/compute_sample.py --path $path_to_visium_folder --output_folder output