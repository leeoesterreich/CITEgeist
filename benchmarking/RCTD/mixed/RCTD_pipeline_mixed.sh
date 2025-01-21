#!/bin/bash
#SBATCH --job-name=RCTD_mx
#SBATCH --output=./logs/RCTD_mixed_%a.out
#SBATCH --error=./logs/RCTD_mixed_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=0-4

cd /bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/RCTD/mixed

# Load R module
module load gcc/12.2.0 r/4.4.0

# Set replicate index from job array ID
REPLICATE_INDEX=$SLURM_ARRAY_TASK_ID

# Define output directory
OUTPUT_DIR="/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/RCTD/mixed"

# Run the R script for the specific replicate
Rscript --vanilla RCTD_pipeline_mixed.R --replicates $REPLICATE_INDEX --output_dir $OUTPUT_DIR

