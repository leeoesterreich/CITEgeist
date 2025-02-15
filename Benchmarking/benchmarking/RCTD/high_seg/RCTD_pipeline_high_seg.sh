#!/bin/bash
#SBATCH --job-name=RCTD_hs
#SBATCH --output=./logs/RCTD_hs_%a.out
#SBATCH --error=./logs/RCTD_hs_%a.err
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=0-4

# Record the start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[$START_TIME] Job started for replicate $SLURM_ARRAY_TASK_ID" | tee -a ./logs/RCTD_hs_runtime.log

cd /bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/RCTD/high_seg

# Load R module
module load gcc/12.2.0 r/4.4.0

# Set replicate index from job array ID
REPLICATE_INDEX=$SLURM_ARRAY_TASK_ID

# Define output directory
OUTPUT_DIR="/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/RCTD/high_seg"

# Run the R script for the specific replicate
Rscript --vanilla RCTD_pipeline_high_seg.R --replicates $REPLICATE_INDEX --output_dir $OUTPUT_DIR

# Record the end time
END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')

# Calculate total runtime
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$END_TIME] Job completed for replicate $SLURM_ARRAY_TASK_ID" | tee -a ./logs/RCTD_hs_runtime.log
echo "RCTD_TOTAL_RUNTIME: Replicate $SLURM_ARRAY_TASK_ID took $RUNTIME seconds ($RUNTIME_MINUTES minutes)" | tee -a ./logs/RCTD_hs_runtime.log
