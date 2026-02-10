#!/bin/bash
#SBATCH --job-name=CARD_hs
#SBATCH --output=./logs/CARD_hs_%a.out
#SBATCH --error=./logs/CARD_hs_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Run CARD deconvolution on high_seg simulation (REFERENCE mode)
# Array job: processes replicates 0-4
#
# For reference-free mode, change --mode to reference_free
# and ensure marker_genes_simulation.csv exists

set -e

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[$START_TIME] Job started for replicate $SLURM_ARRAY_TASK_ID"

# Directories
SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CARD/high_seg"
OUTPUT_DIR="${SCRIPT_DIR}"

# Create log directory
mkdir -p "${SCRIPT_DIR}/logs"

cd "${SCRIPT_DIR}"

# Environment
CARD_ENV="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CARD_env"

# Activate CARD environment
echo "Activating CARD_env..."
export PATH="${CARD_ENV}/bin:$PATH"
source activate "${CARD_ENV}"

echo "R: $(which R)"
echo "R version: $(R --version | head -1)"

# Set replicate index from array ID
REPLICATE_INDEX=$SLURM_ARRAY_TASK_ID

# Run CARD (reference mode)
Rscript --vanilla CARD_pipeline_high_seg.R \
    --replicates $REPLICATE_INDEX \
    --output_dir $OUTPUT_DIR \
    --mode reference

# Record end time
END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$END_TIME] Job completed for replicate $SLURM_ARRAY_TASK_ID"
echo "CARD_TOTAL_RUNTIME: Replicate $SLURM_ARRAY_TASK_ID took $RUNTIME seconds ($RUNTIME_MINUTES minutes)"
