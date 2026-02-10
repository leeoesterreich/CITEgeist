#!/bin/bash
#SBATCH --job-name=CARD_mx_rf
#SBATCH --output=./logs/CARD_mx_reffree_%a.out
#SBATCH --error=./logs/CARD_mx_reffree_%a.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --array=0-4
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Run CARD deconvolution on mixed simulation (REFERENCE-FREE mode)
# Array job: processes replicates 0-4
#
# Prerequisites: marker_genes_simulation.csv must exist in parent directory

set -e

# Record start time
START_TIMESTAMP=$(date +%s)
START_TIME=$(date +'%Y-%m-%d %H:%M:%S')
echo "[$START_TIME] Job started for replicate $SLURM_ARRAY_TASK_ID (ref-free)"

# Directories
SCRIPT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CARD/mixed"
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

# Run CARD (reference-free mode)
Rscript --vanilla CARD_pipeline_mixed.R \
    --replicates $REPLICATE_INDEX \
    --output_dir $OUTPUT_DIR \
    --mode reference_free

# Record end time
END_TIMESTAMP=$(date +%s)
END_TIME=$(date +'%Y-%m-%d %H:%M:%S')
RUNTIME=$((END_TIMESTAMP - START_TIMESTAMP))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$END_TIME] Job completed for replicate $SLURM_ARRAY_TASK_ID (ref-free)"
echo "CARD_REFFREE_TOTAL_RUNTIME: Replicate $SLURM_ARRAY_TASK_ID took $RUNTIME seconds ($RUNTIME_MINUTES minutes)"
