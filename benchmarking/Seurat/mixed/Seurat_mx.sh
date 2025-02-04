#! /usr/bin/bash
#SBATCH --job-name=Seurat_mx_deconv
#SBATCH --output=./logs/Seurat_deconv_mx_%a.out
#SBATCH --error=./logs/Seurat_deconv_mx_%a.err
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --array=0-4  # Adjust this based on the number of replicates

# Load R module
module load gcc/12.2.0 r/4.4.0

# Define input directory
INPUT_DIR="/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/replicates/mixed/h5ad_objects"
OUTPUT_H5SEURAT_DIR="/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/replicates/mixed/h5seurat"
OUTPUT_DIR="/bgfs/alee/LO_LAB/Personal/Brent_Schlegel/Projects/Wu_Visium/Simulations/scCube_12k/Seurat/mixed/"

mkdir -pv "$OUTPUT_H5SEURAT_DIR"

# Get list of h5ad files
FILES=($(ls $INPUT_DIR/*GEX.h5ad))

# Ensure we do not go out of bounds
if [ "$SLURM_ARRAY_TASK_ID" -ge "${#FILES[@]}" ]; then
    echo "Error: SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID is out of range for available files."
    exit 1
fi

# Get the specific file for this array job
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Extract the replicate name from the file
BASENAME=$(basename "$FILE" _GEX.h5ad)

# Capture start time
START_TIME=$(date +%s)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running Seurat deconvolution for replicate: $BASENAME"

# Run the Seurat R script for this replicate
Rscript --vanilla Seurat_mx.R "$FILE" "$OUTPUT_DIR" "$BASENAME"

# Capture end time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_MINUTES=$(echo "scale=2; $RUNTIME / 60" | bc)

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed Seurat deconvolution for replicate: $BASENAME"
echo "SEURAT_TOTAL_RUNTIME: $BASENAME took $RUNTIME seconds ($RUNTIME_MINUTES minutes)"

# Move converted h5Seurat files
mv "$INPUT_DIR"/*.h5seurat "$OUTPUT_H5SEURAT_DIR"
