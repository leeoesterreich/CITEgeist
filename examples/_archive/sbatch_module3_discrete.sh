#!/bin/bash
#SBATCH --job-name=module3_discrete
#SBATCH --output=slurm_log/module3_discrete_%A_%a.out
#SBATCH --error=slurm_log/module3_discrete_%A_%a.err
#SBATCH --array=0-13
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Module 3 Discrete Cell Assignment (IQP + EM)
# Alternative to continuous proportion estimation
# Requires nuclei counts from Cellpose segmentation

# Create log directory
mkdir -p slurm_log

# Patient samples array (14 total)
SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

# Get current sample from array index
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "Module 3 Discrete Cell Assignment"
echo "Processing: ${SAMPLE}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Date: $(date)"
echo "=========================================="

# Change to project directory
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Load required modules
module load gurobi/12.0.3

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

# Output directory
OUTPUT_DIR="output/module3_discrete/${SAMPLE}"
NUCLEI_DIR="output/segmentation/${SAMPLE}"

# Check if nuclei counts exist
if [[ ! -f "${NUCLEI_DIR}/nuclei_counts.csv" ]]; then
    echo "ERROR: Nuclei counts not found at ${NUCLEI_DIR}/nuclei_counts.csv"
    echo "Run Cellpose segmentation first using segmentation.py"
    exit 1
fi

# Run discrete cell assignment
python examples/run_module3_discrete.py \
    --sample "${SAMPLE}" \
    --output-dir "${OUTPUT_DIR}" \
    --nuclei-file "${NUCLEI_DIR}/nuclei_counts.csv"

echo "=========================================="
echo "Completed: ${SAMPLE}"
echo "Date: $(date)"
echo "=========================================="
