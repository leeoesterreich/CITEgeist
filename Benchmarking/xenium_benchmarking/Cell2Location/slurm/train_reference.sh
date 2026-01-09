#!/bin/bash
#SBATCH --job-name=c2l_train_ref
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/Cell2Location/slurm/slurm_log/train_ref_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/Cell2Location/slurm/slurm_log/train_ref_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Train Cell2Location reference model on GSE156632
# This is a one-time job that creates the reference signatures

set -e

# Load CUDA module for GPU support
module load cuda/11.8.0
echo "Loaded CUDA: $(nvcc --version 2>/dev/null | grep release || echo 'CUDA loaded')"

# Directories (use absolute paths - SCRIPT_DIR resolves incorrectly in SLURM)
BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
C2L_DIR="${BASE_DIR}/Cell2Location"
SRC_DIR="${C2L_DIR}/src"
REF_DATA="${BASE_DIR}/reference_data/GSE156632/processed_granular/cell2location/reference.h5ad"
OUTPUT_DIR="${C2L_DIR}/reference_model_granular"

echo "=============================================="
echo "Cell2Location Reference Training"
echo "=============================================="
echo "Start time: $(date)"
echo "Reference data: ${REF_DATA}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if reference data exists
if [ ! -f "${REF_DATA}" ]; then
    echo "ERROR: Reference data not found at ${REF_DATA}"
    echo "Please run reference_data/GSE156632/slurm/process_reference.sh first"
    exit 1
fi

# Activate Cell2Location environment (use export PATH + source activate pattern)
ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/cell2location_env"
MINICONDA_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/miniconda3"
echo "Activating cell2location_env from ${ENV_PATH}..."
export PATH="${ENV_PATH}/bin:${MINICONDA_PATH}/bin:$PATH"
source activate "${ENV_PATH}"

echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo ""

# Check GPU
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}'); print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"
echo ""

# Run training
python "${SRC_DIR}/train_reference.py" \
    --ref-path "${REF_DATA}" \
    --output-dir "${OUTPUT_DIR}" \
    --batch-key "sample" \
    --labels-key "cell_type" \
    --max-epochs 250 \
    --filter-genes

echo ""
echo "=============================================="
echo "Training Complete!"
echo "=============================================="
echo "End time: $(date)"
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/"
