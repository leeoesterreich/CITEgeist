#!/bin/bash
#SBATCH --job-name=c2l_train_p7
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/Cell2Location/slurm/slurm_log/train_ref_protein7_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/Cell2Location/slurm/slurm_log/train_ref_protein7_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=a100
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Train Cell2Location reference model on GSE156632 with 7 protein-GT types

set -e

module load cuda/11.8.0
echo "Loaded CUDA: $(nvcc --version 2>/dev/null | grep release || echo 'CUDA loaded')"

BASE_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking"
C2L_DIR="${BASE_DIR}/Cell2Location"
SRC_DIR="${C2L_DIR}/src"
REF_DATA="${BASE_DIR}/reference_data/GSE156632/processed_protein7/cell2location/reference.h5ad"
OUTPUT_DIR="${C2L_DIR}/reference_model_protein7"

echo "=============================================="
echo "Cell2Location Reference Training (7 Protein-GT Types)"
echo "=============================================="
echo "Start time: $(date)"
echo "Reference data: ${REF_DATA}"
echo "Output directory: ${OUTPUT_DIR}"

if [ ! -f "${REF_DATA}" ]; then
    echo "ERROR: Reference data not found at ${REF_DATA}"
    echo "Please run process_reference_protein7.sh first"
    exit 1
fi

ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/cell2location_env"
MINICONDA_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/miniconda3"
echo "Activating cell2location_env..."
export PATH="${ENV_PATH}/bin:${MINICONDA_PATH}/bin:$PATH"
source activate "${ENV_PATH}"

echo "Python: $(which python)"
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}'); print(f'GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"
echo ""

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
ls -la "${OUTPUT_DIR}/"
