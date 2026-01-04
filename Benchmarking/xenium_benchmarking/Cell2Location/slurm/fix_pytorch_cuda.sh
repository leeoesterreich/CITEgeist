#!/bin/bash
#SBATCH --job-name=fix_pytorch
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/Cell2Location/slurm/slurm_log/fix_pytorch_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/Cell2Location/slurm/slurm_log/fix_pytorch_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc

# Fix PyTorch CUDA support in cell2location_env
# The current PyTorch installation is CPU-only

set -e

ENV_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/cell2location_env"
MINICONDA_PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/miniconda3"

echo "=============================================="
echo "Fixing PyTorch CUDA support in cell2location_env"
echo "=============================================="
echo "Start time: $(date)"

# Activate environment
export PATH="${ENV_PATH}/bin:${MINICONDA_PATH}/bin:$PATH"
source activate "${ENV_PATH}"

echo "Current PyTorch version:"
python -c "import torch; print(f'Version: {torch.__version__}, CUDA: {torch.version.cuda}')"

echo ""
echo "Reinstalling PyTorch with CUDA 11.8 support..."
# Uninstall current CPU-only PyTorch
pip uninstall -y torch torchvision torchaudio 2>/dev/null || true

# Install PyTorch with CUDA 11.8 support
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

echo ""
echo "Verifying PyTorch CUDA support:"
python -c "import torch; print(f'Version: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}'); print(f'CUDA version: {torch.version.cuda}')"

echo ""
echo "=============================================="
echo "Fix Complete!"
echo "=============================================="
echo "End time: $(date)"
