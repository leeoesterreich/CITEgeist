#!/bin/bash
#SBATCH --job-name=fix_pytorch
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --cluster=htc

# Fix PyTorch CUDA support in cell2location_env
# The current PyTorch installation is CPU-only

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
set -e

ENV_PATH="${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/cell2location_env"
MINICONDA_PATH="$(conda info --base)"

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
