#!/bin/bash
#SBATCH --job-name=train_vae
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/train_vae_%j.out
#SBATCH --error=logs/train_vae_%j.err

# VAE Training for Nucleus Patch Representation Learning (Stage 1)
# This script trains a Variational Autoencoder on pre-extracted nucleus patches.
# The trained encoder will be used in Stage 2 for prototype learning.

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

# Create logs directory if not exists
mkdir -p logs

# Load modules and activate environment
eval "$(conda shell.bash hook)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

# Change to project root
cd ${CITEGEIST_REPO_ROOT}

# Define paths (modify as needed)
PATCHES_DIR="${PATCHES_DIR:-${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_patches}"
OUTPUT_DIR="${OUTPUT_DIR:-${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_vae}"

echo "========================================"
echo "VAE Training Script"
echo "========================================"
echo "Patches directory: $PATCHES_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Device: cuda"
echo "========================================"

# Run VAE training
python -m CITEgeist.model.train_vae \
    --patches-dir "$PATCHES_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --in-channels 3 \
    --latent-dim 128 \
    --batch-size 64 \
    --epochs 100 \
    --lr 1e-4 \
    --beta 0.5 \
    --device cuda \
    --checkpoint-interval 10 \
    --num-workers 4

echo "========================================"
echo "VAE training completed successfully"
echo "Model saved to: $OUTPUT_DIR/vae_final.pt"
echo "========================================"
