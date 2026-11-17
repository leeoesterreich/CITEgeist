#!/bin/bash
#SBATCH --job-name=train_proto
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/train_proto_%j.out
#SBATCH --error=logs/train_proto_%j.err

# Prototype Learning Training (Stage 2)
# This script trains projection heads and prototypes for cell type assignment
# using a pre-trained VAE encoder from Stage 1.

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
VAE_CHECKPOINT="${VAE_CHECKPOINT:-${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_vae/vae_final.pt}"
PATCHES_DIR="${PATCHES_DIR:-${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_patches}"
PROPORTIONS_FILE="${PROPORTIONS_FILE:-${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_proportions/proportions.csv}"
OUTPUT_DIR="${OUTPUT_DIR:-${CITEGEIST_REPO_ROOT}/Benchmarking/xenium_benchmarking/CITEgeist/output_prototypes}"

# Number of cell types (default 7 for Xenium benchmark)
N_TYPES="${N_TYPES:-7}"

echo "========================================"
echo "Prototype Learning Training Script"
echo "========================================"
echo "VAE checkpoint: $VAE_CHECKPOINT"
echo "Patches directory: $PATCHES_DIR"
echo "Proportions file: $PROPORTIONS_FILE"
echo "Output directory: $OUTPUT_DIR"
echo "Number of cell types: $N_TYPES"
echo "Device: cuda"
echo "========================================"

# Run prototype training
python -m CITEgeist.model.train_prototypes \
    --vae-checkpoint "$VAE_CHECKPOINT" \
    --patches-dir "$PATCHES_DIR" \
    --proportions-file "$PROPORTIONS_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --n-types "$N_TYPES" \
    --latent-dim 128 \
    --projection-dim 32 \
    --epochs 50 \
    --lr 1e-3 \
    --sinkhorn-temp 0.1 \
    --sinkhorn-iters 50 \
    --device cuda \
    --checkpoint-interval 10 \
    --min-cells-per-spot 3

echo "========================================"
echo "Prototype training completed successfully"
echo "Model saved to: $OUTPUT_DIR/prototypes_final.pt"
echo "========================================"
