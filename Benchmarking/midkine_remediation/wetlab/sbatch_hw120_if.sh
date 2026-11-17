#!/bin/bash
# HW120 midkine immunofluorescence quantification (analyze_hw120_if.py, Cellpose 4.0.9). Needs a GPU for Cellpose.
# Usage: IF_DATA_ROOT=/path/to/raw_zstack_folders [IF_OUT=/path/to/output] gpu-race-submit.sh sbatch_hw120_if.sh
# Do not put --partition or --cluster here; the GPU wrapper sets them.
#SBATCH --job-name=hw120_if
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=hw120_if_%j.out
#SBATCH --error=hw120_if_%j.err

set -euo pipefail

: "${IF_DATA_ROOT:?set IF_DATA_ROOT to the folder of raw z-stack TIFFs (one sub-folder per condition)}"
HERE="${HERE_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}}"  # submit from this folder, or set HERE_DIR
echo "namespace: ${GPU_RACE_NAMESPACE:-none} partition: ${GPU_RACE_PARTITION:-none}"
echo "=== HW120 MDK IF Analysis ==="
echo "Node: $(hostname)"
echo "Start: $(date)"

# environment: Python with cellpose==4.0.9, scikit-image, tifffile, seaborn, scipy, pandas (see README.md)
python -u "$HERE/analyze_hw120_if.py" --data-root "$IF_DATA_ROOT" ${IF_OUT:+--output-dir "$IF_OUT"}

echo "End: $(date)"
