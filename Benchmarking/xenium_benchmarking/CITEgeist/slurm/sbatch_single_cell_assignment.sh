#!/bin/bash
#SBATCH --job-name=sc_assign
#SBATCH --output=logs/sc_assign_%A_%a.out
#SBATCH --error=logs/sc_assign_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell assignment benchmark
# Runs hybrid CITEgeist + constrained assignment to produce single-cell h5ad
#
# Modes:
#   random:     Constrained random shuffle within counts (baseline, ~22% accuracy)
#   morphology: Full 7-class morphology Hungarian (~46% accuracy)
#   xgboost:    VAE + morphology XGBoost classifier (~50%+ accuracy)
#
# Usage:
#   MODE=morphology sbatch sbatch_single_cell_assignment.sh
#   MODE=random sbatch sbatch_single_cell_assignment.sh
#   MODE=xgboost sbatch sbatch_single_cell_assignment.sh

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Default to 'morphology' mode if not specified
MODE=${MODE:-morphology}

# XGBoost model paths (only used for xgboost mode)
XGBOOST_MODEL="Benchmarking/xenium_benchmarking/CITEgeist/output/xgboost_combined/xgboost_model.pkl"
VAE_CHECKPOINT="Benchmarking/xenium_benchmarking/CITEgeist/output/vae_masked/vae/vae_final.pt"

echo "Running single-cell assignment benchmark"
echo "  Region: ${SLURM_ARRAY_TASK_ID}"
echo "  Mode: ${MODE}"
echo "  Date: $(date)"

# Build command arguments
CMD_ARGS="--region ${SLURM_ARRAY_TASK_ID}"
CMD_ARGS="${CMD_ARGS} --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_${MODE}"
CMD_ARGS="${CMD_ARGS} --use-gpu"
CMD_ARGS="${CMD_ARGS} --lambda-laplacian 0.1"
CMD_ARGS="${CMD_ARGS} --spot-diameter-um 55.0"
CMD_ARGS="${CMD_ARGS} --single-cell ${MODE}"
CMD_ARGS="${CMD_ARGS} --patches-dir Benchmarking/xenium_benchmarking/CITEgeist/output/patches_v2"

# Add XGBoost-specific args for xgboost mode
if [ "${MODE}" == "xgboost" ]; then
    CMD_ARGS="${CMD_ARGS} --xgboost-model ${XGBOOST_MODEL}"
    CMD_ARGS="${CMD_ARGS} --vae-checkpoint ${VAE_CHECKPOINT}"
    echo "  XGBoost model: ${XGBOOST_MODEL}"
    echo "  VAE checkpoint: ${VAE_CHECKPOINT}"
fi

# Run hybrid benchmark with single-cell output
python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_hybrid_cellpose.py ${CMD_ARGS}

echo "Benchmark complete: $(date)"
