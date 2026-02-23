#!/bin/bash
#SBATCH --job-name=cellpose_rerun
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=1,3
#SBATCH --output=logs/cellpose_rerun_%A_%a.out
#SBATCH --error=logs/cellpose_rerun_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load CUDA before conda to ensure PyTorch finds it
module load cuda/12.1.1
module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

# Verify GPU is available
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}, Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"

REGION_ID=$SLURM_ARRAY_TASK_ID
OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py \
    --region $REGION_ID \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20

echo "Completed region $REGION_ID"
