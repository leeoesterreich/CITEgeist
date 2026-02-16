#!/bin/bash
#SBATCH --job-name=cellpose_r3
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cellpose_region3_%j.out
#SBATCH --error=logs/cellpose_region3_%j.err
#SBATCH --exclude=gpu-n66
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Load CUDA before conda to ensure PyTorch finds it
module load cuda/12.1.1
module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

# Verify GPU is available
python -c "import torch; print(f'CUDA available: {torch.cuda.is_available()}, Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}')"

# Exit if CUDA not available
python -c "import torch; assert torch.cuda.is_available(), 'CUDA not available!'"

OUTPUT_DIR=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py \
    --region 3 \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20

echo "Completed region 3"
