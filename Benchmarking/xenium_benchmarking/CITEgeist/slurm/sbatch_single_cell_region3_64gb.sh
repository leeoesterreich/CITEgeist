#!/bin/bash
#SBATCH --job-name=sc_res_r3
#SBATCH --output=logs/sc_resolution_r3_%j.out
#SBATCH --error=logs/sc_resolution_r3_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Single-cell resolution benchmark for region 3 with increased memory (64GB)
# Region 3 has a larger image that requires more RAM for Cellpose segmentation

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_single_cell_resolution.py \
    --region 3 \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/single_cell_resolution \
    --use-gpu \
    --spot-diameter-um 55.0
