#!/bin/bash
#SBATCH --job-name=marker_guidance_test
#SBATCH --output=logs/marker_guidance_%j.out
#SBATCH --error=logs/marker_guidance_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Integration test for marker-guided GEX allocation
# Tests the new use_marker_guidance=True (default) behavior on a Xenium pseudo-Visium region

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

# Run on region 2 (a stable region, avoiding region 1 which can be an outlier)
# --run-gex enables GEX deconvolution, which uses use_marker_guidance=True by default
python Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py \
    --region-id 2 \
    --mode manual \
    --run-gex \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output/marker_guidance_test
