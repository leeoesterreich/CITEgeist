#!/bin/bash
#SBATCH --job-name=laplacian_ablation
#SBATCH --output=logs/laplacian_ablation_%A_%a.out
#SBATCH --error=logs/laplacian_ablation_%A_%a.err
#SBATCH --array=0-9
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Array tasks 0-4: lambda_laplacian=0.0 (NO spatial smoothing)
# Array tasks 5-9: lambda_laplacian=0.1 (WITH spatial smoothing)

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

if [ ${SLURM_ARRAY_TASK_ID} -lt 5 ]; then
    REPLICATE_ID=${SLURM_ARRAY_TASK_ID}
    LAMBDA=0.0
else
    REPLICATE_ID=$((SLURM_ARRAY_TASK_ID - 5))
    LAMBDA=0.1
fi

echo "Running replicate ${REPLICATE_ID} with lambda_laplacian=${LAMBDA}"

python Benchmarking/simulation_benchmarking/CITEgeist/src/test_laplacian_ablation.py \
    --replicate-id ${REPLICATE_ID} \
    --condition mixed \
    --lambda-laplacian ${LAMBDA}
