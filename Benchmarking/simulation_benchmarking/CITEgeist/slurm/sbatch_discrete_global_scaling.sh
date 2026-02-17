#!/bin/bash
#SBATCH --job-name=sim_global
#SBATCH --cluster=gpu
#SBATCH --partition=l40s
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/sim_global_%A_%a.out
#SBATCH --error=logs/sim_global_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

# Test global scaling mode on simulation benchmark
# This should dramatically improve proportion correlation for mixed condition

CONDITION=${1:-mixed}
MODE=${2:-dapi}
REPLICATE_ID=$SLURM_ARRAY_TASK_ID

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate ~/alc376_bgfs/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

OUTPUT_DIR=Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global_scaling

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id $REPLICATE_ID \
    --condition $CONDITION \
    --mode $MODE \
    --output-dir $OUTPUT_DIR \
    --use-gpu \
    --max-em-iterations 20 \
    --scale-mode global

echo "Completed replicate $REPLICATE_ID, condition $CONDITION, scale_mode=global"
