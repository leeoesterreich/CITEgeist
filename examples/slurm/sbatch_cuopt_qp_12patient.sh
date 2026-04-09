#!/bin/bash
#SBATCH --job-name=cuopt_qp_12patient
#SBATCH --array=0-11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=output/module3_cuopt_qp/logs/qp_%A_%a.out
#SBATCH --error=output/module3_cuopt_qp/logs/qp_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

# Activate conda
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

set -u

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd "$REPO_ROOT"

# Map array task ID to sample name
SAMPLES=(
    "HCC22-088-P1-S1"
    "HCC22-088-P1-S2"
    "HCC22-088-P2-S1"
    "HCC22-088-P2-S2"
    "HCC22-088-P3-S1_A"
    "HCC22-088-P3-S2"
    "HCC22-088-P4-S1"
    "HCC22-088-P4-S2_1i_rep"
    "HCC22-088-P5-S1"
    "HCC22-088-P5-S2_F_rep"
    "HCC22-088-P6-S1"
    "HCC22-088-P6-S2_D"
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# Create log dir
mkdir -p output/module3_cuopt_qp/logs

echo "=== Array task $SLURM_ARRAY_TASK_ID: $SAMPLE ==="
echo "Node: $(hostname), GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"

python examples/run_cuopt_qp_patient.py --sample "$SAMPLE"

echo "=== Complete: $SAMPLE ==="
