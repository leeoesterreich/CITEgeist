#!/bin/bash
#SBATCH --job-name=stages456_rerun
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --output=output/morphology_assignment_v2/logs/stages456_%j.out
#SBATCH --error=output/morphology_assignment_v2/logs/stages456_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

set -eo pipefail

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

set -u

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd "$REPO_ROOT"

QP_DIR="output/module3_cuopt_qp"
V2_DIR="output/morphology_assignment_v2"
V1_DIR="output/morphology_assignment"

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

POOLED_CKPT="$V2_DIR/pooled_mil_checkpoint.pt"

echo "=== Stages 4-6 rerun (MIL checkpoint already exists) ==="
echo "Checkpoint: $POOLED_CKPT"
ls -la "$POOLED_CKPT"

for SAMPLE in "${SAMPLES[@]}"; do
    echo ""
    echo "--- Processing $SAMPLE ---"

    SAMPLE_OUT="$V2_DIR/$SAMPLE"
    mkdir -p "$SAMPLE_OUT"

    python examples/run_morphology_assignment.py \
        --sample "$SAMPLE" \
        --stages 4,5,6 \
        --module3-dir "$QP_DIR" \
        --embeddings-dir "$V1_DIR/$SAMPLE" \
        --output-dir-override "$SAMPLE_OUT" \
        --pooled-checkpoint "$POOLED_CKPT" \
        --device cuda

    echo "  $SAMPLE complete."
done

echo ""
echo "=== All 12 samples processed ==="
echo "Output: $V2_DIR/"
