#!/bin/bash
#SBATCH --job-name=pooled_mil_rerun
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=output/morphology_assignment_v2/logs/mil_rerun_%j.out
#SBATCH --error=output/morphology_assignment_v2/logs/mil_rerun_%j.err
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

mkdir -p "$V2_DIR/logs"

# ============================================================
# Step 0: Poll for Phase 1 completion (all 12 .done markers)
# ============================================================
echo "=== Polling for Phase 1 completion ==="
MAX_WAIT=7200  # 2 hours
INTERVAL=30
ELAPSED=0

while true; do
    DONE_COUNT=0
    for s in "${SAMPLES[@]}"; do
        if [ -f "$QP_DIR/$s/.done" ]; then
            DONE_COUNT=$((DONE_COUNT + 1))
        fi
    done

    echo "  $(date): $DONE_COUNT/12 samples complete"

    if [ "$DONE_COUNT" -eq 12 ]; then
        echo "All 12 samples complete!"
        break
    fi

    if [ "$ELAPSED" -ge "$MAX_WAIT" ]; then
        echo "ERROR: Timed out after ${MAX_WAIT}s. Only $DONE_COUNT/12 done."
        exit 1
    fi

    sleep "$INTERVAL"
    ELAPSED=$((ELAPSED + INTERVAL))
done

# ============================================================
# Step 1: Train pooled MIL on cuOPT QP proportions
# ============================================================
echo ""
echo "=== Step 1: Pooled MIL training ==="

python examples/train_pooled_mil.py \
    --prop-dir "$QP_DIR" \
    --output-dir "$V2_DIR" \
    --epochs 300 \
    --lr 1e-4 \
    --device cuda

echo "Pooled MIL training complete."

# ============================================================
# Step 2: Per-sample Stages 4-6
# ============================================================
echo ""
echo "=== Step 2: Per-sample Stages 4-6 ==="

POOLED_CKPT="$V2_DIR/pooled_mil_checkpoint.pt"

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
