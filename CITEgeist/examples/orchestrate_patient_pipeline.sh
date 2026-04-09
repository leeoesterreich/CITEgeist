#!/bin/bash
# Orchestrate the 5-phase patient sample single-cell assignment pipeline.
#
# Cross-cluster: GPU (l40s) and HTC are on different clusters at Pitt CRC.
# SLURM --dependency=afterok does NOT work across clusters.
# GPU phases use intra-cluster deps; HTC phases use delayed start + marker checks.

set -eo pipefail

EXAMPLES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${EXAMPLES_DIR}/.."   # cd to CITEgeist/ package root

echo "=== Patient Pipeline Orchestrator ==="
echo "Working directory: $(pwd)"

mkdir -p output/patient_pipeline/logs

# Helper: strip cluster suffix from --parsable output
strip_cluster() { echo "$1" | cut -d';' -f1; }

# --- GPU cluster phases (1, 3, 4) ---

# Phase 1: StarDist segmentation (GPU, array 0-11)
PHASE1_JOB=$(strip_cluster "$(sbatch --parsable examples/sbatch_patient_phase1.sh)")
echo "Phase 1 (segmentation): job ${PHASE1_JOB} [gpu]"

# Phase 3: ViT features (GPU, array 0-11, depends on Phase 1)
PHASE3_JOB=$(strip_cluster "$(sbatch --parsable --dependency=afterok:${PHASE1_JOB} examples/sbatch_patient_phase3.sh)")
echo "Phase 3 (features):     job ${PHASE3_JOB} [gpu, depends on ${PHASE1_JOB}]"

# --- HTC cluster phases (2, 5) with delayed start ---

# Phase 2: Module 3 (HTC, array 0-23, delayed 30min for Phase 1)
PHASE2_JOB=$(strip_cluster "$(sbatch --parsable --begin=now+30 examples/sbatch_patient_phase2.sh)")
echo "Phase 2 (module3):      job ${PHASE2_JOB} [htc, delayed 30min]"

# --- Phase 4: GPU, depends on Phase 3, delayed for Phase 2 ---
PHASE4_JOB=$(strip_cluster "$(sbatch --parsable --dependency=afterok:${PHASE3_JOB} --begin=now+120 examples/sbatch_patient_phase4.sh)")
echo "Phase 4 (MIL+Hungarian):job ${PHASE4_JOB} [gpu, depends on ${PHASE3_JOB}, delayed 2h]"

# --- Phase 5: HTC, delayed enough for everything ---
PHASE5_JOB=$(strip_cluster "$(sbatch --parsable --begin=now+240 examples/sbatch_patient_phase5_validate.sh)")
echo "Phase 5 (validation):   job ${PHASE5_JOB} [htc, delayed 4h]"

echo ""
echo "=== All jobs submitted ==="
echo "  Phase 1 (gpu):  ${PHASE1_JOB}"
echo "  Phase 2 (htc):  ${PHASE2_JOB}  [delayed 30min]"
echo "  Phase 3 (gpu):  ${PHASE3_JOB}"
echo "  Phase 4 (gpu):  ${PHASE4_JOB}  [delayed 2h]"
echo "  Phase 5 (htc):  ${PHASE5_JOB}  [delayed 4h]"
echo ""
echo "Monitor with: squeue -u \$(whoami) -M gpu,htc"
