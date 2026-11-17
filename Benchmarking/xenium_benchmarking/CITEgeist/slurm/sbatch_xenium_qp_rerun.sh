#!/bin/bash
#SBATCH --job-name=xenium_qp_rerun
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
# GPU_RACE_READONLY_OK=1
#
# Xenium benchmark full rerun: 6-type QP, 7-type QP, SACE GEX (regions 0-4)
#
# Run via gpu-race-submit.sh (DO NOT add --partition or --cluster here):
#   gpu-race-submit.sh sbatch_xenium_qp_rerun.sh
#
# Pipeline order:
#   1. run_qp_only_benchmark.py    — 6-type cuOPT QP, all 5 regions
#   2. run_qp_7type_benchmark.py   — 7-type cuOPT QP, all 5 regions
#   3. benchmark_sace_xenium.py    — SACE GEX (own QP props), regions 0-4 sequentially
#   4. Summary printout of canonical_summary JSON files

_cg_find() {
    if [ -n "${CITEGEIST_PATHS_SH:-}" ]; then echo "${CITEGEIST_PATHS_SH}"; return 0; fi
    local p d
    p="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)/../../../_paths.sh"
    if [ -f "$p" ]; then echo "$p"; return 0; fi
    d="${SLURM_SUBMIT_DIR:-}"
    while [ -n "$d" ] && [ "$d" != "/" ] && [ "$d" != "." ]; do
        if [ -f "$d/Benchmarking/_paths.sh" ]; then echo "$d/Benchmarking/_paths.sh"; return 0; fi
        d="$(dirname "$d")"
    done
    echo "cannot find Benchmarking/_paths.sh; set CITEGEIST_PATHS_SH" >&2
    return 1
}
_cg_paths="$(_cg_find)" || exit 1
source "${_cg_paths}"
unset -f _cg_find
unset _cg_paths
set -eo pipefail

eval "$(conda shell.bash hook 2>/dev/null)"
conda activate ${CITEGEIST_ENVS_ROOT:?CITEGEIST_ENVS_ROOT is not set and conda was not found; set it to your conda envs directory}/CITEgeist_env

set -u

export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

REPO="${CITEGEIST_REPO_ROOT}"
SRC="${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/src"

OUTPUT_6TYPE="${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/output_qp_singler_6type"
OUTPUT_7TYPE="${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/output_qp_singler_7type"
OUTPUT_SACE="${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/results_sace_singler_7type"

echo "========================================================"
echo "Xenium benchmark full rerun (6-type QP + 7-type QP + SACE)"
echo "Job ID:   ${SLURM_JOB_ID}"
echo "Node:     $(hostname)"
echo "Start:    $(date)"
echo "GPU:      $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'none')"
echo "Env:      CITEgeist_env"
echo "========================================================"

# ------------------------------------------------------------------
# Stage 1: 6-type QP benchmark (all 5 regions)
# ------------------------------------------------------------------
echo ""
echo "--- Stage 1: 6-type QP benchmark ---"
echo "Output: ${OUTPUT_6TYPE}"
echo "Start:  $(date)"

python "${SRC}/run_qp_only_benchmark.py" \
    --output-dir "${OUTPUT_6TYPE}"

echo "Stage 1 done: $(date)"

# ------------------------------------------------------------------
# Stage 2: 7-type QP benchmark (all 5 regions)
# ------------------------------------------------------------------
echo ""
echo "--- Stage 2: 7-type QP benchmark ---"
echo "Output: ${OUTPUT_7TYPE}"
echo "Start:  $(date)"

python "${SRC}/run_qp_7type_benchmark.py" \
    --output-dir "${OUTPUT_7TYPE}"

echo "Stage 2 done: $(date)"

# ------------------------------------------------------------------
# Stage 3: SACE GEX benchmark (regions 0-4, sequential)
# ------------------------------------------------------------------
echo ""
echo "--- Stage 3: SACE GEX benchmark (regions 0-4) ---"
echo "Output: ${OUTPUT_SACE}"
echo "Start:  $(date)"

for REGION in 0 1 2 3 4; do
    echo ""
    echo "  [SACE] Region ${REGION} start: $(date)"
    python "${SRC}/benchmark_sace_xenium.py" \
        --region "${REGION}" \
        --output_dir "${OUTPUT_SACE}"
    echo "  [SACE] Region ${REGION} done:  $(date)"
done

echo "Stage 3 done: $(date)"

# ------------------------------------------------------------------
# Stage 4: Summary
# ------------------------------------------------------------------
echo ""
echo "========================================================"
echo "RESULTS SUMMARY"
echo "========================================================"

python - <<'PYEOF'
import json
from pathlib import Path

REPO = Path(__import__("os").environ["CITEGEIST_REPO_ROOT"])
BM   = REPO / "Benchmarking" / "xenium_benchmarking" / "CITEgeist"

summaries = [
    ("6-type QP",  BM / "output_qp_singler_6type" / "canonical_summary_6type.json"),
    ("7-type QP",  BM / "output_qp_singler_7type" / "canonical_summary_7type.json"),
]

for label, path in summaries:
    print(f"\n--- {label} ---")
    if not path.exists():
        print(f"  [MISSING] {path}")
        continue
    with open(path) as f:
        s = json.load(f)
    mean_r = s.get("overall_mean_r", float("nan"))
    std_r  = s.get("overall_std_r",  float("nan"))
    print(f"  Regions:        {s.get('n_regions', '?')}/5")
    print(f"  Overall mean r: {mean_r:.4f} ± {std_r:.4f}")
    per_region = s.get("per_region_r", {})
    for rid in sorted(per_region.keys(), key=lambda x: int(x)):
        r = per_region[rid]
        print(f"    Region {rid}: {r:.4f}" if r is not None else f"    Region {rid}: N/A")
    per_type = s.get("per_type_mean_r", {})
    if per_type:
        print("  Per-type mean r:")
        for ct, r in sorted(per_type.items()):
            print(f"    {ct:<20s}: {r:.4f}")

# SACE: check layer outputs exist per region
print("\n--- SACE GEX ---")
sace_dir = BM / "results_sace_singler_7type"
TYPE_FILES = [
    "B_cells_layer.csv", "CD4pos_T_cells_layer.csv", "CD8pos_T_cells_layer.csv",
    "Macrophages_layer.csv", "Endothelial_layer.csv", "Epithelial_layer.csv",
    "Fibroblasts_layer.csv",
]
for rid in range(5):
    region_layers = sace_dir / f"Xenium_region_{rid}" / "layers"
    found = [f for f in TYPE_FILES if (region_layers / f).exists()]
    missing = [f for f in TYPE_FILES if f not in found]
    status = "OK" if not missing else f"MISSING {missing}"
    print(f"  Region {rid}: {len(found)}/7 layer files — {status}")
PYEOF

echo ""
echo "========================================================"
echo "All stages complete: $(date)"
echo "========================================================"
