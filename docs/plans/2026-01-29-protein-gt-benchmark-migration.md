# Protein GT Benchmark Migration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Switch entire Xenium benchmarking framework from RNA-based ground truth (10 granular types with problematic "Mixed Immune" → CD4+ T cells mapping) to protein-gated ground truth with correct achievable-7 cell types, and re-run all 5 methods (CITEgeist, Cell2Location, Tangram, RCTD, Seurat).

**Architecture:** The protein GT uses hierarchical single-cell protein gating (like flow cytometry) BEFORE pseudo-Visium aggregation, avoiding circular logic. Currently outputs 8 types (separate Myofibroblasts and Stromal) — needs collapsing to 7 (Fibroblasts = Myofibroblasts + Stromal). All competitor methods use scRNA-seq reference with 6 types (combined "T cells") — needs reprocessing to 7 types matching the protein GT (split CD4+ T / CD8+ T). The evaluation script needs updating to work with protein GT directly (no more 10→7 collapse mapping).

**Tech Stack:** Python 3.10, scanpy, anndata, Gurobi, Cell2Location, Tangram, RCTD (R), Seurat (R), SLURM HPC

---

## Task 1: Fix protein GT generation to output achievable-7 types

**Files:**
- Modify: `Benchmarking/xenium_pseudovisium/src/create_protein_gt.py:30-39` (CELL_TYPE_ORDER)
- Modify: `Benchmarking/xenium_pseudovisium/src/create_protein_gt.py:132-142` (Gates 7-8 labels)

**Step 1: Update CELL_TYPE_ORDER**

Change lines 30-39 from:

```python
CELL_TYPE_ORDER = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Myofibroblasts",
    "Stromal",
]
```

To:

```python
CELL_TYPE_ORDER = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]
```

**Step 2: Update Gate 7 (Myofibroblasts → Fibroblasts)**

Change line 135 from:
```python
    cell_types[myofib] = 'Myofibroblasts'
```
To:
```python
    cell_types[myofib] = 'Fibroblasts'
```

Update log message on line 136 to say `"7. Fibroblasts (alphaSMA high)"`.

**Step 3: Update Gate 8 (Stromal → also Fibroblasts)**

Change line 141 from:
```python
    cell_types[stromal] = 'Stromal'
```
To:
```python
    cell_types[stromal] = 'Fibroblasts'
```

Update log message on line 142 to say `"8. Fibroblasts/Stromal (Vimentin+ remaining)"`.

Note: Both gates now label as "Fibroblasts", so the output will have 7 types total. `calculate_spot_proportions()` uses `CELL_TYPE_ORDER` (line 166), which now has 7 types — the counts for both alphaSMA-high and Vimentin-positive cells will be summed under "Fibroblasts" automatically since they share the same label.

**Step 4: Commit**

```bash
git add Benchmarking/xenium_pseudovisium/src/create_protein_gt.py
git commit -m "fix: collapse Myofibroblasts+Stromal to Fibroblasts in protein GT generation"
```

---

## Task 2: Re-generate protein GT with achievable-7 types

**Files:**
- Create SLURM script: `Benchmarking/xenium_pseudovisium/slurm/regenerate_protein_gt.sh`
- Output: `Benchmarking/xenium_pseudovisium/data_protein_gt/` (overwritten)

**Step 1: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=regen_protein_gt
#SBATCH --output=Benchmarking/xenium_pseudovisium/slurm/slurm_log/regen_protein_gt_%j.out
#SBATCH --error=Benchmarking/xenium_pseudovisium/slurm/slurm_log/regen_protein_gt_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"

echo "Re-generating protein GT with achievable-7 types..."
python "${REPO}/Benchmarking/xenium_pseudovisium/src/create_protein_gt.py"

echo "Done! Verify output:"
python -c "
import pandas as pd
gt = pd.read_csv('${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt/ground_truth/Xenium_region_0_prop.csv', index_col=0)
print('Columns:', list(gt.columns))
print('Shape:', gt.shape)
for c in [col for col in gt.columns if col not in ['n_cells','spot_x','spot_y']]:
    print(f'  {c}: mean={gt[c].mean():.4f}')
"
```

**Step 2: Submit job**

```bash
mkdir -p Benchmarking/xenium_pseudovisium/slurm/slurm_log
sbatch Benchmarking/xenium_pseudovisium/slurm/regenerate_protein_gt.sh
```

**Step 3: Verify output has 7 cell types**

Expected columns: B cells, CD4+ T cells, CD8+ T cells, Macrophages, Endothelial, Epithelial, Fibroblasts, n_cells, spot_x, spot_y

---

## Task 3: Archive RNA-based GT directories

**Step 1: Create archive and move**

```bash
cd Benchmarking/xenium_pseudovisium
mkdir -p _archive
mv data_granular_gt _archive/data_granular_gt
mv data_rna_gt _archive/data_rna_gt
```

**Step 2: Verify protein GT is the only active GT**

```bash
ls -la Benchmarking/xenium_pseudovisium/
# Should show data_protein_gt/ as the only data_*_gt directory
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_pseudovisium/_archive/
git commit -m "chore: archive RNA-based GT directories (replaced by protein GT)"
```

Note: The `_archive/` directories are large (h5ad files). Consider adding to `.gitignore` if not already tracked. The important thing is that no scripts point to them anymore.

---

## Task 4: Update evaluate_benchmark.py for protein GT

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_benchmark.py:26-76`

**Step 1: Replace GT collapse logic with protein GT passthrough**

The protein GT already has the correct 7 types. Replace the `GT_TO_ACHIEVABLE_7_MAPPING` dict (lines 31-42) and `collapse_gt_to_achievable_7` function (lines 55-76) with a simpler version:

```python
# =============================================================================
# ACHIEVABLE-7 GROUND TRUTH (Protein-Gated)
# =============================================================================
# Protein GT already outputs 7 cell types directly.
# No collapse mapping needed — types match achievable-7 exactly.

ACHIEVABLE_7_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Endothelial",
    "Epithelial",
    "Fibroblasts",
]


def filter_gt_to_achievable_7(gt_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter protein GT to achievable-7 cell type columns.

    Protein GT already has the correct 7 types. This function simply
    selects those columns and drops metadata columns (n_cells, spot_x, etc.).
    """
    present_types = [ct for ct in ACHIEVABLE_7_CELL_TYPES if ct in gt_df.columns]
    return gt_df[present_types].copy()
```

**Step 2: Update evaluate_region() to use filter instead of collapse**

Change line 186 from:
```python
        gt_df = collapse_gt_to_achievable_7(gt_df)
```
To:
```python
        gt_df = filter_gt_to_achievable_7(gt_df)
```

**Step 3: Update module docstring**

Replace lines 1-11 (the module docstring) — remove reference to "10 granular types collapsed to 7" and RNA-based ground truth. Update to describe protein-gated GT.

**Step 4: Commit**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/src/evaluate_benchmark.py
git commit -m "feat: update evaluation to use protein GT directly (no collapse needed)"
```

---

## Task 5: Update benchmark_constants.py GT mapping

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/benchmark_constants.py:98-110`

**Step 1: Replace GT_TO_ACHIEVABLE_7_MAPPING**

The old mapping collapsed 10 RNA types to 7. With protein GT, the 7 types map 1:1. Replace lines 98-110:

```python
# Protein GT → Achievable-7 mapping (identity — protein GT already has correct types)
GT_TO_ACHIEVABLE_7_MAPPING: Dict[str, str] = {
    "B cells": "B cells",
    "CD4+ T cells": "CD4+ T cells",
    "CD8+ T cells": "CD8+ T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/benchmark_constants.py
git commit -m "feat: update GT mapping for protein GT (1:1 identity, no Mixed Immune)"
```

---

## Task 6: Update CITEgeist run_benchmark.py default input dir

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py:543`

**Step 1: Update default --input-dir**

Change line 543 from:
```python
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_granular_gt"),
```
To:
```python
        default=str(REPO_ROOT / "Benchmarking/xenium_pseudovisium/data_protein_gt"),
```

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py
git commit -m "feat: default CITEgeist benchmark to protein GT input"
```

---

## Task 7: Update all competitor SLURM scripts to use protein GT

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/Cell2Location/slurm/xenium_benchmark.sh:31`
- Modify: `Benchmarking/xenium_benchmarking/Tangram/slurm/xenium_benchmark.sh:27`
- Modify: `Benchmarking/xenium_benchmarking/RCTD/slurm/xenium_benchmark.sh:27`
- Modify: `Benchmarking/xenium_benchmarking/Seurat/slurm/xenium_benchmark.sh:27`

**Step 1: Update INPUT_DIR in all 4 scripts**

In each script, change:
```bash
INPUT_DIR="...data_granular_gt"
```
To:
```bash
INPUT_DIR="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/xenium_pseudovisium/data_protein_gt"
```

Also update OUTPUT_DIR to new output directories (e.g., `output_protein_gt` instead of `output_granular`) so we don't overwrite existing results.

And update REF_DIR paths to point to `processed_protein7/` (created in Task 8) instead of `processed_granular/`:
- Cell2Location: `REF_DIR="${C2L_DIR}/reference_model_protein7"` (line 32)
- Tangram: `REF_PATH="${BASE_DIR}/reference_data/GSE156632/processed_protein7/tangram/reference.h5ad"` (line 28)
- RCTD: `REF_DIR="${BASE_DIR}/reference_data/GSE156632/processed_protein7/rctd"` (line 28)
- Seurat: `REF_DIR="${BASE_DIR}/reference_data/GSE156632/processed_protein7/seurat"` (line 28)

**Step 2: Commit**

```bash
git add Benchmarking/xenium_benchmarking/*/slurm/xenium_benchmark.sh
git commit -m "feat: update all competitor SLURM scripts for protein GT"
```

---

## Task 8: Reprocess scRNA-seq reference for 7 achievable cell types

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/reference_data/GSE156632/src/process_reference.py:55-62`
- Create SLURM script: `Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/reprocess_protein7.sh`
- Output: `Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed_protein7/`

**Step 1: Update MARKER_GENES for 7 types**

The current `process_reference.py` uses 6 types with combined "T cells". We need 7 types matching protein GT. Replace lines 55-62:

```python
MARKER_GENES = {
    "B cells": ["CD79A", "CD79B", "MS4A1", "CD19", "PAX5", "BANK1"],
    "CD4+ T cells": ["CD3D", "CD3E", "CD3G", "CD4", "TRAC", "IL7R", "TCF7"],
    "CD8+ T cells": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7"],
    "Macrophages": ["CD68", "CD163", "CSF1R", "MARCO", "MSR1", "CD14", "FCGR3A"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "VIM", "ACTA2", "FAP", "PDGFRA"],
    "Epithelial": ["EPCAM", "KRT18", "KRT8", "CA9", "PAX8", "CDH1", "KRT19"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CD34", "KDR", "FLT1"],
}
```

Key changes:
- "T cells" split into "CD4+ T cells" (CD4, IL7R, TCF7 added) and "CD8+ T cells" (CD8A, CD8B, GZMA, GZMB, PRF1, NKG7 added)
- Names match protein GT exactly

**Step 2: Update TARGET_CELL_TYPES in all competitor benchmark scripts**

In `Cell2Location/src/run_benchmark.py` (lines 35-42) and `Tangram/src/run_benchmark.py` (lines 36-43), change:

```python
TARGET_CELL_TYPES = [
    "B cells",
    "CD4+ T cells",
    "CD8+ T cells",
    "Macrophages",
    "Fibroblasts",
    "Epithelial",
    "Endothelial",
]
```

For RCTD and Seurat (R scripts), they read cell types from the reference CSV, so they'll automatically pick up the new 7 types.

**Step 3: Write SLURM script for reprocessing**

```bash
#!/bin/bash
#SBATCH --job-name=reprocess_ref
#SBATCH --output=Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/slurm_log/reprocess_%j.out
#SBATCH --error=Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/slurm_log/reprocess_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
INPUT_DIR="/ix1/alee/LO_LAB/General/Public_Data/GSE156632"
OUTPUT_DIR="${REPO}/Benchmarking/xenium_benchmarking/reference_data/GSE156632/processed_protein7"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

echo "Reprocessing reference for 7 protein-GT-matched cell types..."
python "${REPO}/Benchmarking/xenium_benchmarking/reference_data/GSE156632/src/process_reference.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}"

echo "Done!"
```

**Step 4: Submit and verify**

```bash
sbatch Benchmarking/xenium_benchmarking/reference_data/GSE156632/slurm/reprocess_protein7.sh
```

Verify output has 7 cell types in `processed_protein7/{cell2location,tangram,rctd,seurat}/reference.h5ad`.

**Step 5: Commit**

```bash
git add Benchmarking/xenium_benchmarking/reference_data/GSE156632/src/process_reference.py
git add Benchmarking/xenium_benchmarking/Cell2Location/src/run_benchmark.py
git add Benchmarking/xenium_benchmarking/Tangram/src/run_benchmark.py
git commit -m "feat: reprocess scRNA-seq reference with 7 cell types matching protein GT"
```

---

## Task 9: Train Cell2Location reference model on new 7-type reference

**Files:**
- Existing: `Benchmarking/xenium_benchmarking/Cell2Location/slurm/train_reference.sh` (if exists)
- Output: `Benchmarking/xenium_benchmarking/Cell2Location/reference_model_protein7/`

Cell2Location requires a separate reference model training step before deconvolution. This needs to be re-run with the new 7-type reference.

**Step 1: Create/update training SLURM script**

Point to `processed_protein7/cell2location/reference.h5ad` as input.

**Step 2: Submit and verify training completes**

**Step 3: Update Cell2Location benchmark SLURM script REF_DIR to use new model**

This should already be done in Task 7, but verify after training.

---

## Task 10: Run CITEgeist benchmark against protein GT

**Files:**
- Existing SLURM script or create new one
- Output: `Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt/`

**Step 1: Create SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=cg_protein_gt
#SBATCH --output=Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/protein_gt_%A_%a.out
#SBATCH --error=Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log/protein_gt_%A_%a.err
#SBATCH --array=0-4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

set -e
module load gurobi/12.0.3

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

REGION_ID=${SLURM_ARRAY_TASK_ID}

python "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py" \
    --region-id ${REGION_ID} \
    --mode manual \
    --input-dir "${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt" \
    --output-dir "${REPO}/Benchmarking/xenium_benchmarking/CITEgeist/output_protein_gt"
```

**Step 2: Submit**

```bash
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_protein_gt.sh
```

---

## Task 11: Run all competitor methods against protein GT

**Prerequisite:** Tasks 8-9 must complete first (reference reprocessing + C2L training).

**Step 1: Submit all 4 competitor benchmarks**

```bash
sbatch Benchmarking/xenium_benchmarking/Cell2Location/slurm/xenium_benchmark.sh
sbatch Benchmarking/xenium_benchmarking/Tangram/slurm/xenium_benchmark.sh
sbatch Benchmarking/xenium_benchmarking/RCTD/slurm/xenium_benchmark.sh
sbatch Benchmarking/xenium_benchmarking/Seurat/slurm/xenium_benchmark.sh
```

**Step 2: Monitor and verify all jobs complete**

---

## Task 12: Evaluate all methods and compare

**Step 1: Run evaluation for each method**

```bash
REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
GT_DIR="${REPO}/Benchmarking/xenium_pseudovisium/data_protein_gt"
EVAL_SCRIPT="${REPO}/Benchmarking/xenium_benchmarking/evaluation/src/evaluate_benchmark.py"

for method in CITEgeist Cell2Location Tangram RCTD Seurat; do
    PRED_DIR="${REPO}/Benchmarking/xenium_benchmarking/${method}/output_protein_gt"
    python "${EVAL_SCRIPT}" \
        --gt-dir "${GT_DIR}" \
        --pred-dir "${PRED_DIR}" \
        --n-regions 5 \
        --output "${REPO}/Benchmarking/xenium_benchmarking/evaluation/results_protein_gt/${method}.json"
done
```

**Step 2: Compare results across methods**

Check per-cell-type Pearson r, JSD, RMSE, MAE. In particular verify:
- CD4+ T cells now have reasonable correlation (should be much better than the -0.05 we saw with RNA GT)
- Fibroblasts maintain the improvement from alphaSMA-only profile

**Step 3: Commit final results**

```bash
git add Benchmarking/xenium_benchmarking/evaluation/results_protein_gt/
git commit -m "feat: benchmark results using protein-gated ground truth"
```

---

## Dependency Graph

```
Task 1 (fix GT script) → Task 2 (re-generate GT)
Task 2 → Task 3 (archive old GT)
Task 2 → Task 4 (update evaluate_benchmark.py)
       → Task 5 (update benchmark_constants.py)
       → Task 6 (update CITEgeist input dir)
       → Task 7 (update competitor SLURM scripts)
Task 2 → Task 10 (run CITEgeist)

Task 8 (reprocess reference) → Task 9 (train C2L model)
Task 7 + Task 8 + Task 9 → Task 11 (run competitors)

Task 10 + Task 11 → Task 12 (evaluate all)
```

Tasks 1-7 are code changes (can be done sequentially fast).
Tasks 8-9 are SLURM jobs (blocking on compute).
Tasks 10-11 are SLURM jobs (blocking on Task 2 + Task 9).
Task 12 is evaluation (blocking on all benchmark jobs completing).
