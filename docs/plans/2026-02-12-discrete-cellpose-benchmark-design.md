# Discrete Cell Assignment + Cellpose Benchmarking Design

**Date:** 2026-02-12
**Branch:** `feat/discrete-cell-assignment` → merge into `dev`
**Status:** Approved

## Overview

Merge the discrete cell assignment feature (IQP solver with nuclei counts) with recent Cellpose fixes on `dev`, then run full benchmarking on Xenium pseudo-Visium regions using the complete pipeline: Cellpose segmentation → discrete assignment → GEX deconvolution → evaluation.

## Phase 1: Merge Strategy

### Branches Involved

- **`feat/discrete-cell-assignment`**: 4 commits implementing IQP discrete cell assignment
  - `973f10f1` docs: add discrete cell assignment design for Module 3
  - `55c3a154` docs: add discrete cell assignment implementation plan
  - `b24ddd6f` feat(model): implement discrete cell assignment with IQP solver
  - `7b43982e` docs: add discrete cell assignment example and documentation

- **`dev`**: 4 commits since branch point
  - `e3273bc2` feat(benchmark): add CARD as benchmarking method
  - `20d21b60` results(benchmark): update Tangram results after preprocessing fix
  - `95d53196` fix(benchmark): exclude VIM from CARD benchmark for fair comparison
  - `049700f2` Merge branch 'feat/cellpose-nuclei-prior' into dev

### Merge Approach

1. In `feat-discrete-cell-assignment` worktree, merge `dev`:
   ```bash
   cd .worktrees/feat-discrete-cell-assignment
   git merge dev
   ```

2. Resolve `segmentation.py` conflict:
   - Both branches modified this file
   - Changes are additive (different functions)
   - Keep Cellpose auto-detection from dev + discrete assignment helpers from feature branch

3. Verify with unit tests:
   ```bash
   module load gurobi/12.0.3
   pytest tests/test_discrete_assignment.py -v
   ```

4. Merge back into `dev`:
   ```bash
   cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
   git merge feat/discrete-cell-assignment
   ```

## Phase 2: Benchmarking Pipeline

### Pipeline Flow

```
Morphology image → Cellpose → nuclei_counts per spot
                                    ↓
Protein data + nuclei_counts → IQP discrete assignment → cell_counts
                                    ↓
GEX data + cell_counts → Phase 2 deconvolution → E[spot, celltype, gene]
                                    ↓
Evaluate proportions (vs GT) + Evaluate GEX (vs GT marker expression)
```

### Configuration

- **Benchmark config**: `achievable_7` (7 cell types with reliable markers)
- **Regions**: 5 Xenium pseudo-Visium regions (Xenium_region_0 through Xenium_region_4)
- **Image source**: `Benchmarking/xenium_benchmarking/scResolve/images/morphology_hires/Xenium_region_X/morphology.png`

### New Benchmark Script

**Location**: `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py`

**Functionality**:
1. Load morphology crop for region
2. Run Cellpose segmentation → nuclei counts per spot
3. Load CITE and GEX h5ad for region
4. Run `preprocess_antibody_discrete()` (CLR without per-spot normalization)
5. Run `run_discrete_cell_assignment()` with nuclei counts
6. Run `run_cell_expression_pass1(use_discrete_mode=True)`
7. Save predictions to output directory

### SLURM Resources

- **Cellpose**: GPU partition, ~5 min/region
- **Discrete + GEX**: CPU partition, ~30 min/region

### Output Location

`Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_cellpose/`

## Phase 3: Evaluation

### Metrics

**Proportion metrics**:
- Pearson correlation
- RMSE
- JSD (Jensen-Shannon Divergence)
- Per cell type breakdown

**GEX metrics**:
- Pearson correlation
- RMSE
- NRMSE (Normalized RMSE) - **NEW, to be added**
- MAE

### Code Changes

Update `compare_all_methods_gex.py` to add NRMSE:
```python
nrmse = rmse / (gt_flat.max() - gt_flat.min())  # or use std normalization
```

## Phase 4: Results Tracking

### Backup Existing Results

Before running new benchmarks:
```bash
cp evaluation/results/method_comparison/full_comparison_gex.json \
   evaluation/results/method_comparison/full_comparison_gex_continuous_baseline.json
```

### Git Track Comparison Files

Add to version control:
```bash
git add evaluation/results/method_comparison/*.json
git commit -m "results(benchmark): track comparison JSONs for metric history"
```

Commit after each benchmark run with descriptive message to track metric changes over time.

## Deliverables

1. **Merged `dev` branch** with discrete cell assignment + Cellpose fixes
2. **Benchmark results** for CITEgeist-discrete on achievable_7 (5 regions)
3. **Updated comparison JSONs** with proportion + GEX metrics (including NRMSE)
4. **Baseline backup** of continuous CITEgeist results for comparison

## Files Modified/Created

| File | Change |
|------|--------|
| `CITEgeist/model/segmentation.py` | MERGE: Combine Cellpose fixes + discrete helpers |
| `Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_cellpose.py` | NEW |
| `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_cellpose.sh` | NEW |
| `Benchmarking/xenium_benchmarking/evaluation/src/compare_all_methods_gex.py` | ADD: NRMSE metric |
| `evaluation/results/method_comparison/full_comparison_gex.json` | UPDATE + git track |
| `evaluation/results/method_comparison/full_comparison.json` | UPDATE + git track |

## Success Criteria

1. Clean merge with no regressions (unit tests pass)
2. Cellpose successfully segments all 5 regions
3. Discrete assignment completes for all regions
4. Metrics computed and comparable to other methods
5. Results JSON committed with full history
