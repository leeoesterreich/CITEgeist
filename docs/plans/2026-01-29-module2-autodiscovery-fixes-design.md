# Module 2 Autodiscovery Fixes Design

Date: 2026-01-29
Branch: `hierarchical_approach`

## Problem

Module 1/2 autodiscovery on the Xenium benchmark produces 14-21 profiles per region, but 70% are single-marker singletons. Three root causes identified through systematic debugging:

1. **Permutation resolution bottleneck**: 199 permutations give minimum p-value = 0.005. BH-FDR at rank 1/351 pairs needs p < 0.00014. Most legitimate pairs can't achieve low enough p-values to survive FDR filtering.

2. **Spatial scale mismatch for stromal markers**: alphaSMA, Vimentin, and CD31 don't colocalize at k=6 neighborhood scale. They may colocalize at tissue-level scales (k=24-48) but multi-scale colocalization is not enabled in the benchmark.

3. **Evaluation uses old 10-type GT**: `evaluate_pipeline_stages.py` scores profiles against 10 granular RNA-based ground truth types that don't include CD4+ T cells. A correctly-discovered CD4+ T cell profile (CD3E, CD4, CD45, CD45RO) gets misattributed to CD8+ T cells.

Additionally, no mechanism exists to distinguish valuable singletons (PanCK covering unique epithelial territory) from noise singletons (PTEN, E-Cadherin).

## Current Performance (5 Xenium regions)

| Metric | Value |
|--------|-------|
| Total profiles discovered | 84 |
| Singletons (1 marker) | 59 (70%) |
| Multi-marker (2+) | 25 (30%) |
| Exact GT matches | 8 (10%) |
| Mean GT coverage | 0.4-0.6 (target: 0.8) |

Multi-marker profiles that do form are often good: Macrophages (perfect in 3/5 regions), B cells (correct in 4/5), T cell core (present in 5/5 but variable composition).

## Design

### Fix 1: Increase Permutation Count

**Goal**: Reduce singleton rate by allowing more pairs to survive FDR.

**Change**: Default `n_permutations` from 199 to 999 in `analyze_marker_colocalization()`.

- Minimum achievable p-value: 0.005 → 0.001
- BH-FDR threshold at rank 1/351: 0.000142. At 999 perms, pairs with very strong colocalization can now reach p=0.001, which passes BH-FDR for ranks >= 7.
- Runtime: ~5x slower for Module 2a (proportional to permutation count). Acceptable on HPC.

**Files**:
- `CITEgeist/model/spatial_colocalization.py`: change default parameter
- `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`: update if passed explicitly
- `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`: update if passed explicitly

### Fix 2: Enable Multi-Scale Colocalization

**Goal**: Allow stromal markers to form edges at tissue-level spatial scales.

**Change**: Pass `multi_scale_k=[6, 12, 24, 48]` and `multi_scale_aggregation="max"` to `analyze_marker_colocalization()` in the benchmark runner and staged evaluation.

The existing multi-scale code computes bivariate Moran's I at each k value, then takes the max across scales per pair. A pair that colocalizes at k=48 but not k=6 still gets a strong score. No changes to the colocalization function internals — just enabling existing parameters.

**Expected effect**: alphaSMA + Vimentin and CD31 + Vimentin edges form, allowing Fibroblast and Endothelial profiles to be discovered.

**Files**:
- `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`: pass multi-scale params
- `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`: pass multi-scale params

### Fix 3: Shared Constants Module + Achievable-7 Evaluation

**Goal**: Single source of truth for achievable-7 definitions; fix evaluation misattribution.

**New file**: `Benchmarking/xenium_benchmarking/benchmark_constants.py`

Contents:
- `ACHIEVABLE_7_CELL_PROFILE_DICT` — Major/Minor marker dict for Module 3
- `ACHIEVABLE_7_MARKER_SIGNATURES` — set-based dict for Jaccard matching
- `GT_TO_ACHIEVABLE_7_MAPPING` — 10→7 collapse mapping
- `CRITICAL_MARKERS` — markers that must be flagged by Module 1
- `ACHIEVABLE_7_GT_MARKERS` — primary/secondary format for profile scoring:

```python
ACHIEVABLE_7_GT_MARKERS = {
    "B cells": {
        "primary": ["CD20"],
        "secondary": ["CD45RA"],
    },
    "CD4+ T cells": {
        "primary": ["CD3E", "CD4"],
        "secondary": ["CD45RO"],
    },
    "CD8+ T cells": {
        "primary": ["CD3E", "CD8A"],
        "secondary": ["GranzymeB"],
    },
    "Macrophages": {
        "primary": ["CD68", "CD163"],
        "secondary": ["CD16"],
    },
    "Endothelial": {
        "primary": ["CD31"],
        "secondary": [],
    },
    "Epithelial": {
        "primary": ["PanCK"],
        "secondary": [],
    },
    "Fibroblasts": {
        "primary": ["alphaSMA", "Vimentin"],
        "secondary": [],
    },
}
```

With this, Profile 0 from Region 0 (CD3E, CD4, CD45, CD45RO) scores:
- CD4+ T cells: `(2 * 2/2 + 1/1) / 3 = 1.0` (exact match)
- CD8+ T cells: `(2 * 1/2 + 0/1) / 3 = 0.333`

**Files**:
- New: `Benchmarking/xenium_benchmarking/benchmark_constants.py`
- `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py`: import from shared module, remove local definitions
- `Benchmarking/xenium_benchmarking/CITEgeist/src/evaluate_pipeline_stages.py`: import from shared module, replace `GT_CELLTYPE_MARKERS`

### Fix 4: Singleton Rescue via Unique Spatial Coverage

**Goal**: Keep singletons that cover unique spatial territory (PanCK for epithelial), drop noise singletons (PTEN, E-Cadherin).

**Why not interest score**: Vimentin has interest score 2.98 (26th of 27) — nearly identical to PTEN (1.80). The multiplicative formula (kurtosis x GMM_SNR x Moran's I) penalizes broadly-expressed markers with flat distributions. Vimentin is spatially important but not peaked.

**Why not reconstruction gain**: PanCK has marginal variance gain of 0.0002 (12th of 15 in greedy selection). Epithelial regions are partially captured by earlier-selected profiles in the variance decomposition.

**Approach**: Unique spatial coverage after Module 2b profile discovery.

Algorithm:
1. Identify all singletons and multi-marker profiles from Module 2b output
2. For each marker, define "on" spots using Module 1 GMM signal component classifications (spots where the marker is in the signal component, not background)
3. Compute union of "on" spots across all multi-marker profiles → "explained territory"
4. For each singleton compute:
   - `unique_coverage`: fraction of its GMM-signal spots outside explained territory
   - `signal_fraction`: fraction of all spots classified as signal by Module 1 GMM (total area)
5. Keep singleton if `unique_coverage >= 0.3` AND `signal_fraction >= 0.05`
6. Drop singletons that fail either criterion

**Data flow change**: Module 1 `MarkerInterestResult` needs to expose per-marker GMM signal spot assignments (binary mask or spot indices). This may require a small addition to the return object — currently it stores aggregate scores but not per-spot classifications.

**Files**:
- `CITEgeist/model/marker_interest.py`: add per-spot GMM signal assignments to `MarkerInterestResult`
- `CITEgeist/model/spatial_colocalization.py`: add singleton rescue logic in or before `select_profiles()`
- Pipeline orchestration code: thread Module 1 results through to Module 2c

## Not Changed

- Module 2b clustering algorithm (hierarchical clustering, dynamic tree cutting, gap-based lineage splitting)
- Module 2a colocalization scoring formula (0.3 Spearman + 0.3 cosine + 0.4 bivariate Moran's I)
- Module 1 marker detection logic
- Module 3 deconvolution
- Core algorithm defaults beyond `n_permutations`

## Validation Plan

Re-run staged evaluation on all 5 Xenium regions after implementation. Success criteria:

1. **Singleton rate drops** from 70% to <40%
2. **Stromal markers group** — alphaSMA + Vimentin in same profile in >=2 regions
3. **Evaluation correctly identifies CD4+ T cells** — profiles containing CD4 match "CD4+ T cells"
4. **No regressions** — Macrophage and B cell profiles remain intact
5. **Noise singletons dropped** — PTEN, E-Cadherin removed; PanCK, Vimentin retained

Expected non-fixes:
- GranzymeB likely remains singleton (biologically sparse, not spatially clustered)
- Region-to-region T cell structure variability may persist (real biology)
- PanCK remains singleton (sole epithelial marker in panel) but is correctly rescued

Run as SLURM array job over 5 regions. Compare before/after on all staged evaluation metrics.
