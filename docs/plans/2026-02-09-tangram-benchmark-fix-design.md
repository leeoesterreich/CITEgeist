# Tangram Xenium Benchmark Fix Design

**Date:** 2026-02-09
**Status:** Ready for Implementation
**Scope:** Fix two bugs in Xenium Tangram benchmark to match simulated benchmark methodology

---

## Problem Statement

The Xenium Tangram benchmark has two critical bugs that make results incomparable to both the simulated benchmark and other methods:

1. **Preprocessing mismatch**: Reference data is log-transformed but spatial data is raw counts
2. **Output normalization bug**: Cell type proportions don't sum to 1

These issues explain why Tangram shows near-uniform predictions (~36% Fibroblasts everywhere) with collapsed variance (std=0.03 vs truth std=0.09).

---

## Issue 1: Preprocessing Mismatch

### Current State (Wrong)

| Data | Processing | Space |
|------|------------|-------|
| Reference | `normalize_total()` + `log1p()` + `scale()` | Log-space |
| Spatial | Raw summed counts | Linear-space |

Tangram tries to learn a mapping between incompatible feature spaces.

### Simulated Benchmark (Correct)

```python
# Reference (lines 50-51 of Tangram.py)
sc.pp.normalize_total(wu_adata)
sc.pp.log1p(wu_adata)

# Spatial (line 62 of Tangram.py)
sc.pp.log1p(ad_sp)
```

Both datasets in same log-space before `tg.pp_adatas()`.

### Fix

Add normalization + log-transformation to spatial data in `run_benchmark.py`:

```python
# After loading spatial data (around line 347)
adata_sp = sc.read_h5ad(gex_path)
logger.info(f"  Spatial shape: {adata_sp.shape}")

# ADD: Match preprocessing to reference
sc.pp.normalize_total(adata_sp)
sc.pp.log1p(adata_sp)
logger.info("  Applied normalize_total + log1p to spatial data")
```

---

## Issue 2: Output Normalization Bug

### Current State (Wrong)

```python
# run_benchmark.py lines 169-175
if 'tangram_ct_pred' in adata_sp.obsm:
    proportions = pd.DataFrame(
        adata_sp.obsm['tangram_ct_pred'],  # Raw weights from project_cell_annotations
        ...
    )
# No normalization - proportions sum to ~39 on average
```

### Simulated Benchmark (Correct)

```python
# Tangram.py line 105
proportions = tangram_ct_pred.div(tangram_ct_pred.sum(axis=1), axis=0)
proportions = proportions.fillna(0)
```

Explicit row-normalization to sum to 1.

### Fix

Add normalization after extracting proportions in `extract_proportions()`:

```python
# After creating proportions DataFrame (around line 175)
proportions = pd.DataFrame(
    adata_sp.obsm['tangram_ct_pred'],
    index=adata_sp.obs_names,
    columns=adata_sp.uns['tangram_ct_pred_names'] if 'tangram_ct_pred_names' in adata_sp.uns else None
)

# ADD: Normalize to sum to 1 (matching simulated benchmark)
proportions = proportions.div(proportions.sum(axis=1), axis=0)
proportions = proportions.fillna(0)
logger.info(f"  Normalized proportions (sum={proportions.sum(axis=1).mean():.3f})")
```

---

## Files to Modify

1. `Benchmarking/xenium_benchmarking/Tangram/src/run_benchmark.py`
   - Add spatial preprocessing after line 347
   - Add proportion normalization in `extract_proportions()` around line 175

---

## Validation

After fix, verify:

1. **Preprocessing**: Both reference and spatial in log-space before `tg.pp_adatas()`
2. **Proportions**: All rows sum to 1.0 (tolerance 1e-6)
3. **Variance**: Predictions should have realistic variance (closer to ground truth std)

Run on one region first:
```bash
python run_benchmark.py --region-id 0 --ref-path /path/to/reference.h5ad --output-dir test_output
```

Check output:
```python
import pandas as pd
props = pd.read_csv("test_output/Xenium_region_0/Xenium_region_0_cell_prop_finetuned_results.csv", index_col=0)
print(f"Row sums: mean={props.sum(axis=1).mean():.4f}, std={props.sum(axis=1).std():.6f}")
print(f"Fibroblast std: {props['Fibroblasts'].std():.4f}")  # Should be > 0.03
```

---

## Re-run Plan

After fixing, re-run Tangram on all 14 regions:

```bash
# Submit batch job for all regions
sbatch sbatch_run_all_regions.sh
```

Then re-run evaluation to update comparison results.

---

## Risk Assessment

**Low risk**: These are clear bugs where the Xenium implementation doesn't match the working simulated benchmark. The fixes align with Tangram's documented usage.

**Potential concern**: If the reference is already heavily scaled (z-scored), the spatial data should match. Check that `process_reference.py` export for Tangram uses log-space, not scaled space.

Looking at `process_reference.py` line 367-368:
```python
# Tangram needs normalized counts + cell type labels
adata_export = adata.copy()  # This has log-transformed .X (from line 215)
```

The reference export uses log-transformed data (post `log1p()`, pre `scale()`), so log-transforming spatial data is correct.
