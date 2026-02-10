# Gene-Parallel U-Net Spatial Smoothing for Module 3 Pass 2

**Date:** 2026-02-10
**Status:** Design Complete
**Goal:** Improve gene expression deconvolution accuracy by leveraging spatial neighborhood structure

---

## Problem Statement

Current Module 3 Pass 2 optimizes each spot independently. Gene expression deconvolution has poor correlation (Pearson r = 0.05–0.35 for some cell types) because low-abundance cell types lack sufficient signal to reliably assign gene counts.

**Current benchmark performance (Xenium, Pearson r):**
- Macrophages: 0.52–0.62 (best)
- CD8+ T cells: 0.42–0.48
- Endothelial/Epithelial: 0.33–0.49
- B cells, CD4+ T, Fibroblasts: 0.05–0.35 (worst)

---

## Solution Overview

Add proportion-weighted spatial Laplacian regularization via a hierarchical U-Net architecture. High-confidence spots (high cell type proportion) anchor the solution; low-confidence neighbors get pulled toward them. The hierarchy discovers regional expression patterns while preserving biological heterogeneity (e.g., cancer subclones).

**Key insight:** A spot with 5% macrophages can't reliably estimate macrophage expression from its own counts. But if neighboring spots have 30-40% macrophages, their expression estimates can propagate to the low-proportion spot, improving the overall solution.

---

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Information source | Spatial neighborhood structure | Nearby spots should have correlated expression |
| Smoothing mechanism | Laplacian regularization | Penalizes expression differences between neighbors |
| Confidence metric | Cell type proportion | High-proportion spots more reliably estimate that cell type's expression |
| Integration | Joint optimization | Single objective with smoothing built in |
| Pull direction | Asymmetric (high → low) | High-proportion anchors, low-proportion gets pulled |
| Heterogeneity handling | U-Net hierarchy | Discovers regional patterns without over-smoothing |
| λ_spatial tuning | Cross-validation | Held-out spots, optimize Pearson r |
| Computational strategy | Gene-parallel decomposition | Genes don't couple; 9,800 vars per gene is tractable |

---

## Architecture

### Gene-Parallel Decomposition

Genes are independent in the objective (no cross-gene coupling). Each gene solved as a separate QP:

```
For each gene g in parallel (15,000 genes):

    Variables: X[spot, celltype] for all spots
               1400 spots × 7 cell types = 9,800 variables

    Solve U-Net hierarchy for gene g

    Output: X[spot, celltype] for gene g
```

### U-Net Hierarchy (Per Gene)

```
UPWARD PASS (encode regional patterns):

  Level 0: Raw observations per spot
      ↓ solve local QP (radius 1-2, ~20 spots)
  Level 1: Local estimates + confidence
      ↓ pool (proportion-weighted spatial regions)
  Level 2: Regional estimates (~100 spots)
      ↓ pool
  Level 3: Global cell-type profile for gene g

DOWNWARD PASS (refine with regional context):

  Level 3: Global profile
      ↓ broadcast as weak prior
  Level 2: Re-solve regional QP (regional prior)
      ↓ broadcast as prior
  Level 1: Re-solve local QP (local prior from region)
      ↓
  Level 0: Final spot estimates for gene g
```

**Why U-Net (not simple Laplacian):**
- Simple Laplacian would over-smooth regional heterogeneity
- Example: ESR1+ vs ESR1- cancer subclones would be averaged out
- U-Net learns regional patterns at coarse levels
- Downward pass refines within regions without forcing cross-region homogeneity

---

## Mathematical Formulation

### Upward Pass (Level L, Gene g)

```
max  Σ_i∈batch Σ_t [ enrichment[i,t,g] × prop[i,t] × X[i,t] ]
   - λ_reg × Σ_i Σ_t X[i,t]²
   - λ_spatial × Σ_(i~j∈batch) Σ_t [ w_ij[t] × (X[j,t] - X[i,t])² ]

subject to:
   Σ_t X[i,t] = observed[i,g]   ∀ spots i in batch  (count conservation)
   X[i,t] ≥ 0
```

### Asymmetric Proportion-Weighted Laplacian

```
w_ij[t] = prop[i,t] / (prop[i,t] + prop[j,t] + ε)
```

Effect for neighboring spots i (40% cell type t) and j (5% cell type t):
- Weight ≈ 0.89: strong pull from i to j
- Spot i anchors; spot j gets pulled toward i's expression

### Downward Pass (Level L, Gene g)

Adds prior term pulling toward coarser level estimate:

```
max  [ same terms as upward ]
   - λ_prior × Σ_i Σ_t [ (1 - prop[i,t]) × (X[i,t] - X_prior[i,t,g])² ]
```

Where:
- `X_prior[i,t,g]` = estimate from Level L+1 for cell type t, gene g
- `(1 - prop[i,t])` = low-confidence spots shrink more toward prior

### Pooling Between Levels (Upward)

```python
# For each region R at level L+1, for each cell type t:
X_region[R,t,g] = Σ_{i∈R} prop[i,t] × X[i,t,g] / Σ_{i∈R} prop[i,t]
```

### Broadcasting Prior (Downward)

```python
# Each spot i in region R gets the regional estimate as prior:
X_prior[i,t,g] = X_region[R,t,g]
```

---

## Cross-Validation

### Parameters to Tune

| Parameter | Role | Search range |
|-----------|------|--------------|
| `λ_spatial` | Laplacian smoothing strength | [0.01, 0.05, 0.1, 0.25, 0.5] |
| `λ_prior` | Downward pass prior strength | [0.1, 0.25, 0.5, 1.0] |
| `n_levels` | Hierarchy depth | [2, 3, 4] |

### Hold-Out Strategy

1. Hold out 15% of spots randomly (stratified by spatial region)
2. For held-out spots: mask gene expression, keep proportions
3. Run full U-Net pipeline on training spots
4. Predict held-out spots via:
   - Regional prior (from upward pass pooling)
   - Interpolation from trained neighbors
5. Evaluate: Pearson r, RMSE vs held-out ground truth

### Selection Criterion

- Primary: mean Pearson r across cell types (weighted by proportion)
- Secondary: RMSE

---

## Computational Considerations

### Parallelization Strategy

```
Outer parallelism: 15,000 genes
    │
    ├─ Gene 1: U-Net (upward → downward)
    ├─ Gene 2: U-Net (upward → downward)
    ├─ ...
    └─ Gene 15,000: U-Net (upward → downward)

Per-gene U-Net:
    Level 1 batches: ~70 local neighborhoods (20 spots each)
    Level 2 batches: ~14 regions (100 spots each)
    Level 3: 1 global solve (1400 spots)
```

### Gurobi Settings for Large Batches

```python
model.setParam("Method", 2)           # Barrier method
model.setParam("BarConvTol", 1e-6)    # Slightly relaxed tolerance
model.setParam("Threads", 8)          # Parallel barrier
model.setParam("Presolve", 2)         # Aggressive presolve
model.setParam("NumericFocus", 1)     # Better numerical stability
```

### Per-Gene Solve Time (Estimated)

| Level | Spots | Variables | Gurobi time |
|-------|-------|-----------|-------------|
| 1 (local) | 20 | 140 | <0.1s |
| 2 (region) | 100 | 700 | ~0.5s |
| 3 (global) | 1400 | 9,800 | ~2-5s |

**Total per gene:** ~10-15s (both passes)

### Wall Clock with Parallelism

- 15,000 genes ÷ 64 cores = ~235 batches
- 235 × 15s = **~1 hour total**

### Memory per Worker

- One gene's data: 1400 spots × 7 types × 8 bytes = ~80 KB
- Gurobi model overhead: ~10-50 MB
- Comfortable with 2-4 GB per worker

---

## Implementation Outline

### New Module Structure

```
CITEgeist/model/
├── gurobi_impl.py              # Existing (keep for reference)
├── spatial_unet_gex.py         # NEW: U-Net GEX deconvolution
│   ├── build_spatial_hierarchy()
│   ├── solve_gene_unet()
│   ├── upward_pass_level()
│   ├── downward_pass_level()
│   ├── pool_to_region()
│   └── run_unet_gex_deconvolution()  # Main entry point
└── citegeist_model.py          # Update run_cell_expression_pass1() to use U-Net
```

### Key Functions

```python
def run_unet_gex_deconvolution(
    adata: sc.AnnData,
    cell_type_proportions: np.ndarray,
    lambda_spatial: float = 0.1,
    lambda_prior: float = 0.5,
    n_levels: int = 3,
    max_workers: int = 64,
) -> Dict[str, np.ndarray]:
    """
    Run gene-parallel U-Net GEX deconvolution.

    Returns:
        Dict with keys per cell type, values are (n_spots, n_genes) arrays
    """

def solve_gene_unet(
    gene_idx: int,
    observed: np.ndarray,          # (n_spots,) observed counts
    proportions: np.ndarray,       # (n_spots, n_celltypes)
    enrichment: np.ndarray,        # (n_spots, n_celltypes)
    hierarchy: SpatialHierarchy,
    lambda_spatial: float,
    lambda_prior: float,
) -> np.ndarray:
    """
    Solve U-Net for a single gene.

    Returns:
        (n_spots, n_celltypes) deconvolved expression
    """
```

---

## Testing Strategy

### Unit Tests

1. `test_asymmetric_laplacian_weights()` - verify w_ij formula
2. `test_count_conservation()` - verify Σ_t X[i,t] = observed[i,g]
3. `test_hierarchy_pooling()` - verify proportion-weighted pooling
4. `test_prior_shrinkage()` - verify low-confidence spots shrink more

### Integration Tests

1. Run on simulated data with known ground truth
2. Compare Pearson r before/after U-Net
3. Verify regional heterogeneity preserved (inject two subclones)

### Benchmark Validation

1. Run on Xenium benchmark (5 regions, 7 cell types)
2. Compare to current Module 3 Pass 2 results
3. Target: improve mean Pearson r from ~0.35 to ~0.50+

---

## Risk Assessment

| Risk | Mitigation |
|------|------------|
| Over-smoothing destroys signal | CV tunes λ_spatial; hierarchy preserves regional patterns |
| Computational cost too high | Gene-parallel keeps each QP small (~10k vars) |
| Memory blowup | Stream results to disk; 2-4 GB per worker sufficient |
| CV overfits to hold-out | Use spatial stratification; test on separate regions |
| Regional boundaries misaligned | Pooling is soft (proportion-weighted), not hard partitions |

---

## Success Criteria

- [ ] Mean Pearson r improves by ≥0.10 across cell types
- [ ] Rare cell types (B cells, CD4+ T) improve most (currently worst)
- [ ] Regional heterogeneity preserved (visual inspection of subclones)
- [ ] Runtime ≤2 hours on 64 cores for full gene set
- [ ] Memory ≤4 GB per worker

---

## Next Steps

1. Implement `spatial_unet_gex.py` with core functions
2. Add unit tests for Laplacian weights and pooling
3. Run on Xenium benchmark region 0 as proof-of-concept
4. Tune hyperparameters via CV
5. Full benchmark comparison across all regions
