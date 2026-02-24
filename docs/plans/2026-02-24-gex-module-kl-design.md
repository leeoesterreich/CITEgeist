# GEX Deconvolution: Module-Aware Enrichment + Softmax KL Regularization

**Date:** 2026-02-24
**Status:** Approved
**Author:** Claude + Alex

## Problem Statement

Current GEX deconvolution uses L2 regularization which creates uniform spreading of gene counts across cell types, resulting in:
- 3x spurious variance (variance ratio = 3.02)
- Poor spatial pattern preservation for general genes (per-gene r = 0.273)
- Large gap between marker genes (r = 0.64) and general genes (r = 0.27)

## Solution Overview

Two integrated improvements to the GEX deconvolution pipeline:

1. **Module-Aware Enrichment**: Use spatial correlation with data-driven anchor genes to boost enrichment for weak genes
2. **Softmax KL Regularization**: Replace L2 penalty with KL-divergence from softmax target distribution

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                    GEX Deconvolution Pipeline                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────────┐                                           │
│  │ Cell Proportions │ (from Pass 1)                             │
│  │   (N × T)        │                                           │
│  └────────┬─────────┘                                           │
│           │                                                      │
│           ▼                                                      │
│  ┌──────────────────────────────────────┐                       │
│  │  STAGE 1: Anchor Gene Discovery      │  (run once globally)  │
│  │  ─────────────────────────────────   │                       │
│  │  For each cell type t:               │                       │
│  │    anchor_genes[t] = top-K genes     │                       │
│  │    correlated with proportion[t]     │                       │
│  └────────┬─────────────────────────────┘                       │
│           │                                                      │
│           ▼                                                      │
│  ┌──────────────────────────────────────┐                       │
│  │  STAGE 2: Module-Aware Enrichment    │  (per spot)           │
│  │  ─────────────────────────────────   │                       │
│  │  For each gene g in neighborhood:    │                       │
│  │    base_enrichment = existing method │                       │
│  │    module_boost = corr(g, anchors)   │                       │
│  │    adjusted_enrichment = combine()   │                       │
│  └────────┬─────────────────────────────┘                       │
│           │                                                      │
│           ▼                                                      │
│  ┌──────────────────────────────────────┐                       │
│  │  STAGE 3: Softmax KL Allocation      │  (per spot)           │
│  │  ─────────────────────────────────   │                       │
│  │  target = softmax(adjusted_enrich/T) │                       │
│  │  objective = enrichment_term         │                       │
│  │             - λ_KL × KL(X || target) │                       │
│  └────────┬─────────────────────────────┘                       │
│           │                                                      │
│           ▼                                                      │
│  ┌──────────────────┐                                           │
│  │  Deconvolved GEX │                                           │
│  │   (N × T × M)    │                                           │
│  └──────────────────┘                                           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Stage 1: Anchor Gene Discovery

**Purpose:** Identify genes whose spatial expression pattern reliably tracks each cell type's proportion map.

**Algorithm:**

```python
def discover_anchor_genes(
    gene_expression: np.ndarray,    # (N_spots × M_genes)
    cell_proportions: np.ndarray,   # (N_spots × T_celltypes)
    min_anchors: int = 5,
    max_anchors: int = 10,
    initial_min_correlation: float = 0.3,
    min_expressing_spots: float = 0.1,
) -> Tuple[Dict[int, List[int]], Dict[int, float]]:
    """
    Returns:
        anchors: {cell_type_idx: [gene_idx1, gene_idx2, ...]}
        thresholds_used: {cell_type_idx: correlation_threshold_used}
    """

    anchors = {}
    thresholds_used = {}

    # Threshold levels: 0.30 → 0.25 → 0.20 → 0.15 → 0.10
    threshold_sequence = [0.30, 0.25, 0.20, 0.15, 0.10]

    for t in range(T):
        prop_vector = cell_proportions[:, t]

        # Compute all correlations once
        all_correlations = []
        for g in range(M):
            gene_vector = gene_expression[:, g]

            if (gene_vector > 0).mean() < min_expressing_spots:
                continue

            r, p = pearsonr(gene_vector, prop_vector)
            if p < 0.05:
                all_correlations.append((g, r))

        all_correlations.sort(key=lambda x: -x[1])

        # Step through thresholds until we have min_anchors
        selected_threshold = threshold_sequence[-1]  # default to floor
        for threshold in threshold_sequence:
            candidates = [g for g, r in all_correlations if r >= threshold]
            if len(candidates) >= min_anchors:
                selected_threshold = threshold
                break

        # Select anchors at chosen threshold (capped at max)
        final_candidates = [g for g, r in all_correlations if r >= selected_threshold]
        anchors[t] = final_candidates[:max_anchors]
        thresholds_used[t] = selected_threshold

    return anchors, thresholds_used
```

**Parameters:**

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `min_anchors` | 5 | Minimum for robust module signal |
| `max_anchors` | 10 | Cap to avoid dilution |
| `initial_min_correlation` | 0.3 | Start with strong requirement |
| `threshold_sequence` | [0.30, 0.25, 0.20, 0.15, 0.10] | Step down in 0.05 increments |

**Adaptive behavior:** Threshold drops in 0.05 increments until cell type has at least 5 anchors.

## Stage 2: Module-Aware Enrichment

**Purpose:** For each gene at each spot, boost enrichment toward cell types whose anchor genes are spatially co-expressed in that neighborhood.

**Algorithm:**

```python
def compute_module_aware_enrichment(
    spot_idx: int,
    neighborhood_expression: np.ndarray,  # (N_neighbors × M_genes)
    base_enrichment: np.ndarray,           # (M_genes × T_celltypes)
    anchor_genes: Dict[int, List[int]],    # {celltype: [gene_indices]}
    module_weight: float = 0.5,
    min_neighbors_for_corr: int = 10,
) -> np.ndarray:
    """
    Returns: adjusted_enrichment (M_genes × T_celltypes)
    """
    M, T = base_enrichment.shape
    adjusted = base_enrichment.copy()

    # Skip module adjustment if neighborhood too small
    if neighborhood_expression.shape[0] < min_neighbors_for_corr:
        return adjusted

    for g in range(M):
        gene_expr = neighborhood_expression[:, g]

        if np.std(gene_expr) < 1e-6:
            continue

        # Compute correlation with each cell type's anchor genes
        module_scores = np.zeros(T)
        for t in range(T):
            if not anchor_genes.get(t):
                continue

            anchor_corrs = []
            for anchor_idx in anchor_genes[t]:
                anchor_expr = neighborhood_expression[:, anchor_idx]
                if np.std(anchor_expr) > 1e-6:
                    r, _ = pearsonr(gene_expr, anchor_expr)
                    if not np.isnan(r):
                        anchor_corrs.append(max(0, r))  # only positive correlations

            if anchor_corrs:
                module_scores[t] = np.mean(anchor_corrs)

        # Normalize module scores
        if module_scores.sum() > 0:
            module_scores = module_scores / module_scores.sum()
        else:
            module_scores = np.ones(T) / T

        # Combine base enrichment with module signal
        adjusted[g, :] = (
            (1 - module_weight) * base_enrichment[g, :] +
            module_weight * module_scores
        )
        adjusted[g, :] = adjusted[g, :] / (adjusted[g, :].sum() + 1e-10)

    return adjusted
```

**Parameters:**

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `module_weight` | 0.5 | Equal blend of base enrichment and module signal |
| `min_neighbors_for_corr` | 10 | Need enough points for meaningful correlation |

**Key design choice:** Only positive correlations boost enrichment. Negative correlation means the gene anti-correlates with that cell type's anchors.

## Stage 3: Softmax KL-Divergence Allocation

**Purpose:** Replace L2 regularization with KL-divergence that pulls allocations toward the adjusted enrichment distribution.

**Current objective (L2):**
```
maximize: Σ_j [ enrichment[j] × prop[j] × X[j] ] - λ_L2 × Σ_j X[j]²
```

**New objective (KL):**
```
maximize: Σ_j [ enrichment[j] × prop[j] × X[j] ] - λ_KL × KL(X_norm || p_target)
```

Where: `p_target = softmax(adjusted_enrichment / temperature)`

**Algorithm:**

```python
def build_gex_objective_with_kl(
    model: gp.Model,
    X: Dict[Tuple[int,int], gp.Var],
    gene_idx: int,
    total_counts: int,
    adjusted_enrichment: np.ndarray,   # (T,) for this gene
    center_props: np.ndarray,          # (T,) cell type proportions
    temperature: float = 0.3,
    lambda_kl: float = 0.1,
) -> List:
    """
    Returns list of objective terms for one gene.
    """
    T = len(center_props)
    obj_terms = []

    # Compute softmax target distribution
    logits = adjusted_enrichment / temperature
    logits = logits - logits.max()  # numerical stability
    p_target = np.exp(logits) / np.exp(logits).sum()
    p_target = np.clip(p_target, 1e-6, 1.0)
    p_target = p_target / p_target.sum()

    for j in range(T):
        # Term 1: Enrichment × proportion (maximize)
        base_term = adjusted_enrichment[j] * center_props[j] * X[j, gene_idx]
        obj_terms.append(base_term)

        # Term 2: KL divergence approximation
        # Penalize deviation from target allocation
        target_count = p_target[j] * total_counts
        deviation = X[j, gene_idx] - target_count
        kl_penalty = lambda_kl * deviation * deviation / (total_counts + 1)
        obj_terms.append(-kl_penalty)

    return obj_terms
```

**Parameters:**

| Parameter | Default | Rationale |
|-----------|---------|-----------|
| `temperature` | 0.3 | Sharper target - favors top enriched cell types |
| `lambda_kl` | 0.1 | Strength of pull toward target |

**Temperature effect:**
```
adjusted_enrichment = [0.4, 0.3, 0.1, 0.1, 0.1]

T=0.3 (sharp):  p_target = [0.47, 0.35, 0.06, 0.06, 0.06]
T=0.5 (medium): p_target = [0.39, 0.32, 0.10, 0.10, 0.10]
T=1.0 (soft):   p_target = [0.30, 0.27, 0.14, 0.14, 0.14]
```

## Integration

**Files to modify:**

| File | Changes |
|------|---------|
| `CITEgeist/model/gurobi_impl.py` | Add `discover_anchor_genes()`, `compute_module_aware_enrichment()`, `build_gex_objective_with_kl()`. Modify `deconvolute_spot_with_neighbors_with_prior()` |
| `CITEgeist/model/citegeist_model.py` | Add parameters to `run_cell_expression_pass1()` |

**New parameters for `run_cell_expression_pass1()`:**

```python
def run_cell_expression_pass1(
    self,
    # ... existing params ...

    # NEW parameters
    use_module_enrichment: bool = True,
    module_weight: float = 0.5,
    use_kl_regularization: bool = True,
    kl_temperature: float = 0.3,
    lambda_kl: float = 0.1,
    n_anchor_genes: Tuple[int, int] = (5, 10),  # (min, max)

    # DEPRECATED (kept for backward compat)
    lambda_gex_reg: float = 0.01,
):
```

**Backward compatibility:** Setting `use_module_enrichment=False` and `use_kl_regularization=False` reproduces old L2 behavior.

## Testing Strategy

**Unit tests:**

| Test | Purpose |
|------|---------|
| `test_anchor_discovery_finds_known_markers` | Verify recovery of canonical markers |
| `test_anchor_threshold_adaptation` | Verify threshold drops until min anchors |
| `test_module_enrichment_boosts_correlated_genes` | Verify module boost mechanism |
| `test_kl_objective_concentrates_allocation` | Verify KL vs L2 behavior |
| `test_backward_compatibility` | Verify old behavior when disabled |

**Ablation tests:**

| Variant | Module | KL | Purpose |
|---------|--------|-----|---------|
| Baseline | No | No | Current L2 (control) |
| KL-only | No | Yes | Isolate KL effect |
| Module-only | Yes | No | Isolate module effect |
| Full | Yes | Yes | Combined (expected best) |

## Success Criteria

| Metric | Baseline L2 | Target | Stretch Goal |
|--------|-------------|--------|--------------|
| Per-gene spatial r | 0.273 | > 0.30 | > 0.35 |
| Marker gene r | 0.636 | >= 0.636 | > 0.70 |
| General gene r | ~0.20 | > 0.25 | > 0.30 |
| Variance ratio | 3.02 | 1.0 - 1.5 | ~1.0 |

Metric: Pearson correlation between predicted and ground truth per spot per cell type.

## Implementation Plan

1. Implement `discover_anchor_genes()` in `gurobi_impl.py`
2. Implement `compute_module_aware_enrichment()` in `gurobi_impl.py`
3. Implement `build_gex_objective_with_kl()` in `gurobi_impl.py`
4. Modify `deconvolute_spot_with_neighbors_with_prior()` to use new components
5. Add parameters to `run_cell_expression_pass1()` in `citegeist_model.py`
6. Write unit tests
7. Run Xenium benchmark ablation (4 variants)
8. Evaluate against success criteria
