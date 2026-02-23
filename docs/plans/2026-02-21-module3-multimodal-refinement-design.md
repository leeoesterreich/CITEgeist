# Module 3 Multimodal Refinement Design

**Date:** 2026-02-21
**Status:** Implemented (2026-02-23)
**Author:** Brainstorming session

## Problem Statement

CITEgeist's current Module 3 uses protein markers exclusively for cell type proportion estimation (Pass 1), then deconvolves gene expression (Pass 2) using those proportions. This works well when protein markers are reliable, but fails for ~22% of cells that have clear RNA signatures but low protein expression.

**Analysis showed:**
- 54% of RNA-Epithelial cells have no detectable epithelial protein markers (PanCK/E-Cadherin)
- 39% of RNA-Fibroblasts have no alphaSMA protein signal
- These cells contribute RNA signal to spots but CITEgeist can't assign them

**Result:** CITEgeist underperforms RNA-based methods (Cell2Location, RCTD) when evaluated against RNA-defined ground truth.

## Design Goals

1. **Reference-free**: No external scRNA-seq reference required (core CITEgeist value)
2. **Protein-first**: Preserve protein's role as the primary signal; RNA refines, not replaces
3. **Modular**: Minimal changes to existing Pass 1; enhancements in Pass 2
4. **Principled**: Use EM framework already present in codebase

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│ PASS 1 (existing, unchanged)                                    │
│   Input: Protein markers                                        │
│   EM: Protein signal ↔ Y_protein                                │
│   Output: Y_protein (N_spots × T_celltypes)                     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ PASS 1.5 (new): Learn Anchor Gene Signatures                    │
│                                                                 │
│   1. Identify confident spots per cell type:                    │
│      confident_spots[t] = spots where Y_protein[:, t] > 0.3     │
│                                                                 │
│   2. For each cell type t, rank genes by correlation:           │
│      r[g, t] = pearson(GEX[:, g], Y_protein[:, t])              │
│                                                                 │
│   3. Compute specificity (distinctiveness):                     │
│      specificity[g, t] = r[g, t] - max(r[g, other_types])       │
│                                                                 │
│   4. Select top-N anchor genes (N=20) per cell type:            │
│      anchors[t] = top_N_genes_by(r * specificity)               │
│      weights[g, t] = r[g, t]  (correlation as weight)           │
│                                                                 │
│   Output: anchor_genes (T × N), weights (T × N)                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ PASS 2 (enhanced): EM with RNA Refinement                       │
│                                                                 │
│   Initialize:                                                   │
│     Y = Y_protein                                               │
│     E = zeros(T_celltypes × G_genes)                            │
│                                                                 │
│   EM Loop until convergence:                                    │
│   ┌─────────────────────────────────────────────────────────┐   │
│   │ E-STEP: Estimate gene expression profiles E             │   │
│   │                                                         │   │
│   │   For anchor genes (LOCKED):                            │   │
│   │     E[t, g] = weighted_mean(GEX[:, g], weights=Y[:, t]) │   │
│   │     (gene g is assigned to cell type t, fixed)          │   │
│   │                                                         │   │
│   │   For non-anchor genes (FREE):                          │   │
│   │     E[:, g] = solve least squares:                      │   │
│   │       minimize ||GEX[:, g] - Y @ E[:, g]||²             │   │
│   │     (gene g can load on any cell type)                  │   │
│   └─────────────────────────────────────────────────────────┘   │
│                              │                                  │
│                              ▼                                  │
│   ┌─────────────────────────────────────────────────────────┐   │
│   │ M-STEP: Refine proportions Y                            │   │
│   │                                                         │   │
│   │   Objective (per spot i):                               │   │
│   │     minimize:                                           │   │
│   │       Σ_g w[g] * (GEX[i,g] - Y[i,:] @ E[:,g])²         │   │
│   │       + λ * ||Y[i,:] - Y_protein[i,:]||²    (prior)     │   │
│   │                                                         │   │
│   │   Constraints:                                          │   │
│   │     Y[i, :] >= 0                                        │   │
│   │     sum(Y[i, :]) <= 1                                   │   │
│   │                                                         │   │
│   │   Where:                                                │   │
│   │     w[g] = anchor_weight[g] if anchor, else 1.0         │   │
│   │     λ = prior strength (trust in protein)               │   │
│   └─────────────────────────────────────────────────────────┘   │
│                                                                 │
│   Convergence: ||Y_new - Y_old|| < tolerance                    │
│                                                                 │
│   Output: Y_refined, E_final                                    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ PASS 2b (existing): Final GEX Deconvolution                     │
│   Use Y_refined to produce cell-type-specific expression layers │
└─────────────────────────────────────────────────────────────────┘
```

## Key Design Decisions

### 1. Anchor Gene Selection (Pass 1.5)

**Method:** Pearson correlation between gene expression and cell type proportion

```python
def select_anchor_genes(GEX, Y_protein, n_anchors=20, min_correlation=0.3):
    """
    Select top-N anchor genes per cell type based on correlation and specificity.

    Args:
        GEX: Gene expression matrix (N_spots × G_genes)
        Y_protein: Cell type proportions from Pass 1 (N_spots × T_types)
        n_anchors: Number of anchor genes per cell type
        min_correlation: Minimum correlation to be considered

    Returns:
        anchors: Dict[cell_type] -> List[gene_names]
        weights: Dict[cell_type] -> Dict[gene] -> weight
    """
    n_spots, n_genes = GEX.shape
    n_types = Y_protein.shape[1]

    # Compute correlation matrix (genes × cell types)
    correlations = np.zeros((n_genes, n_types))
    for t in range(n_types):
        for g in range(n_genes):
            correlations[g, t] = pearsonr(GEX[:, g], Y_protein[:, t])[0]

    # Compute specificity: how much better is best type vs second best
    specificity = np.zeros((n_genes, n_types))
    for g in range(n_genes):
        for t in range(n_types):
            other_max = np.max(np.delete(correlations[g, :], t))
            specificity[g, t] = correlations[g, t] - other_max

    # Combined score
    score = correlations * np.clip(specificity, 0, None)

    # Select top-N per cell type
    anchors = {}
    weights = {}
    for t in range(n_types):
        valid = correlations[:, t] > min_correlation
        ranked = np.argsort(score[valid, t])[::-1][:n_anchors]
        anchors[t] = gene_names[valid][ranked]
        weights[t] = {g: correlations[gene_idx, t] for g, gene_idx in zip(anchors[t], ranked)}

    return anchors, weights
```

### 2. Anchor Locking Strategy

**Decision:** Anchors are LOCKED during EM iterations.

**Rationale:**
- Anchors come from confident (high-proportion) spots - we trust them
- Locking provides stability during optimization
- Non-anchor genes provide flexibility to capture heterogeneity
- Post-hoc validation can flag problematic anchors

### 3. Weighting Scheme

**Anchor gene weight = Pearson correlation coefficient**

- High r (e.g., 0.8) → gene strongly tracks cell type → high weight in objective
- Low r (e.g., 0.3) → weak marker → low weight
- Non-anchor genes → weight = 1.0 (baseline)

### 4. Prior Strength (λ)

The parameter λ controls trust in protein vs RNA:

| λ value | Behavior |
|---------|----------|
| λ = 0 | Pure RNA refinement, ignore protein |
| λ = 1 | Balance protein and RNA equally |
| λ = 10 | Strong protein prior, RNA makes small corrections |
| λ = ∞ | No refinement, Y = Y_protein |

**Default recommendation:** λ = 1.0, tunable per dataset.

## Implementation Plan

### New Functions

```python
# In gurobi_impl.py or new file: multimodal_refinement.py

def select_anchor_genes(GEX, Y_protein, n_anchors=20, min_correlation=0.3):
    """Pass 1.5: Learn anchor gene signatures from confident spots."""
    ...

def compute_expression_profiles(GEX, Y, anchors, weights):
    """E-step: Estimate E given current Y. Anchors locked, others free."""
    ...

def refine_proportions(GEX, Y_current, E, Y_protein, anchors, weights, lambda_prior=1.0):
    """M-step: Update Y given current E, with prior toward Y_protein."""
    ...

def multimodal_em_refinement(GEX, Y_protein, n_anchors=20, lambda_prior=1.0,
                              max_iterations=20, tolerance=1e-4):
    """Full Pass 2 EM with RNA refinement."""
    # Pass 1.5
    anchors, weights = select_anchor_genes(GEX, Y_protein, n_anchors)

    # Initialize
    Y = Y_protein.copy()

    # EM loop
    for iteration in range(max_iterations):
        # E-step
        E = compute_expression_profiles(GEX, Y, anchors, weights)

        # M-step
        Y_new = refine_proportions(GEX, Y, E, Y_protein, anchors, weights, lambda_prior)

        # Check convergence
        if np.max(np.abs(Y_new - Y)) < tolerance:
            break
        Y = Y_new

    return Y, E, anchors
```

### Integration with CitegeistModel

```python
# In citegeist_model.py

def run_cell_proportion_model(self, ...):
    """Pass 1: Existing protein-based optimization."""
    # ... existing code ...

def run_multimodal_refinement(self, n_anchors=20, lambda_prior=1.0):
    """Pass 1.5 + Pass 2 EM: RNA-based refinement."""
    from .multimodal_refinement import multimodal_em_refinement

    Y_protein = self.cell_proportions  # From Pass 1
    GEX = self.gene_expression_adata.X

    Y_refined, E, anchors = multimodal_em_refinement(
        GEX, Y_protein, n_anchors, lambda_prior
    )

    self.cell_proportions_refined = Y_refined
    self.expression_profiles = E
    self.anchor_genes = anchors

def run_cell_expression_pass1(self, use_refined=True, ...):
    """Pass 2b: Final GEX deconvolution."""
    if use_refined and hasattr(self, 'cell_proportions_refined'):
        proportions = self.cell_proportions_refined
    else:
        proportions = self.cell_proportions
    # ... existing deconvolution code ...
```

## Expected Outcomes

### Performance Improvements

| Scenario | Current | Expected |
|----------|---------|----------|
| Protein GT benchmark | 0.73 | 0.73 (unchanged) |
| RNA GT benchmark | 0.51 | 0.60-0.65 |
| Cell types with high Unknown % | Poor | Improved |

### Discovered Gene Markers

The EM process will identify genes that track with anchor genes, effectively discovering:
- Additional markers for each cell type
- Cell-type-specific expression signatures
- Potential novel markers not in standard panels

## Testing Strategy

1. **Unit tests:** Each new function with synthetic data
2. **Simulation benchmark:** Run on existing mixed/high_seg simulations
3. **Xenium benchmark:** Compare against protein GT and RNA GT
4. **Ablation studies:**
   - Effect of n_anchors (5, 10, 20, 50)
   - Effect of lambda_prior (0.1, 1.0, 10.0)
   - Locked vs unlocked anchors

## Open Questions

1. **Convergence:** Will EM always converge? May need damping or max iterations.
2. **Scalability:** With ~400 genes × ~1000 spots, should be tractable. May need batching for larger datasets.
3. **Negative correlations:** Should genes anti-correlated with a cell type be used as negative markers?

## Summary

This design introduces **multimodal refinement** to CITEgeist Module 3:

1. **Pass 1.5** learns anchor gene signatures from protein-confident spots
2. **Pass 2 EM** iteratively refines proportions using both protein prior and RNA reconstruction
3. **Anchor locking** with **correlation-based weighting** provides stability

The approach is reference-free, preserves protein-first philosophy, and addresses the gap where cells have RNA signal but low protein expression.
