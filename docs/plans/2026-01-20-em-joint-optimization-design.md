# EM-Style Joint Optimization for CITEgeist

**Date:** 2026-01-20
**Status:** Design approved, ready for implementation
**Branch:** `hierarchical_approach`

## Problem Statement

CITEgeist currently performs two sequential passes:
1. **Pass 1:** Protein markers → Cell type proportions
2. **Pass 2:** Proportions → Gene expression deconvolution

This sequential approach has limitations:
- Gene expression information is not used to refine proportions
- Errors in Pass 1 propagate to Pass 2
- No feedback loop between modalities

**Current benchmark performance (Xenium, 5 regions):**

| Metric | CITEgeist | Cell2Location | RCTD |
|--------|-----------|---------------|------|
| Proportion Pearson r | 0.600 | 0.611 | 0.619 |
| Proportion RMSE | 0.167 | 0.179 | 0.177 |
| GEX RMSE | 0.836 | 0.477 | N/A |

Cell2Location's GEX advantage comes from external scRNA-seq reference signatures. CITEgeist must work without external references (protein + gene data from same slide only).

## Solution: EM-Style Joint Optimization

Learn gene-cell type associations from **spatial co-variation** with protein markers, then use these to refine both proportions and gene assignments iteratively.

### Key Insight

Genes that spatially co-vary with cell-type-specific proteins are likely expressed by those cell types. This lets us "learn" pseudo-signatures from the slide itself, varying by local neighborhood (capturing microenvironment heterogeneity).

### Algorithm Overview

```
┌─────────────────────────────────────────────────────────┐
│                    INITIALIZATION                        │
│  Run current protein-based Pass 1 → initial proportions │
└─────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────┐
│                 E-STEP (per neighborhood)                │
│  For each spot i and its neighborhood N(i):             │
│    R[i,t,g] = "responsibility of type t for gene g      │
│                in neighborhood N(i)"                     │
│  Based on LOCAL correlation between gene g and type t   │
└─────────────────────────────────────────────────────────┘
                           ↓
┌─────────────────────────────────────────────────────────┐
│                 M-STEP (per neighborhood)                │
│  For each spot i, refine proportions P[i,:] using:      │
│    - Protein signal (anchors cell type identity)        │
│    - Gene evidence weighted by LOCAL responsibilities   │
└─────────────────────────────────────────────────────────┘
                           ↓
                    Converged? ──No──→ (back to E-step)
                           ↓
┌─────────────────────────────────────────────────────────┐
│                   FINAL OUTPUT                           │
│  Spatially-refined proportions P*                       │
│  Local gene assignments R* (varies across tissue)       │
└─────────────────────────────────────────────────────────┘
```

**Key difference from Cell2Location:** They use fixed global signatures. We learn **spatially-varying** signatures, capturing microenvironment heterogeneity.

## E-Step: Local Gene Responsibility Estimation

**Goal:** For each spot's neighborhood, estimate how much each cell type "owns" each gene based on local spatial co-variation.

**Inputs:**
- Current proportions `P[i, t]` for all spots i and cell types t
- Gene expression `G[i, g]` for all spots i and genes g
- Neighborhood structure `N(i)` for each spot

**Algorithm:**

```python
def compute_local_responsibilities(
    spot_idx: int,
    P: np.ndarray,           # (N_spots, T) proportions
    G: np.ndarray,           # (N_spots, M) gene expression
    neighborhood: List[int], # indices of neighbors
    temperature: float,
) -> np.ndarray:
    """
    Compute R[t, g] for this spot's neighborhood.

    Returns:
        R: (T, M) responsibility matrix
    """
    # 1. Extract local data
    P_local = P[neighborhood, :]  # (n_neighbors, T)
    G_local = G[neighborhood, :]  # (n_neighbors, M)

    T, M = P_local.shape[1], G_local.shape[1]
    corr = np.zeros((T, M))

    # 2. Compute correlation between each gene and cell type
    for g in range(M):
        for t in range(T):
            corr[t, g] = pearsonr(G_local[:, g], P_local[:, t])[0]

    # Handle NaN correlations (constant expression)
    corr = np.nan_to_num(corr, nan=0.0)

    # 3. Convert to responsibilities via softmax with temperature
    # Negative correlations → near-zero responsibility
    corr_clipped = np.maximum(corr, 0)
    R = softmax(corr_clipped / temperature, axis=0)

    return R
```

**Temperature annealing schedule:**

| Iteration | Temperature | Behavior |
|-----------|-------------|----------|
| 1 | 1.0 | Soft assignments, explore broadly |
| 2 | 0.7 | Start sharpening |
| 3 | 0.49 | More confident |
| ... | ... | ... |
| N | ~0.1 | Near winner-take-all |

Formula: `T(k) = T_initial * (cooling_rate ^ k)`

## M-Step: Proportion Refinement with Dual Evidence

**Goal:** Update cell type proportions using BOTH protein markers AND gene expression evidence.

**Combined objective:**

```python
def optimize_with_dual_evidence(
    spot_idx: int,
    A: np.ndarray,              # (N_spots, P) protein data
    G: np.ndarray,              # (N_spots, M) gene expression
    R: np.ndarray,              # (T, M) local responsibilities
    protein_signatures: np.ndarray,  # (T, P) protein profiles
    alpha: float = 0.7,         # protein weight
    beta: float = 0.3,          # gene weight
    lambda_spatial: float = 0.1,
) -> np.ndarray:
    """
    Optimize proportions for one spot using dual evidence.

    Returns:
        P_opt: (T,) optimized proportions
    """
    # PROTEIN TERM: ||A[i,:] - P @ ProteinSignatures||²
    # (Current CITEgeist approach)

    # GENE TERM: For each gene g
    #   Expected[g] = Σ_t P[t] * R[t,g] * scale[t,g]
    #   Error = ||G[i,g] - Expected[g]||²

    # COMBINED OBJECTIVE
    # loss = α * protein_error + β * gene_error + λ * spatial_smoothness

    # CONSTRAINTS
    # P >= 0, sum(P) = 1

    # Solve via Gurobi QP (existing infrastructure)
```

**Weighting rationale:**
- `α = 0.7` (protein): Direct cell type markers, high confidence
- `β = 0.3` (gene): Indirect evidence via learned associations, noisier
- Proteins anchor identity; genes provide supporting evidence

## Convergence & Iteration Control

```python
def run_em_refinement(
    P_init: np.ndarray,
    A: np.ndarray,
    G: np.ndarray,
    neighborhoods: List[List[int]],
    max_iterations: int = 10,
    initial_temperature: float = 1.0,
    cooling_rate: float = 0.7,
    alpha: float = 0.7,
    beta: float = 0.3,
    tol_p: float = 0.01,
    tol_r: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Main EM loop.
    """
    P = P_init.copy()
    temperature = initial_temperature
    prev_loss = float('inf')

    for iteration in range(max_iterations):
        # E-step: compute responsibilities for all spots
        R = np.zeros((N_spots, T, M))
        for i in range(N_spots):
            R[i] = compute_local_responsibilities(
                i, P, G, neighborhoods[i], temperature
            )

        # M-step: optimize proportions for all spots
        P_new = np.zeros_like(P)
        for i in range(N_spots):
            P_new[i] = optimize_with_dual_evidence(
                i, A, G, R[i], protein_signatures, alpha, beta
            )

        # Check convergence
        proportion_change = np.mean(np.abs(P_new - P))

        if proportion_change < tol_p:
            logging.info(f"Converged at iteration {iteration}")
            break

        # Early stopping: revert if loss increases
        current_loss = compute_total_loss(P_new, A, G, R, alpha, beta)
        if current_loss > prev_loss * 1.01:  # 1% tolerance
            logging.warning("Loss increased, reverting and stopping")
            break

        # Anneal temperature
        temperature = temperature * cooling_rate
        P = P_new
        prev_loss = current_loss

    return P, R
```

**Default parameters:**

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| max_iterations | 10-15 | Usually converges in 5-8 |
| initial_temperature | 1.0 | Soft start |
| cooling_rate | 0.7 | Geometric decay |
| tol_p | 0.01 | 1% proportion change threshold |
| alpha (protein) | 0.7 | Protein-anchored |
| beta (gene) | 0.3 | Gene as supporting evidence |

## Integration with Existing Codebase

**Pipeline integration:**

```
Current:
  Module 1-2 → Profiles → Pass 1 (Proportions) → Pass 2 (GEX) → Done

Proposed:
  Module 1-2 → Profiles → Pass 1 (Init) → EM Joint Refinement → Done
```

**New module structure:**

```
CITEgeist/model/
├── gurobi_impl.py          # Existing - keep Pass 1 for initialization
├── joint_optimization.py   # NEW - EM algorithm
│   ├── compute_local_responsibilities()   # E-step
│   ├── optimize_with_dual_evidence()      # M-step
│   └── run_em_refinement()                # Main loop
└── citegeist_model.py      # Add: run_joint_optimization()
```

**API addition:**

```python
class CitegeistModel:
    def run_joint_optimization(
        self,
        max_iterations: int = 10,
        initial_temperature: float = 1.0,
        cooling_rate: float = 0.7,
        alpha_protein: float = 0.7,
        beta_gene: float = 0.3,
        convergence_tol: float = 0.01,
        radius: float = 4.0,
    ) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
        """
        Run EM-style joint optimization of proportions and gene assignments.

        Returns:
            proportions: Refined cell type proportions (N_spots × T)
            gene_assignments: Local gene responsibilities (N_spots × T × M)
        """
```

**Backward compatibility:**
- Existing `run_cell_proportion_model()` and `run_cell_expression_pass1()` unchanged
- New method is opt-in addition
- Benchmark can compare both approaches

## Expected Outcomes

**Target improvements:**

| Metric | Current | Target | Mechanism |
|--------|---------|--------|-----------|
| Proportion Pearson r | 0.60 | 0.62-0.65 | Gene evidence reinforces protein signal |
| GEX RMSE | 0.84 | 0.55-0.65 | Learned local signatures vs none |
| JSD | 0.355 | 0.33-0.34 | Better calibrated distributions |

**Why this should help GEX:**
- Current: Assigns genes using only enrichment heuristics
- Proposed: Learns gene-cell type associations from spatial co-variation with proteins
- Mimics Cell2Location's reference advantage, but learned locally

## Risks & Mitigations

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Slower runtime | Medium | Parallelize E-step; limit iterations |
| Overfitting in small neighborhoods | Medium | Min neighborhood size; regularization |
| Divergence/instability | Low | Early stopping; loss monitoring |
| No improvement | Low | Keep original pipeline as fallback |

## Validation Plan

1. Run on Xenium benchmark (5 regions, 6 cell types)
2. Compare metrics against:
   - Current CITEgeist (baseline)
   - Cell2Location (GEX target)
   - RCTD (proportion target)
3. Verify convergence (should stabilize by iteration 5-8)
4. Check spatial coherence of learned gene assignments

## Implementation Checklist

- [ ] Create `joint_optimization.py` module
- [ ] Implement `compute_local_responsibilities()` (E-step)
- [ ] Implement `optimize_with_dual_evidence()` (M-step)
- [ ] Implement `run_em_refinement()` (main loop)
- [ ] Add `run_joint_optimization()` to CitegeistModel
- [ ] Add unit tests for new module
- [ ] Update benchmark script with new flag
- [ ] Run Xenium benchmark comparison
- [ ] Document results and update CLAUDE.md if successful
