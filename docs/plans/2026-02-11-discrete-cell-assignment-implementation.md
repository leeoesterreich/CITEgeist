# Discrete Cell Assignment Implementation Plan for Module 3

**Date:** 2026-02-11
**Branch:** `feat/discrete-cell-assignment`
**Parent Branch:** `feat/cellpose-nuclei-prior`
**Design Document:** `docs/plans/2026-02-11-discrete-cell-assignment-design.md`
**Status:** Implementation Plan

---

## Overview

This plan details the implementation of discrete cell assignment for CITEgeist Module 3, replacing continuous proportion estimation with integer cell count assignment based on nuclei detection from Cellpose segmentation.

### Key Changes from Current Implementation

| Aspect | Current (Continuous) | New (Discrete) |
|--------|---------------------|----------------|
| **Variables** | Y[i,t] ∈ [0,1] (proportions) | c[i,t] ∈ Z≥0 (integer counts) |
| **Constraint** | Σ Y ∈ [0.9, 1.2] | Σ c = N_i (exact nuclei count) |
| **Preprocessing** | CLR normalization | Winsorization only, no per-spot norm |
| **Optimization** | QP (quadratic program) | IQP (integer quadratic program) |
| **Phase 2** | Proportions as soft weights | Fixed counts as hard multipliers |

---

## Task Breakdown

### Phase 0: Prerequisites

#### Task 0.1: Verify Cellpose Integration Branch
**Dependencies:** None
**Priority:** Prerequisite

The discrete cell assignment feature depends on Cellpose segmentation from `feat/cellpose-nuclei-prior`. The parent branch has `segmentation.py` (uncommitted) providing:

```python
def compute_spot_nuclei_counts_from_adata(adata, resolution_mode='hires', ...) -> SegmentationResult
```

**Acceptance Criteria:**
- [ ] Parent branch has working `segmentation.py`
- [ ] Nuclei counts extractable per spot as pd.Series

---

### Phase 1: Preprocessing Changes

#### Task 1.1: Create Alternative Antibody Preprocessing
**Dependencies:** None
**Priority:** High
**File:** `CITEgeist/model/citegeist_model.py`

Add `preprocess_antibody_discrete()` that:
1. Winsorizes extreme values (5th-95th percentile)
2. Optional per-marker column scaling to [0, 1]
3. Does NOT apply CLR or per-spot normalization
4. Preserves cellularity signal (more cells = more signal)

```python
def preprocess_antibody_discrete(
    self,
    winsorize_lower: int = 5,
    winsorize_upper: int = 95,
    scale_per_marker: bool = True,
) -> None:
    """
    Preprocess antibody data for discrete cell assignment.

    Unlike preprocess_antibody(), does NOT apply per-spot normalization,
    preserving the cellularity signal.
    """
```

**Acceptance Criteria:**
- [ ] Method exists and runs without error
- [ ] Row sums vary (no per-spot normalization)
- [ ] Per-marker max is 1.0 when scale_per_marker=True

---

### Phase 2: IQP Solver Implementation

#### Task 2.1: Create Integer Quadratic Program Solver
**Dependencies:** Task 0.1
**Priority:** Critical
**File:** `CITEgeist/model/gurobi_impl.py`

```python
def solve_discrete_cell_counts(
    marker_level_data: np.ndarray,      # (N, M) antibody data
    marker_names: List[str],
    assignment_matrix: np.ndarray,       # (M, T) marker-to-type mapping
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,           # (N,) integer nuclei per spot
    beta_values: np.ndarray,             # (M,) per-marker scaling
    alpha_values: Optional[np.ndarray] = None,  # (M,) baselines
    max_nuclei_cap: int = 30,            # Above this, use relaxation
    timeout_per_spot: float = 60.0,
) -> np.ndarray:
    """
    Solve IQP for discrete cell counts given fixed beta (E-step).

    Returns:
        c_values: (N, T) integer cell counts per type per spot
    """
```

**Mathematical Formulation:**
```
minimize    Σᵢ Σₘ (X[i,m] - Σₜ c[i,t] × profile[t,m] × β[m])²
subject to  Σₜ c[i,t] = N_i     ∀ spots i with N_i > 0
            c[i,t] ∈ Z≥0        ∀ i, t
```

**Implementation Notes:**
- Use `GRB.INTEGER` for c variables when N_i ≤ max_nuclei_cap
- Use continuous relaxation + rounding for N_i > max_nuclei_cap
- Skip spots with N_i = 0 (return zeros)
- Handle infeasible/timeout with uniform fallback

**Acceptance Criteria:**
- [ ] Returns integer counts
- [ ] Sum equals nuclei count per spot
- [ ] Spots with N=0 get all zeros
- [ ] High-N spots use relaxation

---

#### Task 2.2: Implement EM Algorithm for Beta
**Dependencies:** Task 2.1
**Priority:** Critical
**File:** `CITEgeist/model/gurobi_impl.py`

```python
def optimize_discrete_cell_assignment_em(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,
    max_em_iterations: int = 20,
    beta_convergence_tol: float = 1e-3,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    max_nuclei_cap: int = 30,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    """
    EM algorithm for discrete cell assignment with per-marker beta.

    E-step: Solve IQP for cell counts given beta
    M-step: Update beta via OLS given cell counts

    Returns:
        c_values: (N, T) integer cell counts
        beta_values: (M,) final per-marker scaling
        marker_beta_dict: {marker_name: beta_value}
    """
```

**M-step Update:**
```
β[m] = Σᵢ X[i,m] × pred[i,m] / Σᵢ pred[i,m]²
where pred[i,m] = Σₜ c[i,t] × profile[t,m]
```

**Acceptance Criteria:**
- [ ] Converges within max_iterations
- [ ] Beta within [beta_min, beta_max]
- [ ] Loss decreases monotonically

---

### Phase 3: Model Integration

#### Task 3.1: Add run_discrete_cell_assignment Method
**Dependencies:** Tasks 2.2, 1.1
**Priority:** High
**File:** `CITEgeist/model/citegeist_model.py`

```python
def run_discrete_cell_assignment(
    self,
    nuclei_counts: Optional[pd.Series] = None,
    max_em_iterations: int = 20,
    beta_convergence_tol: float = 1e-3,
    max_nuclei_cap: int = 30,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
) -> pd.DataFrame:
    """
    Phase 1: Assign discrete cell identities using IQP with EM.

    Returns:
        DataFrame with cell type columns, integer count values.
    """
```

**Stores in self.results:**
- `"marker_beta"`: Dict of per-marker betas
- `"discrete_cell_counts"`: (N, T) array
- `"nuclei_counts"`: (N,) array
- `"cell_prop"`: Derived proportions for Phase 2 compatibility

**Acceptance Criteria:**
- [ ] Runs without error on valid input
- [ ] Output has integer values
- [ ] Saves CSV to output_folder
- [ ] Sets self.results for Phase 2

---

#### Task 3.2: Modify run_cell_expression_pass1 for Fixed Counts
**Dependencies:** Task 3.1
**Priority:** High
**File:** `CITEgeist/model/citegeist_model.py`

Add parameters:
```python
def run_cell_expression_pass1(
    self,
    # ... existing params ...
    cell_counts: Optional[pd.DataFrame] = None,  # NEW
    use_discrete_mode: bool = False,             # NEW
)
```

**Changes:**
1. If discrete mode: use integer counts as multipliers
2. Skip spots with 0 total cells
3. GEX model becomes: `G[i,g] = Σₜ c[i,t] × E[i,t,g]`

**Acceptance Criteria:**
- [ ] Accepts cell_counts parameter
- [ ] Discrete mode uses integer counts
- [ ] Zero-cell spots handled gracefully

---

### Phase 4: Testing

#### Task 4.1: Unit Tests for IQP Solver
**Priority:** High
**File:** `tests/test_discrete_assignment.py` (NEW)

Test cases:
- Known ground truth recovery
- Zero nuclei handling
- Single cell assignment
- High nuclei relaxation

---

#### Task 4.2: Unit Tests for EM Algorithm
**Priority:** High
**File:** `tests/test_discrete_assignment.py`

Test cases:
- Convergence within iterations
- Beta recovery from synthetic data
- Monotonic loss decrease

---

#### Task 4.3: Integration Tests
**Priority:** High
**File:** `tests/test_discrete_assignment.py`

Test cases:
- End-to-end simulated data pipeline
- Discrete vs continuous comparison

---

### Phase 5: Documentation

#### Task 5.1: Example Script
**Priority:** Medium
**File:** `examples/run_discrete_assignment.py` (NEW)

Demonstrates full pipeline:
1. Cellpose segmentation
2. preprocess_antibody_discrete()
3. run_discrete_cell_assignment()
4. run_cell_expression_pass1(use_discrete_mode=True)

---

#### Task 5.2: Update CLAUDE.md
**Priority:** Medium
**File:** `CLAUDE.md`

Document new methods and discrete mode.

---

## Dependency Graph

```
Task 0.1 (Cellpose)
    │
    ├──► Task 1.1 (Preprocessing) ──► Task 3.1 (Model) ──► Task 3.2 (Phase 2)
    │                                      │                    │
    └──► Task 2.1 (IQP) ──► Task 2.2 (EM) ─┘                    │
              │                                                  │
              ▼                                                  ▼
         Task 4.1 (Tests)                              Task 4.3 (Integration)
              │                                                  │
              ▼                                                  ▼
         Task 4.2 (Tests)                              Task 5.1 (Example)
                                                                 │
                                                                 ▼
                                                       Task 5.2 (Docs)
```

---

## Implementation Order

### Week 1: Core Solver
1. Task 0.1: Verify Cellpose (ensure parent branch ready)
2. Task 1.1: Preprocessing changes
3. Task 2.1: IQP solver

### Week 2: EM and Integration
4. Task 2.2: EM algorithm
5. Task 3.1: Model integration
6. Task 3.2: Phase 2 integration

### Week 3: Testing and Documentation
7. Tasks 4.1-4.3: All testing
8. Tasks 5.1-5.2: Documentation

---

## Acceptance Criteria Summary

### Functional
- [ ] Discrete counts sum to nuclei counts per spot
- [ ] EM converges within 20 iterations
- [ ] Phase 2 works with fixed counts
- [ ] Zero-nuclei spots handled

### Performance
- [ ] Phase 1 completes in < 2x continuous time
- [ ] Per-spot IQP < 60 seconds

### Quality
- [ ] All tests pass
- [ ] PEP 8 compliant
- [ ] Docstrings complete

---

## Critical Files

| File | Changes |
|------|---------|
| `CITEgeist/model/gurobi_impl.py` | Add `solve_discrete_cell_counts()`, `optimize_discrete_cell_assignment_em()` |
| `CITEgeist/model/citegeist_model.py` | Add `preprocess_antibody_discrete()`, `run_discrete_cell_assignment()`, modify `run_cell_expression_pass1()` |
| `tests/test_discrete_assignment.py` | NEW: All unit and integration tests |
| `examples/run_discrete_assignment.py` | NEW: Example script |
