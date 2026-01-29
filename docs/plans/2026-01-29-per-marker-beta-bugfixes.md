# Per-Marker Beta Optimization Bug Fixes

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix three bugs in `optimize_cell_proportions_per_marker()` that cause uniform proportion spreading: marker-count-weighted loss asymmetry, exclusive shared-marker assignment via `argmax`, and insufficient beta range for cross-marker signal compensation.

**Architecture:** All three bugs live in `gurobi_impl.py`. The loss function currently sums raw error terms per marker, giving 2-marker cell types (Fibroblasts, Macrophages) double the loss weight vs 1-marker types (Epithelial, Endothelial). Shared markers (CD3E in both CD4+ T and CD8+ T) are assigned exclusively to whichever cell type appears first via `argmax`, starving the second type. Beta clipping at `[0.1, 2.0]` can't compensate for 42x signal dynamic range (PanCK mean=88 vs Vimentin mean=3728). These bugs affect the global optimizer, the local finetuning solver, and the marker mapping function.

**Tech Stack:** Python 3.10, gurobipy 11.0.2, numpy 1.26.4, pytest

**Evidence (from debugging session):**
- PanCK (Epithelial): mean=88, 55% spots >0, ratio-to-VIM=1:42
- Vimentin (Fibroblasts): mean=3728, near-universal expression
- Benchmark Epithelial Pearson r = 0.478, Fibroblasts r = 0.360
- All predicted proportions cluster at 0.09-0.15 (uniform spreading)

---

## Bug Inventory

| # | Bug | Location | Root Cause | Fix |
|---|-----|----------|------------|-----|
| **1** | Marker-count loss asymmetry | `gurobi_impl.py:625-637` (global), `:1298-1306` (local) | Loss sums per-marker terms without normalizing by markers-per-celltype. 2-marker types get 2× gradient. | Divide each marker's error contribution by the number of markers assigned to its owner cell type. |
| **2** | Exclusive shared-marker assignment | `gurobi_impl.py:592` (global), `:1236` (local) | `owners = np.argmax(assignment_matrix, axis=1)` picks first cell type for shared markers (CD3E). Second type loses that marker entirely. | Replace `argmax` with multi-owner support: iterate over all owners for shared markers, contributing fractional weight to each. |
| **3** | Beta range too narrow | `gurobi_impl.py:689` (global), `:1337` (local), defaults in `citegeist_model.py:409-410` | `beta_min=0.1, beta_max=2.0` gives 20× range; actual signal dynamic range is 42×+. | Widen default to `beta_min=0.01, beta_max=100.0`. |

---

### Task 1: Write Failing Tests for All Three Bugs

**Files:**
- Create: `tests/test_per_marker_beta_fixes.py`

**Step 1: Write the test file**

```python
"""
Tests for per-marker beta bug fixes.

Bug 1: Marker-count loss asymmetry
Bug 2: Exclusive shared-marker assignment
Bug 3: Beta range too narrow for signal dynamic range
"""

import numpy as np
import pytest

from CITEgeist.model.gurobi_impl import (
    map_antibodies_to_profiles_v2,
    optimize_cell_proportions_per_marker,
)


def make_synthetic_marker_data(N=100, seed=42):
    """
    Create synthetic marker data with known ground truth.

    Layout (N spots):
      - Spots 0..49: 80% TypeA (markers: M1+M2), 20% TypeB (marker: M3)
      - Spots 50..99: 20% TypeA, 80% TypeB

    TypeA has 2 markers (M1, M2), TypeB has 1 marker (M3).
    If the optimizer is fair, TypeB should still get ~80% in spots 50..99.
    """
    rng = np.random.default_rng(seed)
    N_half = N // 2

    # Ground truth proportions
    gt = np.zeros((N, 2))
    gt[:N_half, 0] = 0.8   # TypeA dominant
    gt[:N_half, 1] = 0.2
    gt[N_half:, 0] = 0.2   # TypeB dominant
    gt[N_half:, 1] = 0.8

    # Generate marker signals: signal = proportion * scale + noise
    # M1, M2 belong to TypeA; M3 belongs to TypeB
    # All markers have similar scale (no signal bias)
    scale = 1.0
    noise = 0.05
    M1 = gt[:, 0] * scale + rng.normal(0, noise, N)
    M2 = gt[:, 0] * scale + rng.normal(0, noise, N)
    M3 = gt[:, 1] * scale + rng.normal(0, noise, N)

    marker_data = np.column_stack([M1, M2, M3]).clip(0, None)
    # Normalize per column
    col_max = marker_data.max(axis=0)
    col_max[col_max == 0] = 1e-6
    marker_data = marker_data / col_max

    marker_names = ["M1", "M2", "M3"]
    cell_type_names = ["TypeA", "TypeB"]
    # Assignment: M1->TypeA, M2->TypeA, M3->TypeB
    assignment_matrix = np.array([
        [1.0, 0.0],  # M1 -> TypeA
        [1.0, 0.0],  # M2 -> TypeA
        [0.0, 1.0],  # M3 -> TypeB
    ])

    return marker_data, marker_names, assignment_matrix, cell_type_names, gt


def make_shared_marker_data(N=100, seed=42):
    """
    Create synthetic data where one marker (Mshared) belongs to two cell types.

    Layout (N spots):
      - Spots 0..33: 80% TypeA (markers: Mshared + Ma)
      - Spots 34..66: 80% TypeB (markers: Mshared + Mb)
      - Spots 67..99: 50/50 TypeA/TypeB

    Both TypeA and TypeB use Mshared. If argmax assigns Mshared exclusively
    to TypeA, TypeB loses half its signal.
    """
    rng = np.random.default_rng(seed)
    N3 = N // 3

    gt = np.zeros((N, 2))
    gt[:N3, 0] = 0.8         # TypeA dominant
    gt[:N3, 1] = 0.2
    gt[N3:2*N3, 0] = 0.2     # TypeB dominant
    gt[N3:2*N3, 1] = 0.8
    gt[2*N3:, 0] = 0.5       # Mixed
    gt[2*N3:, 1] = 0.5

    scale = 1.0
    noise = 0.05
    # Mshared reflects sum of both types
    Mshared = (gt[:, 0] + gt[:, 1]) * scale * 0.5 + rng.normal(0, noise, N)
    Ma = gt[:, 0] * scale + rng.normal(0, noise, N)
    Mb = gt[:, 1] * scale + rng.normal(0, noise, N)

    marker_data = np.column_stack([Mshared, Ma, Mb]).clip(0, None)
    col_max = marker_data.max(axis=0)
    col_max[col_max == 0] = 1e-6
    marker_data = marker_data / col_max

    marker_names = ["Mshared", "Ma", "Mb"]
    cell_type_names = ["TypeA", "TypeB"]
    # Mshared belongs to BOTH types
    assignment_matrix = np.array([
        [1.0, 1.0],  # Mshared -> TypeA AND TypeB
        [1.0, 0.0],  # Ma -> TypeA only
        [0.0, 1.0],  # Mb -> TypeB only
    ])

    return marker_data, marker_names, assignment_matrix, cell_type_names, gt


class TestMarkerCountAsymmetry:
    """Bug 1: Cell types with more markers shouldn't get more weight."""

    def test_single_marker_type_not_underestimated(self):
        """TypeB (1 marker) should reach ~0.8 in its dominant spots,
        not be suppressed by TypeA (2 markers) getting double loss weight."""
        marker_data, names, assign, ct_names, gt = make_synthetic_marker_data()

        Y, beta, beta_dict = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=names,
            assignment_matrix=assign,
            cell_type_names=ct_names,
            lambda_reg=0.1,
            alpha=0.5,
            max_iterations=10,
            warn_only=True,
            lambda_laplacian=0.0,
        )

        # In TypeB-dominant spots (50..99), TypeB should be >= 0.5
        typeb_in_dominant = Y[50:, 1].mean()
        typea_in_dominant = Y[50:, 0].mean()
        assert typeb_in_dominant > 0.5, (
            f"TypeB (1 marker) mean={typeb_in_dominant:.3f} in its dominant region. "
            f"Should be >0.5 but marker-count asymmetry is suppressing it."
        )
        assert typeb_in_dominant > typea_in_dominant, (
            f"TypeB={typeb_in_dominant:.3f} should exceed TypeA={typea_in_dominant:.3f} "
            f"in TypeB-dominant spots"
        )


class TestSharedMarkerAssignment:
    """Bug 2: Shared markers should contribute to all owner cell types."""

    def test_shared_marker_not_exclusive(self):
        """When Mshared belongs to both TypeA and TypeB, both types should
        benefit from it. argmax assigns it to only one."""
        marker_data, names, assign, ct_names, gt = make_shared_marker_data()

        Y, beta, beta_dict = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=names,
            assignment_matrix=assign,
            cell_type_names=ct_names,
            lambda_reg=0.1,
            alpha=0.5,
            max_iterations=10,
            warn_only=True,
            lambda_laplacian=0.0,
        )

        # In TypeB-dominant region (34..66), TypeB should be dominant
        typeb_dominant = Y[34:67, 1].mean()
        typea_dominant = Y[34:67, 0].mean()
        assert typeb_dominant > typea_dominant, (
            f"In TypeB-dominant region: TypeB={typeb_dominant:.3f}, TypeA={typea_dominant:.3f}. "
            f"TypeB should be higher but shared marker exclusivity is breaking it."
        )

        # Both types should be reasonable in mixed region
        mixed_a = Y[67:, 0].mean()
        mixed_b = Y[67:, 1].mean()
        assert abs(mixed_a - mixed_b) < 0.2, (
            f"Mixed region: TypeA={mixed_a:.3f}, TypeB={mixed_b:.3f}. "
            f"Should be roughly equal but shared marker bias creates asymmetry."
        )


class TestBetaRangeCompensation:
    """Bug 3: Beta range must accommodate signal dynamic range."""

    def test_weak_signal_marker_compensated(self):
        """A marker with 50x weaker signal should still drive proportions
        correctly if beta can scale up enough."""
        rng = np.random.default_rng(42)
        N = 100

        gt = np.zeros((N, 2))
        gt[:50, 0] = 0.8   # TypeA dominant
        gt[:50, 1] = 0.2
        gt[50:, 0] = 0.2   # TypeB dominant
        gt[50:, 1] = 0.8

        # M1 (TypeA): STRONG signal (like Vimentin)
        # M2 (TypeB): WEAK signal (like PanCK) - 50x weaker
        M1 = gt[:, 0] * 50.0 + rng.normal(0, 1, N)
        M2 = gt[:, 1] * 1.0 + rng.normal(0, 0.05, N)

        marker_data = np.column_stack([M1, M2]).clip(0, None)
        col_max = marker_data.max(axis=0)
        col_max[col_max == 0] = 1e-6
        marker_data = marker_data / col_max

        marker_names = ["Mstrong", "Mweak"]
        cell_type_names = ["TypeA", "TypeB"]
        assignment_matrix = np.array([
            [1.0, 0.0],  # Mstrong -> TypeA
            [0.0, 1.0],  # Mweak -> TypeB
        ])

        Y, beta, beta_dict = optimize_cell_proportions_per_marker(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
            lambda_reg=0.1,
            alpha=0.5,
            max_iterations=10,
            warn_only=True,
            lambda_laplacian=0.0,
        )

        # TypeB should be dominant in spots 50..99 even with weak signal
        typeb_dominant = Y[50:, 1].mean()
        assert typeb_dominant > 0.5, (
            f"TypeB (weak signal) mean={typeb_dominant:.3f} in dominant region. "
            f"Beta range [{beta_dict['Mweak']:.2f}] can't compensate for signal weakness."
        )
```

**Step 2: Run tests to verify they fail**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_per_marker_beta_fixes.py -v`

Expected: At least one test should FAIL (validating the bugs exist). If all pass, the synthetic data needs adjusting to trigger the bugs — but given the code analysis, the marker-count and shared-marker tests should fail.

**Step 3: Commit failing tests**

```bash
git add tests/test_per_marker_beta_fixes.py
git commit -m "test: add failing tests for per-marker beta optimization bugs

Three bugs identified in optimize_cell_proportions_per_marker():
1. Marker-count loss asymmetry (2-marker types get 2x loss weight)
2. Exclusive shared-marker assignment via argmax
3. Beta range [0.1, 2.0] too narrow for 42x signal dynamic range"
```

---

### Task 2: Fix Bug 1 — Marker-Count Loss Normalization

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:590-637` (global optimizer)
- Modify: `CITEgeist/model/gurobi_impl.py:1235-1306` (local finetuner)

**The problem:** The loss function at line 625-637 sums `(S[i,m] - β[m] * Y[i,j])²` for every marker `m`. A cell type with 2 Major markers contributes 2N error terms; a cell type with 1 marker contributes N terms. The optimizer naturally prioritizes reducing error for multi-marker types.

**The fix:** Precompute `markers_per_celltype[j]` = number of markers assigned to cell type `j`. Weight each marker's error terms by `1.0 / markers_per_celltype[owner(m)]`.

**Step 1: Add marker-count weight computation (global optimizer)**

After line 593 (`marker_has_owner = ...`), add:

```python
# Compute markers-per-celltype for loss normalization (Bug 1 fix)
# Cell types with more markers should not get more total loss weight
markers_per_celltype = np.zeros(T, dtype=np.float64)
for m in range(M):
    if marker_has_owner[m]:
        for j in range(T):
            if assignment_matrix[m, j] > 0:
                markers_per_celltype[j] += 1
# Avoid division by zero for cell types with no markers (Unknown)
markers_per_celltype = np.maximum(markers_per_celltype, 1.0)
```

**Step 2: Apply weight in loss construction (global optimizer)**

Replace lines 625-637:

```python
# Objective: sum_m sum_i w[m] * (S[i,m] - beta[m] * Y[i, owner(m)])^2
# where w[m] = 1 / markers_per_celltype[owner(m)]  (marker-count normalization)
error_terms = []
for m in range(M):
    if not marker_has_owner[m]:
        continue

    j = owners[m]
    beta_m = beta_values[m]
    weight = 1.0 / markers_per_celltype[j]

    for i in range(N):
        S_im = marker_level_data[i, m]
        Y_ij = Y[i, j]
        error_terms.append(weight * (S_im - beta_m * Y_ij) * (S_im - beta_m * Y_ij))
```

**Step 3: Apply same fix in local finetuner**

In `deconvolute_local_cell_proportions_per_marker()`, after line 1237 (`marker_has_owner = ...`), add the same `markers_per_celltype` computation.

Then replace lines 1298-1306:

```python
error_terms = []
for m in range(M):
    if not marker_has_owner[m]:
        continue
    j = owners[m]
    beta_m = local_beta[m]
    weight = 1.0 / markers_per_celltype[j]
    for i in range(local_N):
        S_im = local_marker_data[i, m]
        error_terms.append(weight * (S_im - beta_m * Y_vars[i, j]) ** 2)
```

**Step 4: Run tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_per_marker_beta_fixes.py::TestMarkerCountAsymmetry -v`

Expected: `test_single_marker_type_not_underestimated` should now PASS.

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "fix: normalize per-marker loss by markers-per-celltype

Cell types with N markers were getting N× the loss weight, causing the
optimizer to prioritize multi-marker types (Fibroblasts, Macrophages)
over single-marker types (Epithelial, Endothelial). Now each marker's
error contribution is divided by the number of markers for its cell type,
giving equal total loss weight regardless of marker count.

Applied to both global optimizer and local finetuner."
```

---

### Task 3: Fix Bug 2 — Multi-Owner Shared Marker Support

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:590-637` (global optimizer)
- Modify: `CITEgeist/model/gurobi_impl.py:1235-1306` (local finetuner)
- Modify: `CITEgeist/model/gurobi_impl.py:676-689` (beta update)
- Modify: `CITEgeist/model/gurobi_impl.py:1326-1337` (local beta update)

**The problem:** Line 592 uses `owners = np.argmax(assignment_matrix, axis=1)` which picks ONE cell type per marker. When CD3E is Major for both CD4+ T and CD8+ T, `argmax` picks whichever is first. The second type loses CD3E entirely.

**The fix:** Replace the exclusive `owners` array with iteration over all assignments in the assignment matrix. For shared markers, contribute error terms to ALL owner cell types (weighted by `1/n_owners` so total contribution stays the same).

**Step 1: Replace exclusive owner with multi-owner iteration (global optimizer)**

Replace the `owners` computation and loss construction (lines 590-637). The new approach iterates over `assignment_matrix` entries directly:

```python
# Precompute marker-to-celltype assignments (supports shared markers)
# marker_owners[m] = list of (celltype_idx, weight) tuples
marker_owners = []
for m in range(M):
    owners_for_m = []
    for j in range(T):
        if assignment_matrix[m, j] > 0:
            owners_for_m.append(j)
    marker_owners.append(owners_for_m)

marker_has_owner = np.array([len(o) > 0 for o in marker_owners])

# Compute markers-per-celltype for loss normalization (Bug 1 fix)
markers_per_celltype = np.zeros(T, dtype=np.float64)
for m in range(M):
    for j in marker_owners[m]:
        markers_per_celltype[j] += 1
markers_per_celltype = np.maximum(markers_per_celltype, 1.0)
```

Then replace the loss construction:

```python
# Objective: for each marker m and each owner j of m,
# add (1/n_owners) * (1/markers_per_celltype[j]) * (S[i,m] - beta[m,j] * Y[i,j])^2
error_terms = []
for m in range(M):
    if not marker_has_owner[m]:
        continue

    owners_m = marker_owners[m]
    n_owners = len(owners_m)
    beta_m = beta_values[m]

    for j in owners_m:
        weight = 1.0 / (n_owners * markers_per_celltype[j])
        for i in range(N):
            S_im = marker_level_data[i, m]
            Y_ij = Y[i, j]
            error_terms.append(weight * (S_im - beta_m * Y_ij) * (S_im - beta_m * Y_ij))
```

**Step 2: Update beta computation for multi-owner (global optimizer)**

Replace lines 676-689. For shared markers, beta is updated using the sum of all owner proportions:

```python
# Update beta (per-marker closed-form solution)
beta_new = np.zeros(M, dtype=np.float64)
for m in range(M):
    if not marker_has_owner[m]:
        beta_new[m] = 1.0
        continue

    # For shared markers, use sum of owner proportions
    owners_m = marker_owners[m]
    Y_combined = np.zeros(N, dtype=np.float64)
    for j in owners_m:
        Y_combined += Y_values[:, j]

    S_m = marker_level_data[:, m]
    denominator = np.dot(Y_combined, Y_combined) + 1e-9
    beta_new[m] = np.dot(S_m, Y_combined) / denominator
    beta_new[m] = np.clip(beta_new[m], beta_min, beta_max)
```

**Step 3: Apply same changes to local finetuner**

In `deconvolute_local_cell_proportions_per_marker()`:

Replace line 1236 (`owners = np.argmax(...)`) with the same `marker_owners` list construction.

Replace lines 1298-1306 with the multi-owner loss construction (same pattern as Step 1 but using `local_marker_data` and `Y_vars`).

Replace lines 1326-1337 beta update with the multi-owner beta update (same pattern as Step 2 but using local data).

**Step 4: Run tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_per_marker_beta_fixes.py::TestSharedMarkerAssignment -v`

Expected: `test_shared_marker_not_exclusive` should PASS.

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py
git commit -m "fix: support shared markers across multiple cell types

Replaced argmax-based exclusive marker assignment with multi-owner
iteration. Shared markers (e.g. CD3E in both CD4+ T and CD8+ T) now
contribute error terms to ALL owner cell types, weighted by 1/n_owners.
Beta learning uses combined owner proportions for shared markers.

Applied to both global optimizer, local finetuner, and beta updates."
```

---

### Task 4: Fix Bug 3 — Widen Beta Range Defaults

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:538` (global function signature)
- Modify: `CITEgeist/model/gurobi_impl.py:1200` (local function signature)
- Modify: `CITEgeist/model/citegeist_model.py:409-410` (orchestrator defaults)

**The problem:** Default `beta_min=0.1, beta_max=2.0` gives only 20× range. PanCK-to-Vimentin signal ratio is 42×. After column-max normalization the ratio shrinks, but beta still needs more room to compensate for markers with very different normalized distributions (sparse vs uniform).

**Step 1: Update defaults in all three locations**

In `optimize_cell_proportions_per_marker()` signature (line 538):
```python
beta_min: float = 0.01,
beta_max: float = 100.0,
```

In `deconvolute_local_cell_proportions_per_marker()` signature (line 1200):
```python
beta_min: float = 0.01,
beta_max: float = 100.0,
```

In `citegeist_model.py` `run_cell_proportion_model()` signature (lines 409-410):
```python
beta_min=0.01,
beta_max=100.0,
```

**Step 2: Run tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_per_marker_beta_fixes.py::TestBetaRangeCompensation -v`

Expected: `test_weak_signal_marker_compensated` should PASS.

**Step 3: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/model/citegeist_model.py
git commit -m "fix: widen per-marker beta range to [0.01, 100.0]

Previous range [0.1, 2.0] provided only 20x compensation range,
insufficient for antibody panels with 42x+ signal dynamic range
(PanCK mean=88 vs Vimentin mean=3728). New range [0.01, 100.0]
allows beta to fully compensate for signal strength differences
after column-max normalization."
```

---

### Task 5: Run Full Test Suite and Benchmark Validation

**Files:**
- Test: `tests/test_per_marker_beta_fixes.py`
- Test: `tests/test_citegeist_simulated.py`

**Step 1: Run all new tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_per_marker_beta_fixes.py -v`

Expected: All 3 tests PASS.

**Step 2: Run existing simulated data tests (regression check)**

Run: Submit via SLURM (requires Gurobi license):
```bash
sbatch tests/sbatch_module12_benchmark.sh  # or appropriate existing test script
```

Expected: No regressions in existing tests.

**Step 3: Rerun Xenium benchmark with fixes**

Use the existing benchmark scripts. Submit:
```bash
cd Benchmarking/xenium_benchmarking/CITEgeist/slurm
sbatch run_benchmark.sh  # default manual mode runs achievable-7
```

Then evaluate:
```bash
cd Benchmarking/xenium_benchmarking/evaluation/slurm
sbatch evaluate_achievable_7.sh
```

Compare Epithelial and Fibroblast Pearson r against previous baseline:
- Baseline Epithelial r = 0.478
- Baseline Fibroblasts r = 0.360
- Target: meaningful improvement in both, especially Epithelial

**Step 4: Commit final state**

```bash
git add -A
git commit -m "test: validate per-marker beta fixes on simulated and benchmark data"
```

---

## Affected Functions Summary

All changes are in `CITEgeist/model/gurobi_impl.py` except the default parameter change in `citegeist_model.py`:

| Function | Lines | Bug 1 | Bug 2 | Bug 3 |
|----------|-------|-------|-------|-------|
| `optimize_cell_proportions_per_marker()` | 527-729 | ✓ loss weight | ✓ multi-owner | ✓ defaults |
| `deconvolute_local_cell_proportions_per_marker()` | 1183-1378 | ✓ loss weight | ✓ multi-owner | ✓ defaults |
| `finetune_cell_proportions_per_marker()` | 1386-1537 | — (delegates) | — (delegates) | — |
| `map_antibodies_to_profiles_v2()` | 266-342 | — | — (already stores multi-assignment) | — |
| `run_cell_proportion_model()` (citegeist_model.py) | 388-577 | — | — | ✓ defaults |

Note: `map_antibodies_to_profiles_v2()` already builds the assignment matrix with multi-assignment (line 338: `assignment_matrix[m, j] = 1.0` for each cell type). The bug is only in how the optimizer **consumes** it.
