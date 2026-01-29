# Marker Exclusivity Weighting Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add reconstruction exclusivity weighting to the finetuning pass so non-specific markers (e.g., VIM) contribute less to the loss function, improving cell type discrimination.

**Architecture:** After the global EM pass converges, compute per-marker exclusivity scores by correlating each marker's signal with its owner cell type vs. other cell types. Pass these scores into the finetuning pass as loss weights. No changes to the global EM pass.

**Tech Stack:** numpy (correlation computation), gurobipy (loss weighting in finetuning QP)

**Design doc:** `docs/plans/2026-01-29-marker-exclusivity-weighting-design.md`

---

### Task 1: Add `compute_marker_exclusivity()` function with tests

**Files:**
- Create: `tests/test_marker_exclusivity.py`
- Modify: `CITEgeist/model/gurobi_impl.py` (add new function after line 757)

**Step 1: Write the failing test**

Create `tests/test_marker_exclusivity.py`:

```python
"""Tests for marker exclusivity weighting."""
import numpy as np
import pytest

from CITEgeist.model.gurobi_impl import compute_marker_exclusivity


class TestComputeMarkerExclusivity:
    """Test exclusivity score computation."""

    def _make_data(self, N=200, seed=42):
        """Create synthetic marker data with known specificity patterns."""
        rng = np.random.RandomState(seed)

        # 3 cell types: A, B, C
        # Ground truth proportions
        Y = np.zeros((N, 3))
        Y[:70, 0] = 0.8   # Region 1: mostly A
        Y[:70, 1] = 0.1
        Y[:70, 2] = 0.1
        Y[70:140, 0] = 0.1
        Y[70:140, 1] = 0.8  # Region 2: mostly B
        Y[70:140, 2] = 0.1
        Y[140:, 0] = 0.1
        Y[140:, 1] = 0.1
        Y[140:, 2] = 0.8   # Region 3: mostly C

        # Marker 0: specific to A (tracks Y[:,0] closely)
        # Marker 1: specific to B (tracks Y[:,1] closely)
        # Marker 2: ubiquitous (high everywhere, assigned to C)
        # Marker 3: specific to C (tracks Y[:,2] closely)
        S = np.zeros((N, 4))
        S[:, 0] = Y[:, 0] * 0.9 + rng.normal(0, 0.05, N)  # specific to A
        S[:, 1] = Y[:, 1] * 0.9 + rng.normal(0, 0.05, N)  # specific to B
        S[:, 2] = 0.5 + rng.normal(0, 0.05, N)             # ubiquitous, assigned to C
        S[:, 3] = Y[:, 2] * 0.9 + rng.normal(0, 0.05, N)  # specific to C
        S = np.clip(S, 0, 1)

        # Assignment: marker 0->A, 1->B, 2->C, 3->C
        assignment = np.zeros((4, 3))
        assignment[0, 0] = 1  # marker 0 -> cell type A
        assignment[1, 1] = 1  # marker 1 -> cell type B
        assignment[2, 2] = 1  # marker 2 -> cell type C (ubiquitous)
        assignment[3, 2] = 1  # marker 3 -> cell type C (specific)

        marker_owners = [[0], [1], [2], [2]]

        return S, Y, assignment, marker_owners

    def test_specific_markers_score_higher(self):
        """Specific markers should have higher exclusivity than ubiquitous ones."""
        S, Y, assignment, marker_owners = self._make_data()
        exclusivity = compute_marker_exclusivity(S, Y, marker_owners, assignment)

        assert exclusivity.shape == (4,)
        # Specific markers (0, 1, 3) should score higher than ubiquitous (2)
        assert exclusivity[0] > exclusivity[2], f"Specific marker 0 ({exclusivity[0]:.3f}) should > ubiquitous ({exclusivity[2]:.3f})"
        assert exclusivity[1] > exclusivity[2], f"Specific marker 1 ({exclusivity[1]:.3f}) should > ubiquitous ({exclusivity[2]:.3f})"
        assert exclusivity[3] > exclusivity[2], f"Specific marker 3 ({exclusivity[3]:.3f}) should > ubiquitous ({exclusivity[2]:.3f})"

    def test_exclusivity_range(self):
        """Exclusivity scores should be in [0.3, 1.0] (floored at 0.3)."""
        S, Y, assignment, marker_owners = self._make_data()
        exclusivity = compute_marker_exclusivity(S, Y, marker_owners, assignment)

        assert np.all(exclusivity >= 0.3), f"Min exclusivity {exclusivity.min():.3f} < 0.3 floor"
        assert np.all(exclusivity <= 1.0), f"Max exclusivity {exclusivity.max():.3f} > 1.0"

    def test_unowned_markers_get_default(self):
        """Markers with no owner should get exclusivity = 1.0 (neutral)."""
        S, Y, assignment, marker_owners = self._make_data()
        # Add an unowned marker
        S_extended = np.column_stack([S, np.random.rand(S.shape[0])])
        assignment_extended = np.vstack([assignment, np.zeros((1, 3))])
        marker_owners_extended = marker_owners + [[]]

        exclusivity = compute_marker_exclusivity(S_extended, Y, marker_owners_extended, assignment_extended)
        assert exclusivity[4] == 1.0, f"Unowned marker should get 1.0, got {exclusivity[4]:.3f}"

    def test_shared_marker_uses_combined_owners(self):
        """Shared markers should correlate with combined owner proportions."""
        rng = np.random.RandomState(42)
        N = 200
        Y = np.zeros((N, 3))
        Y[:100, 0] = 0.4
        Y[:100, 1] = 0.4
        Y[:100, 2] = 0.2
        Y[100:, 0] = 0.1
        Y[100:, 1] = 0.1
        Y[100:, 2] = 0.8

        # Marker 0 shared between A and B, tracks their sum
        S = np.zeros((N, 2))
        S[:, 0] = (Y[:, 0] + Y[:, 1]) * 0.5 + rng.normal(0, 0.05, N)
        S[:, 1] = Y[:, 2] * 0.9 + rng.normal(0, 0.05, N)
        S = np.clip(S, 0, 1)

        assignment = np.zeros((2, 3))
        assignment[0, 0] = 1  # shared: A
        assignment[0, 1] = 1  # shared: B
        assignment[1, 2] = 1  # exclusive: C

        marker_owners = [[0, 1], [2]]

        exclusivity = compute_marker_exclusivity(S, Y, marker_owners, assignment)
        # Shared marker should still get reasonable score since it tracks A+B vs C
        assert exclusivity[0] > 0.5, f"Shared marker should be discriminative, got {exclusivity[0]:.3f}"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_marker_exclusivity.py -v`
Expected: FAIL with `ImportError: cannot import name 'compute_marker_exclusivity'`

**Step 3: Write the implementation**

Add to `CITEgeist/model/gurobi_impl.py` after line 757 (after `optimize_cell_proportions_per_marker()` returns):

```python
def compute_marker_exclusivity(
    marker_level_data: np.ndarray,
    Y_values: np.ndarray,
    marker_owners: List[List[int]],
    assignment_matrix: np.ndarray,
    floor: float = 0.3,
    epsilon: float = 1e-9,
) -> np.ndarray:
    """
    Compute per-marker exclusivity scores measuring discriminative power.

    For each marker, measures how exclusively it correlates with its assigned
    cell type(s) versus the best non-owner cell type. Markers that track many
    cell types equally (e.g., VIM) get low scores; markers specific to their
    owner (e.g., CD68) get high scores.

    Args:
        marker_level_data: (N, M) normalized marker signals.
        Y_values: (N, T) cell type proportions from global EM pass.
        marker_owners: List of lists, marker_owners[m] = indices of owner cell types.
        assignment_matrix: (M, T) binary matrix mapping markers to cell types.
        floor: Minimum exclusivity score (default: 0.3).
        epsilon: Small constant to prevent division by zero.

    Returns:
        (M,) array of exclusivity scores in [floor, 1.0].
    """
    N, M = marker_level_data.shape
    T = Y_values.shape[1]
    exclusivity = np.ones(M, dtype=np.float64)

    for m in range(M):
        owners_m = marker_owners[m]
        if not owners_m:
            # Unowned markers: neutral weight (1.0)
            continue

        S_m = marker_level_data[:, m]

        # Owner correlation: correlate with combined owner proportions
        Y_owner = np.zeros(N, dtype=np.float64)
        for j in owners_m:
            Y_owner += Y_values[:, j]

        r_owner = np.corrcoef(S_m, Y_owner)[0, 1]
        if np.isnan(r_owner):
            r_owner = 0.0
        r_owner = max(r_owner, 0.0)

        # Best non-owner correlation
        owner_set = set(owners_m)
        r_best_other = 0.0
        for k in range(T):
            if k in owner_set:
                continue
            r_k = np.corrcoef(S_m, Y_values[:, k])[0, 1]
            if np.isnan(r_k):
                r_k = 0.0
            r_best_other = max(r_best_other, max(r_k, 0.0))

        # Exclusivity ratio
        denom = r_owner + r_best_other + epsilon
        exclusivity[m] = r_owner / denom

    # Apply floor
    exclusivity = np.clip(exclusivity, floor, 1.0)

    return exclusivity
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_marker_exclusivity.py -v`
Expected: All 5 tests PASS

**Step 5: Commit**

```bash
git add tests/test_marker_exclusivity.py CITEgeist/model/gurobi_impl.py
git commit -m "feat: add compute_marker_exclusivity() for discriminative marker weighting"
```

---

### Task 2: Wire exclusivity into finetuning loss function

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:1213-1231` (add parameter to `deconvolute_local_cell_proportions_per_marker`)
- Modify: `CITEgeist/model/gurobi_impl.py:1339-1351` (apply weight in loss construction)
- Modify: `CITEgeist/model/gurobi_impl.py:1429-1450` (add parameter to `finetune_cell_proportions_per_marker`)
- Modify: `CITEgeist/model/gurobi_impl.py:1516-1537` (pass through in executor.submit)

**Step 1: Write the failing test**

Add to `tests/test_marker_exclusivity.py`:

```python
class TestExclusivityInFinetuning:
    """Test that exclusivity weights are applied in finetuning loss."""

    def test_exclusivity_parameter_accepted(self):
        """finetune and deconvolute functions should accept marker_exclusivity parameter."""
        import inspect
        from CITEgeist.model.gurobi_impl import (
            deconvolute_local_cell_proportions_per_marker,
            finetune_cell_proportions_per_marker,
        )

        deconv_sig = inspect.signature(deconvolute_local_cell_proportions_per_marker)
        assert "marker_exclusivity" in deconv_sig.parameters, \
            "deconvolute_local_cell_proportions_per_marker missing marker_exclusivity parameter"

        finetune_sig = inspect.signature(finetune_cell_proportions_per_marker)
        assert "marker_exclusivity" in finetune_sig.parameters, \
            "finetune_cell_proportions_per_marker missing marker_exclusivity parameter"

    def test_none_exclusivity_is_backward_compatible(self):
        """Passing marker_exclusivity=None should not change behavior."""
        import inspect
        from CITEgeist.model.gurobi_impl import (
            deconvolute_local_cell_proportions_per_marker,
            finetune_cell_proportions_per_marker,
        )

        deconv_sig = inspect.signature(deconvolute_local_cell_proportions_per_marker)
        param = deconv_sig.parameters["marker_exclusivity"]
        assert param.default is None, f"Default should be None, got {param.default}"

        finetune_sig = inspect.signature(finetune_cell_proportions_per_marker)
        param = finetune_sig.parameters["marker_exclusivity"]
        assert param.default is None, f"Default should be None, got {param.default}"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_marker_exclusivity.py::TestExclusivityInFinetuning -v`
Expected: FAIL with `AssertionError: ... missing marker_exclusivity parameter`

**Step 3: Modify `deconvolute_local_cell_proportions_per_marker()`**

In `gurobi_impl.py`, add `marker_exclusivity` parameter to the signature at line 1231:

```python
def deconvolute_local_cell_proportions_per_marker(
    spot_idx: int,
    adata: sc.AnnData,
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    radius: float = 2.0,
    tolerance: float = 1e-4,
    lambda_reg: float = 1.0,
    alpha: float = 0.5,
    beta_values: Optional[np.ndarray] = None,
    beta_vary: bool = True,
    normalize_beta: bool = True,
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    max_iterations: int = 20,
    max_y_change: float = 0.4,
    marker_exclusivity: Optional[np.ndarray] = None,
) -> Optional[np.ndarray]:
```

Then modify the loss weight at line 1348 from:

```python
                    weight = 1.0 / (n_owners * markers_per_celltype[j])
```

to:

```python
                    excl = marker_exclusivity[m] if marker_exclusivity is not None else 1.0
                    weight = excl / (n_owners * markers_per_celltype[j])
```

**Step 4: Modify `finetune_cell_proportions_per_marker()`**

Add `marker_exclusivity` parameter to the signature at line 1449 (before `max_workers`):

```python
    beta_min: float = 0.1,
    beta_max: float = 2.0,
    marker_exclusivity: Optional[np.ndarray] = None,
    max_workers: Optional[int] = None,
```

Pass it through in `executor.submit` at line 1536, adding after `max_y_change=max_y_change,`:

```python
                        marker_exclusivity=marker_exclusivity,
```

**Step 5: Run test to verify it passes**

Run: `python -m pytest tests/test_marker_exclusivity.py -v`
Expected: All tests PASS (both Task 1 and Task 2 tests)

**Step 6: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py tests/test_marker_exclusivity.py
git commit -m "feat: wire marker_exclusivity parameter into finetuning functions"
```

---

### Task 3: Wire exclusivity into orchestration layer

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:494-538` (compute and pass exclusivity between global pass and finetuning)

**Step 1: Write the failing test**

Add to `tests/test_marker_exclusivity.py`:

```python
class TestOrchestrationIntegration:
    """Test that CitegeistModel computes and stores exclusivity."""

    def test_model_stores_exclusivity_in_results(self):
        """After run_cell_proportion_model, results should contain marker_exclusivity."""
        import inspect
        from CITEgeist.model.citegeist_model import CitegeistModel

        # Verify the orchestration code references compute_marker_exclusivity
        source = inspect.getsource(CitegeistModel.run_cell_proportion_model)
        assert "compute_marker_exclusivity" in source, \
            "run_cell_proportion_model should call compute_marker_exclusivity"
        assert "marker_exclusivity" in source, \
            "run_cell_proportion_model should pass marker_exclusivity to finetuning"
```

**Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_marker_exclusivity.py::TestOrchestrationIntegration -v`
Expected: FAIL with `AssertionError: run_cell_proportion_model should call compute_marker_exclusivity`

**Step 3: Modify `run_cell_proportion_model()` in `citegeist_model.py`**

Add import at top of file (with other gurobi_impl imports):

```python
from .gurobi_impl import compute_marker_exclusivity
```

After line 497 (`self.results["marker_beta"] = marker_beta_dict`), add:

```python
                # Compute marker exclusivity scores for finetuning
                marker_owners = []
                for m_idx in range(assignment_matrix.shape[0]):
                    owners = [j for j in range(assignment_matrix.shape[1]) if assignment_matrix[m_idx, j] > 0]
                    marker_owners.append(owners)

                marker_exclusivity = compute_marker_exclusivity(
                    marker_level_data=marker_level_data,
                    Y_values=Y_values,
                    marker_owners=marker_owners,
                    assignment_matrix=assignment_matrix,
                )

                # Log exclusivity scores
                for m_idx, m_name in enumerate(marker_names):
                    if marker_owners[m_idx]:
                        logging.info(f"  Marker exclusivity: {m_name} = {marker_exclusivity[m_idx]:.3f}")
                self.results["marker_exclusivity"] = {
                    marker_names[i]: marker_exclusivity[i] for i in range(len(marker_names))
                }
```

Then modify the `finetune_cell_proportions_per_marker()` call at line 517-538 to pass exclusivity. Add after `beta_max=beta_max,` (line 533):

```python
                    marker_exclusivity=marker_exclusivity,
```

**Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_marker_exclusivity.py -v`
Expected: All tests PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/citegeist_model.py tests/test_marker_exclusivity.py
git commit -m "feat: compute marker exclusivity after global EM and pass to finetuning"
```

---

### Task 4: Run benchmark and validate

**Files:**
- Modify: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py` (if needed)
- Use existing: evaluation scripts in `Benchmarking/xenium_benchmarking/evaluation/`

**Step 1: Run existing unit tests to verify no regressions**

Run via SLURM:
```bash
sbatch tests/sbatch_module2c_test.sh  # or equivalent test runner
```
Expected: All existing tests PASS

**Step 2: Run Xenium benchmark (5 regions)**

Submit benchmark job using existing SLURM scripts. This runs CITEgeist with the new exclusivity weighting on the Xenium pseudo-Visium data.

**Step 3: Evaluate results**

Compare against:
- Pre-fix baseline (output_achievable_7): Overall r=0.412, Fibroblasts r=0.360
- Bug1+Bug2 without exclusivity (output/manual): Overall r=0.382, Fibroblasts r=0.156

**Success criteria:**
- Fibroblast r improves from 0.156 toward 0.360
- Other 6 cell types maintain their Bug1+Bug2 improvements
- Overall r improves from 0.382

**Step 4: Commit benchmark results**

```bash
git add -A  # benchmark outputs
git commit -m "bench: validate marker exclusivity weighting on Xenium data"
```
