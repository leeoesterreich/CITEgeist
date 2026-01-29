# Adaptive Per-Marker Baseline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add per-marker intercept (alpha) to the reconstruction model so ubiquitous markers like VIM learn a baseline, preventing Fibroblast over-prediction.

**Architecture:** Extend the EM loop to jointly learn alpha (baseline) and beta (slope) per marker. Alpha is learned in global EM and passed as fixed input to finetuning. The QP structure is unchanged — alpha just baseline-subtracts the signal.

**Tech Stack:** numpy, gurobipy, pytest

---

### Task 1: Add per-marker alpha to global EM

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:527-766` (`optimize_cell_proportions_per_marker`)
- Test: `tests/test_marker_baseline.py` (create)

**Step 1: Write failing tests**

Create `tests/test_marker_baseline.py`:

```python
"""Tests for adaptive per-marker baseline (alpha/intercept)."""
import numpy as np
import pytest

from CITEgeist.model.gurobi_impl import optimize_cell_proportions_per_marker


class TestMarkerBaseline:
    """Test that per-marker alpha (baseline) is learned correctly."""

    def _make_data(self, N=200, seed=42):
        """Create synthetic data with one ubiquitous marker and one specific marker.

        Markers:
        - Marker 0: specific to cell type 0 (signal ~ Y[:,0], baseline ~0)
        - Marker 1: ubiquitous with variation (baseline=0.6, signal += 0.3*Y[:,1])
        - Marker 2: specific to cell type 2 (signal ~ Y[:,2], baseline ~0)
        """
        rng = np.random.RandomState(seed)
        T = 3

        # Ground truth proportions
        Y_true = rng.dirichlet([1, 1, 1], size=N)

        # Marker signals
        M = 3
        S = np.zeros((N, M))
        S[:, 0] = Y_true[:, 0] * 0.8 + rng.normal(0, 0.02, N)  # specific
        S[:, 1] = 0.6 + Y_true[:, 1] * 0.3 + rng.normal(0, 0.02, N)  # ubiquitous with baseline
        S[:, 2] = Y_true[:, 2] * 0.8 + rng.normal(0, 0.02, N)  # specific
        S = np.clip(S, 0, 1)

        # Normalize per column (like map_antibodies_to_profiles_v2)
        col_max = np.max(S, axis=0)
        S = S / col_max

        assignment = np.eye(M, T)
        cell_type_names = ["TypeA", "TypeB", "TypeC"]
        marker_names = ["specific_A", "ubiquitous_B", "specific_C"]

        return S, assignment, cell_type_names, marker_names, Y_true

    def test_returns_alpha_values(self):
        """optimize_cell_proportions_per_marker should return alpha array."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        result = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=5, warn_only=True,
        )
        # Should now return (Y_values, beta_values, marker_beta_dict, alpha_values)
        assert len(result) == 4, f"Expected 4 return values, got {len(result)}"
        alpha_values = result[3]
        assert alpha_values.shape == (3,), f"Alpha shape should be (3,), got {alpha_values.shape}"

    def test_ubiquitous_marker_learns_nonzero_alpha(self):
        """Ubiquitous marker should learn alpha > 0."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        _, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=10, warn_only=True,
        )
        assert alpha_values[1] > 0.1, (
            f"Ubiquitous marker alpha={alpha_values[1]:.3f}, expected > 0.1"
        )

    def test_specific_markers_learn_near_zero_alpha(self):
        """Specific markers should learn alpha near 0."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        _, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=10, warn_only=True,
        )
        assert alpha_values[0] < 0.15, (
            f"Specific marker 0 alpha={alpha_values[0]:.3f}, expected < 0.15"
        )
        assert alpha_values[2] < 0.15, (
            f"Specific marker 2 alpha={alpha_values[2]:.3f}, expected < 0.15"
        )

    def test_alpha_clipped_to_range(self):
        """Alpha should be in [0, alpha_max]."""
        S, assignment, ct_names, m_names, _ = self._make_data()
        _, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=10, warn_only=True,
            alpha_max=0.8,
        )
        assert np.all(alpha_values >= 0), f"Alpha min={alpha_values.min():.3f} < 0"
        assert np.all(alpha_values <= 0.8), f"Alpha max={alpha_values.max():.3f} > 0.8"

    def test_proportions_improve_with_baseline(self):
        """With baseline, proportions for ubiquitous marker's cell type should be more accurate."""
        S, assignment, ct_names, m_names, Y_true = self._make_data(N=300)
        Y_values, _, _, alpha_values = optimize_cell_proportions_per_marker(
            marker_level_data=S, marker_names=m_names,
            assignment_matrix=assignment, cell_type_names=ct_names,
            max_iterations=15, warn_only=True,
        )
        # TypeB (ubiquitous marker owner) should correlate with ground truth
        r = np.corrcoef(Y_true[:, 1], Y_values[:, 1])[0, 1]
        assert r > 0.3, f"TypeB correlation with GT = {r:.3f}, expected > 0.3"
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_marker_baseline.py -v`
Expected: FAIL (wrong number of return values)

**Step 3: Implement alpha in global EM**

In `optimize_cell_proportions_per_marker()`:

1. Add parameters to signature (after `beta_max`):
```python
    alpha_max: float = 0.8,
    lambda_alpha: float = 1.0,
```

2. Initialize alpha after beta initialization (after line 628):
```python
    alpha_values = np.zeros(M, dtype=np.float64)  # per-marker baseline
```

3. In E-step QP (lines 656-659), baseline-subtract the signal:
```python
                    S_im = marker_level_data[i, m] - alpha_values[m]  # baseline-subtracted
```

4. In M-step (lines 708-722), replace simple beta OLS with joint alpha+beta:
```python
        beta_new = np.zeros(M, dtype=np.float64)
        alpha_new = np.zeros(M, dtype=np.float64)
        for m in range(M):
            if not marker_has_owner[m]:
                beta_new[m] = 1.0
                alpha_new[m] = 0.0
                continue

            owners_m = marker_owners[m]
            Y_combined = np.zeros(N, dtype=np.float64)
            for j in owners_m:
                Y_combined += Y_values[:, j]

            S_m = marker_level_data[:, m]

            # OLS: S_m = alpha_m + beta_m * Y_combined
            Y_mean = np.mean(Y_combined)
            S_mean = np.mean(S_m)
            Y_var = np.dot(Y_combined - Y_mean, Y_combined - Y_mean)

            if Y_var > 1e-9:
                beta_new[m] = np.dot(S_m - S_mean, Y_combined - Y_mean) / Y_var
            else:
                beta_new[m] = beta_values[m]  # keep previous

            beta_new[m] = np.clip(beta_new[m], beta_min, beta_max)

            # Alpha with L2 regularization toward zero
            # Regularized: alpha = (S_mean - beta * Y_mean) / (1 + lambda_alpha / N)
            raw_alpha = S_mean - beta_new[m] * Y_mean
            alpha_new[m] = raw_alpha / (1.0 + lambda_alpha / N)
            alpha_new[m] = np.clip(alpha_new[m], 0.0, alpha_max)

        alpha_values = alpha_new.copy()
```

5. Update convergence check to include alpha:
```python
        alpha_diff = np.linalg.norm(alpha_new - alpha_values) if iteration > 0 else float('inf')
```
(Track alpha_prev alongside beta_prev)

6. Update return value (line 766):
```python
    return Y_values, beta_new, marker_beta_dict, alpha_values
```

7. Add alpha logging before return:
```python
    logging.info(f"Alpha (baseline) range: [{alpha_values.min():.3f}, {alpha_values.max():.3f}], mean: {alpha_values.mean():.3f}")
    for m in range(M):
        if marker_has_owner[m] and alpha_values[m] > 0.05:
            logging.info(f"  Marker '{marker_names[m]}': alpha={alpha_values[m]:.3f}, beta={beta_new[m]:.3f}")
```

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_marker_baseline.py -v`
Expected: PASS (all 5 tests)

**Step 5: Commit**

```bash
git add tests/test_marker_baseline.py CITEgeist/model/gurobi_impl.py
git commit -m "feat: add per-marker baseline (alpha) to global EM"
```

---

### Task 2: Wire alpha into finetuning

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:1292-1311` (`deconvolute_local_cell_proportions_per_marker` signature)
- Modify: `CITEgeist/model/gurobi_impl.py:1430-1434` (E-step in finetuning)
- Modify: `CITEgeist/model/gurobi_impl.py:1453-1467` (M-step in finetuning)
- Modify: `CITEgeist/model/gurobi_impl.py:1512-1534` (`finetune_cell_proportions_per_marker` signature)
- Modify: `CITEgeist/model/gurobi_impl.py:1603-1622` (executor.submit call)
- Test: `tests/test_marker_baseline.py` (add tests)

**Step 1: Write failing tests**

Add to `tests/test_marker_baseline.py`:

```python
import inspect
from CITEgeist.model.gurobi_impl import (
    deconvolute_local_cell_proportions_per_marker,
    finetune_cell_proportions_per_marker,
)


class TestBaselineInFinetuning:
    """Test that finetuning functions accept marker_alpha."""

    def test_deconvolute_accepts_marker_alpha(self):
        sig = inspect.signature(deconvolute_local_cell_proportions_per_marker)
        assert "marker_alpha" in sig.parameters
        assert sig.parameters["marker_alpha"].default is None

    def test_finetune_accepts_marker_alpha(self):
        sig = inspect.signature(finetune_cell_proportions_per_marker)
        assert "marker_alpha" in sig.parameters
        assert sig.parameters["marker_alpha"].default is None
```

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_marker_baseline.py::TestBaselineInFinetuning -v`
Expected: FAIL

**Step 3: Implement**

1. Add `marker_alpha: Optional[np.ndarray] = None` parameter to `deconvolute_local_cell_proportions_per_marker()` (after `marker_exclusivity`).

2. In the E-step QP (line 1432-1434), baseline-subtract:
```python
                        S_im = local_marker_data[i, m]
                        if marker_alpha is not None:
                            S_im = S_im - marker_alpha[m]
```

3. In the M-step (line 1464-1466), account for alpha when updating beta:
```python
                    S_m = local_marker_data[:, m]
                    if marker_alpha is not None:
                        S_m = S_m - marker_alpha[m]  # baseline-subtract for beta update
                    denominator = np.dot(Y_combined, Y_combined) + 1e-9
                    new_beta[m] = np.dot(S_m, Y_combined) / denominator
```

4. Add `marker_alpha: Optional[np.ndarray] = None` to `finetune_cell_proportions_per_marker()` (after `marker_exclusivity`).

5. Pass through in executor.submit (after `marker_exclusivity=marker_exclusivity`):
```python
                        marker_alpha=marker_alpha,
```

**Step 4: Run tests**

Run: `pytest tests/test_marker_baseline.py -v`
Expected: PASS (all 7 tests)

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py tests/test_marker_baseline.py
git commit -m "feat: wire marker_alpha into finetuning functions"
```

---

### Task 3: Wire alpha through orchestration

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:521-606`
- Test: `tests/test_marker_baseline.py` (add test)

**Step 1: Write failing test**

```python
class TestOrchestration:
    def test_model_passes_alpha_to_finetuning(self):
        source = inspect.getsource(
            __import__("CITEgeist.model.citegeist_model", fromlist=["CitegeistModel"]).CitegeistModel.run_cell_proportion_model
        )
        assert "alpha_values" in source, "run_cell_proportion_model should handle alpha_values"
        assert "marker_alpha" in source, "run_cell_proportion_model should pass marker_alpha to finetuning"
```

**Step 2: Implement**

1. Update the call to `optimize_cell_proportions_per_marker` (line 521) to unpack 4 values:
```python
                Y_values, beta_values, marker_beta_dict, alpha_values = optimize_cell_proportions_per_marker(
```

2. After storing marker_beta (line 543), log and store alpha:
```python
                # Store marker baselines
                marker_alpha_dict = {marker_names[i]: alpha_values[i] for i in range(len(marker_names))}
                self.results["marker_alpha"] = marker_alpha_dict
                for m_idx, m_name in enumerate(marker_names):
                    if alpha_values[m_idx] > 0.05:
                        logging.info(f"  Marker baseline: {m_name} = {alpha_values[m_idx]:.3f}")
```

3. Pass alpha to finetuning (after `marker_exclusivity=marker_exclusivity`, line 601):
```python
                    marker_alpha=alpha_values,
```

**Step 3: Run tests**

Run: `pytest tests/test_marker_baseline.py -v`
Expected: PASS (all 8 tests)

**Step 4: Commit**

```bash
git add CITEgeist/model/citegeist_model.py tests/test_marker_baseline.py
git commit -m "feat: wire marker alpha through orchestration layer"
```

---

### Task 4: Benchmark and validate

**Step 1: Submit benchmark SLURM job**

Create and submit SLURM script running 5 regions with the updated code, then evaluate against:
1. Current code (WITH alpha baseline)
2. Bug1+Bug2 baseline (output/manual)
3. PRE-FIX baseline (output_achievable_7)

**Step 2: Analyze results**

Compare per-cell-type Pearson r, RMSE, MAE, and bias. Success criteria:
- Fibroblast r improves significantly from 0.156 (Bug1+Bug2)
- Fibroblast bias decreases from +0.205
- Other 6 cell types maintain or improve their metrics
- Overall mean r at least matches pre-fix baseline (0.412)
