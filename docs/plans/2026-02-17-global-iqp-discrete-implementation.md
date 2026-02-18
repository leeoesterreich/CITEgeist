# Global IQP Discrete Cell Assignment Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace per-spot IQP with global joint IQP to improve mixed dataset proportion correlation from ~0.43 to ≥0.65.

**Architecture:** Add `solve_discrete_cell_counts_global()` function that formulates a single Gurobi IQP over all N×T integer variables. Integrate into existing EM loop by replacing per-spot solve with global solve. API changes propagate through `optimize_discrete_cell_assignment_em()` and `CitegeistModel.run_discrete_cell_assignment()`.

**Tech Stack:** Python 3.10, Gurobi 12.0.3, NumPy, pytest

---

## Task 1: Add Global IQP Solver Function

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py` (add after line 3355)
- Test: `CITEgeist/tests/test_discrete_global_solve.py` (new file)

### Step 1.1: Write the basic test for global solve

```python
# CITEgeist/tests/test_discrete_global_solve.py
"""Tests for global IQP discrete cell assignment."""

import numpy as np
import pytest

# Will import after implementation
# from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts_global


def test_global_solve_returns_valid_counts():
    """Global IQP returns integer counts summing to nuclei per spot."""
    # Small test case: 10 spots, 3 markers, 3 cell types
    N, M, T = 10, 6, 3

    # Marker data: each cell type has 2 markers
    np.random.seed(42)
    marker_level_data = np.random.rand(N, M).astype(np.float64)

    # Assignment matrix: markers 0-1 -> type 0, 2-3 -> type 1, 4-5 -> type 2
    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    assignment_matrix[0:2, 0] = 1  # Type 0
    assignment_matrix[2:4, 1] = 1  # Type 1
    assignment_matrix[4:6, 2] = 1  # Type 2

    marker_names = [f"marker_{i}" for i in range(M)]
    cell_type_names = ["TypeA", "TypeB", "TypeC"]
    nuclei_counts = np.array([5, 3, 7, 4, 6, 2, 8, 5, 3, 4], dtype=np.int64)
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)

    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts_global

    c_values = solve_discrete_cell_counts_global(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        alpha_values=alpha_values,
        time_limit=60.0,
        mip_gap=0.05,
    )

    # Verify output shape and type
    assert c_values.shape == (N, T)
    assert c_values.dtype == np.int64

    # Verify non-negative
    assert (c_values >= 0).all()

    # Verify sum equals nuclei per spot
    row_sums = c_values.sum(axis=1)
    np.testing.assert_array_equal(row_sums, nuclei_counts)
```

### Step 1.2: Run test to verify it fails

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py::test_global_solve_returns_valid_counts -v`

Expected: FAIL with `ImportError: cannot import name 'solve_discrete_cell_counts_global'`

### Step 1.3: Implement the global solve function

Add to `CITEgeist/model/gurobi_impl.py` after line 3355 (after `solve_discrete_cell_counts` function ends):

```python
def solve_discrete_cell_counts_global(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,
    beta_values: np.ndarray,
    alpha_values: Optional[np.ndarray] = None,
    time_limit: float = 300.0,
    mip_gap: float = 0.05,
    prev_c_values: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Solve global IQP for discrete cell counts across all spots jointly.

    Unlike solve_discrete_cell_counts() which optimizes each spot independently,
    this function formulates a single IQP over all N×T integer variables. This
    enforces globally consistent marker-celltype behavior across the tissue.

    Mathematical formulation:
        minimize    Σᵢ Σₘ (X[i,m] - α[m] - β[m] × Σₜ c[i,t] × profile[t,m])²
        subject to  Σₜ c[i,t] = Nᵢ     ∀ spots i
                    c[i,t] ∈ Z≥0       ∀ i, t
                    c[i,t] ≤ Nᵢ        ∀ i, t

    Args:
        marker_level_data: (N, M) antibody data (preprocessed for discrete mode)
        marker_names: List of marker names (length M)
        assignment_matrix: (M, T) binary matrix where A[m,t]=1 if marker m belongs to type t
        cell_type_names: List of cell type names (length T)
        nuclei_counts: (N,) integer nuclei count per spot
        beta_values: (M,) per-marker scaling factors
        alpha_values: (M,) per-marker baselines (optional, defaults to zeros)
        time_limit: Maximum solver time in seconds (default: 300)
        mip_gap: Acceptable optimality gap (default: 0.05 = 5%)
        prev_c_values: (N, T) previous cell counts for warm-start (optional)

    Returns:
        c_values: (N, T) integer cell counts per type per spot
    """
    N, M = marker_level_data.shape
    T = len(cell_type_names)

    # Input validation
    if assignment_matrix.shape != (M, T):
        raise ValueError(f"assignment_matrix shape {assignment_matrix.shape} != expected ({M}, {T})")
    if nuclei_counts.shape != (N,):
        raise ValueError(f"nuclei_counts shape {nuclei_counts.shape} != expected ({N},)")
    if beta_values.shape != (M,):
        raise ValueError(f"beta_values shape {beta_values.shape} != expected ({M},)")

    if alpha_values is None:
        alpha_values = np.zeros(M, dtype=np.float64)

    # Build profile matrix: profile[t, m] = assignment_matrix[m, t]
    profile_matrix = assignment_matrix.T  # Shape: (T, M)

    logging.info(f"Global IQP: {N} spots × {T} cell types = {N * T} integer variables")
    logging.info(f"Time limit: {time_limit}s, MIP gap: {mip_gap:.1%}")

    # Create Gurobi model
    model = gp.Model("global_discrete_cell_assignment")
    model.setParam("OutputFlag", 1)  # Show progress for long solves
    model.setParam("TimeLimit", time_limit)
    model.setParam("MIPGap", mip_gap)
    model.setParam("Threads", 0)  # Use all available cores

    # Create integer variables: c[i, t] = count of cell type t at spot i
    c = {}
    for i in range(N):
        N_i = int(nuclei_counts[i])
        for t in range(T):
            c[i, t] = model.addVar(
                lb=0,
                ub=N_i,
                vtype=GRB.INTEGER,
                name=f"c_{i}_{t}"
            )

    model.update()

    # Warm-start from previous solution or uniform distribution
    for i in range(N):
        N_i = int(nuclei_counts[i])
        if prev_c_values is not None:
            for t in range(T):
                c[i, t].Start = int(prev_c_values[i, t])
        else:
            # Uniform distribution
            base = N_i // T
            remainder = N_i % T
            for t in range(T):
                c[i, t].Start = base + (1 if t < remainder else 0)

    # Constraint: sum of counts equals nuclei per spot
    for i in range(N):
        N_i = int(nuclei_counts[i])
        if N_i > 0:
            model.addConstr(
                quicksum(c[i, t] for t in range(T)) == N_i,
                name=f"nuclei_{i}"
            )

    # Objective: minimize reconstruction error across all spots and markers
    # For each (spot, marker): error = (X[i,m] - α[m] - β[m] × pred[i,m])²
    # where pred[i,m] = Σₜ c[i,t] × profile[t,m]
    error_terms = []
    for i in range(N):
        X_i = marker_level_data[i, :]  # Shape: (M,)
        for m in range(M):
            # Predicted signal from cell counts
            pred_im = quicksum(c[i, t] * profile_matrix[t, m] for t in range(T))
            # Residual: observed - baseline - beta * predicted
            residual = X_i[m] - alpha_values[m] - beta_values[m] * pred_im
            error_terms.append(residual * residual)

    model.setObjective(quicksum(error_terms), GRB.MINIMIZE)

    # Solve
    logging.info("Starting global IQP optimization...")
    model.optimize()

    # Check solution status
    if model.status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL]:
        # Extract solution
        c_values = np.zeros((N, T), dtype=np.int64)
        for i in range(N):
            for t in range(T):
                c_values[i, t] = int(round(c[i, t].X))

        # Verify and fix sum constraint after rounding (should be minimal adjustment)
        for i in range(N):
            N_i = int(nuclei_counts[i])
            current_sum = c_values[i, :].sum()
            if current_sum != N_i:
                diff = N_i - current_sum
                if diff > 0:
                    # Add to largest value
                    max_idx = np.argmax(c_values[i, :])
                    c_values[i, max_idx] += diff
                else:
                    # Remove from largest non-zero
                    sorted_idx = np.argsort(c_values[i, :])[::-1]
                    for idx in sorted_idx:
                        if c_values[i, idx] > 0:
                            remove = min(-diff, c_values[i, idx])
                            c_values[i, idx] -= remove
                            diff += remove
                            if diff == 0:
                                break

        gap = model.MIPGap if hasattr(model, 'MIPGap') and model.SolCount > 0 else float('inf')
        logging.info(f"Global IQP solved: status={model.status}, "
                     f"gap={gap:.2%}, time={model.Runtime:.1f}s")

        return c_values

    else:
        raise RuntimeError(f"Global IQP failed: status={model.status}")
```

### Step 1.4: Run test to verify it passes

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py::test_global_solve_returns_valid_counts -v`

Expected: PASS

### Step 1.5: Commit

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_discrete_global_solve.py
git commit -m "feat(discrete): add solve_discrete_cell_counts_global() function

Implements global IQP solver that optimizes all spots jointly in a single
optimization problem, enforcing globally consistent marker-celltype behavior.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 2: Add Time Limit Test

**Files:**
- Modify: `CITEgeist/tests/test_discrete_global_solve.py`

### Step 2.1: Write time limit test

Add to `CITEgeist/tests/test_discrete_global_solve.py`:

```python
def test_global_solve_respects_time_limit():
    """Global IQP respects time limit and returns feasible solution."""
    # Larger case to potentially hit time limit
    N, M, T = 100, 18, 9

    np.random.seed(123)
    marker_level_data = np.random.rand(N, M).astype(np.float64)

    # 2 markers per cell type
    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    for t in range(T):
        assignment_matrix[t * 2, t] = 1
        assignment_matrix[t * 2 + 1, t] = 1

    marker_names = [f"marker_{i}" for i in range(M)]
    cell_type_names = [f"Type_{t}" for t in range(T)]
    nuclei_counts = np.random.randint(3, 15, size=N).astype(np.int64)
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)

    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts_global

    # Short time limit
    c_values = solve_discrete_cell_counts_global(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        alpha_values=alpha_values,
        time_limit=5.0,  # Very short
        mip_gap=0.10,
    )

    # Should still return valid solution
    assert c_values.shape == (N, T)
    assert (c_values >= 0).all()
    row_sums = c_values.sum(axis=1)
    np.testing.assert_array_equal(row_sums, nuclei_counts)
```

### Step 2.2: Run test

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py::test_global_solve_respects_time_limit -v`

Expected: PASS

### Step 2.3: Commit

```bash
git add CITEgeist/tests/test_discrete_global_solve.py
git commit -m "test(discrete): add time limit test for global IQP

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 3: Add Warm-Start Test

**Files:**
- Modify: `CITEgeist/tests/test_discrete_global_solve.py`

### Step 3.1: Write warm-start test

Add to `CITEgeist/tests/test_discrete_global_solve.py`:

```python
def test_global_solve_with_warm_start():
    """Global IQP accepts warm-start from previous iteration."""
    N, M, T = 20, 6, 3

    np.random.seed(456)
    marker_level_data = np.random.rand(N, M).astype(np.float64)

    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    assignment_matrix[0:2, 0] = 1
    assignment_matrix[2:4, 1] = 1
    assignment_matrix[4:6, 2] = 1

    marker_names = [f"marker_{i}" for i in range(M)]
    cell_type_names = ["TypeA", "TypeB", "TypeC"]
    nuclei_counts = np.random.randint(3, 10, size=N).astype(np.int64)
    beta_values = np.ones(M, dtype=np.float64)
    alpha_values = np.zeros(M, dtype=np.float64)

    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts_global

    # First solve without warm-start
    c_values_1 = solve_discrete_cell_counts_global(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        alpha_values=alpha_values,
        time_limit=30.0,
        mip_gap=0.05,
        prev_c_values=None,
    )

    # Second solve with warm-start from first
    c_values_2 = solve_discrete_cell_counts_global(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        alpha_values=alpha_values,
        time_limit=30.0,
        mip_gap=0.05,
        prev_c_values=c_values_1,
    )

    # Both should be valid
    assert c_values_2.shape == (N, T)
    np.testing.assert_array_equal(c_values_2.sum(axis=1), nuclei_counts)

    # With same inputs and warm-start, should get same or similar solution
    # (high correlation)
    corr = np.corrcoef(c_values_1.flatten(), c_values_2.flatten())[0, 1]
    assert corr > 0.95, f"Expected high correlation, got {corr}"
```

### Step 3.2: Run test

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py::test_global_solve_with_warm_start -v`

Expected: PASS

### Step 3.3: Commit

```bash
git add CITEgeist/tests/test_discrete_global_solve.py
git commit -m "test(discrete): add warm-start test for global IQP

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 4: Integrate Global Solve into EM

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py` (lines 3358-3540, `optimize_discrete_cell_assignment_em`)

### Step 4.1: Write integration test

Add to `CITEgeist/tests/test_discrete_global_solve.py`:

```python
def test_em_with_global_solve():
    """EM algorithm works with global_solve=True."""
    N, M, T = 50, 6, 3

    np.random.seed(789)
    marker_level_data = np.random.rand(N, M).astype(np.float64)

    assignment_matrix = np.zeros((M, T), dtype=np.float64)
    assignment_matrix[0:2, 0] = 1
    assignment_matrix[2:4, 1] = 1
    assignment_matrix[4:6, 2] = 1

    marker_names = [f"marker_{i}" for i in range(M)]
    cell_type_names = ["TypeA", "TypeB", "TypeC"]
    nuclei_counts = np.random.randint(3, 12, size=N).astype(np.int64)

    from CITEgeist.model.gurobi_impl import optimize_discrete_cell_assignment_em

    c_values, beta_values, marker_beta_dict, alpha_values = optimize_discrete_cell_assignment_em(
        marker_level_data=marker_level_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        max_em_iterations=3,
        global_solve=True,
        global_time_limit=30.0,
        global_mip_gap=0.10,
    )

    # Verify outputs
    assert c_values.shape == (N, T)
    assert c_values.dtype == np.int64
    assert (c_values >= 0).all()
    np.testing.assert_array_equal(c_values.sum(axis=1), nuclei_counts)

    assert beta_values.shape == (M,)
    assert len(marker_beta_dict) == M
    assert alpha_values.shape == (M,)
```

### Step 4.2: Run test to verify it fails

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py::test_em_with_global_solve -v`

Expected: FAIL with `TypeError: optimize_discrete_cell_assignment_em() got an unexpected keyword argument 'global_solve'`

### Step 4.3: Modify EM function signature and E-step

Edit `CITEgeist/model/gurobi_impl.py` at line 3358. Update function signature:

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
    timeout_per_spot: float = 60.0,
    lambda_sparse: float = 0.0,
    prior_proportions: Optional[np.ndarray] = None,
    lambda_prior: float = 0.0,
    global_solve: bool = True,
    global_time_limit: float = 300.0,
    global_mip_gap: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray]:
```

Update docstring to include new parameters:

```python
        global_solve: If True (default), use global IQP solver. If False, use per-spot IQP.
        global_time_limit: Time limit for global solver in seconds (default: 300).
        global_mip_gap: Acceptable MIP gap for global solver (default: 0.05 = 5%).
```

Replace E-step (around lines 3450-3465) with conditional:

```python
        # ==================== E-Step ====================
        if global_solve:
            # Global solve: single IQP over all spots
            c_values = solve_discrete_cell_counts_global(
                marker_level_data=marker_level_data,
                marker_names=marker_names,
                assignment_matrix=assignment_matrix,
                cell_type_names=cell_type_names,
                nuclei_counts=nuclei_counts,
                beta_values=beta_values,
                alpha_values=alpha_values,
                time_limit=global_time_limit,
                mip_gap=global_mip_gap,
                prev_c_values=c_values if iteration > 0 else None,
            )
        else:
            # Per-spot solve (original behavior)
            c_values = solve_discrete_cell_counts(
                marker_level_data=marker_level_data,
                marker_names=marker_names,
                assignment_matrix=assignment_matrix,
                cell_type_names=cell_type_names,
                nuclei_counts=nuclei_counts,
                beta_values=beta_values,
                alpha_values=alpha_values,
                max_nuclei_cap=max_nuclei_cap,
                timeout_per_spot=timeout_per_spot,
                lambda_sparse=lambda_sparse,
                prior_proportions=prior_proportions,
                lambda_prior=lambda_prior,
            )
```

### Step 4.4: Run test to verify it passes

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py::test_em_with_global_solve -v`

Expected: PASS

### Step 4.5: Run all tests

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_global_solve.py -v`

Expected: All 4 tests PASS

### Step 4.6: Commit

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_discrete_global_solve.py
git commit -m "feat(discrete): integrate global IQP into EM algorithm

Add global_solve, global_time_limit, global_mip_gap parameters to
optimize_discrete_cell_assignment_em(). Default is global_solve=True
for joint optimization across all spots.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 5: Update CitegeistModel API

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py` (lines 942-1056, `run_discrete_cell_assignment`)

### Step 5.1: Update function signature

Edit `CITEgeist/model/citegeist_model.py` at line 942. Update signature:

```python
    def run_discrete_cell_assignment(
        self,
        nuclei_counts: Optional[pd.Series] = None,
        max_em_iterations: int = 20,
        beta_convergence_tol: float = 1e-3,
        max_nuclei_cap: int = 30,
        beta_min: float = 0.1,
        beta_max: float = 2.0,
        timeout_per_spot: float = 60.0,
        lambda_sparse: float = 0.0,
        prior_proportions: Optional[np.ndarray] = None,
        lambda_prior: float = 0.0,
        global_solve: bool = True,
        global_time_limit: float = 300.0,
        global_mip_gap: float = 0.05,
    ) -> pd.DataFrame:
```

Update docstring to add:

```python
            global_solve: If True (default), use global IQP solver for joint
                optimization across all spots. If False, use per-spot IQP
                (original behavior). Global solve enforces globally consistent
                marker-celltype relationships.
            global_time_limit: Time limit for global solver in seconds (default: 300).
            global_mip_gap: Acceptable MIP gap for global solver (default: 0.05).
```

### Step 5.2: Pass new parameters to EM call

Update the call to `optimize_discrete_cell_assignment_em()` around line 1041:

```python
        # Run EM optimization
        c_values, beta_values, marker_beta_dict, alpha_values = optimize_discrete_cell_assignment_em(
            marker_level_data=marker_level_data,
            marker_names=marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_array,
            max_em_iterations=max_em_iterations,
            beta_convergence_tol=beta_convergence_tol,
            beta_min=beta_min,
            beta_max=beta_max,
            max_nuclei_cap=max_nuclei_cap,
            timeout_per_spot=timeout_per_spot,
            lambda_sparse=lambda_sparse,
            prior_proportions=prior_proportions,
            lambda_prior=lambda_prior,
            global_solve=global_solve,
            global_time_limit=global_time_limit,
            global_mip_gap=global_mip_gap,
        )
```

### Step 5.3: Commit

```bash
git add CITEgeist/model/citegeist_model.py
git commit -m "feat(discrete): expose global_solve params in CitegeistModel API

Add global_solve, global_time_limit, global_mip_gap to
run_discrete_cell_assignment() method.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 6: Add CLI Flag to Benchmark Script

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py`

### Step 6.1: Add argument parser flag

Find the argparse section (around line 517-534) and add:

```python
    parser.add_argument("--global-solve", action="store_true", default=True,
                        help="Use global IQP solver (default: True)")
    parser.add_argument("--no-global-solve", action="store_false", dest="global_solve",
                        help="Use per-spot IQP solver instead of global")
    parser.add_argument("--global-time-limit", type=float, default=300.0,
                        help="Time limit for global IQP solver in seconds (default: 300)")
    parser.add_argument("--global-mip-gap", type=float, default=0.05,
                        help="MIP gap tolerance for global solver (default: 0.05)")
```

### Step 6.2: Pass to model call

Find where `model.run_discrete_cell_assignment()` is called and add the new parameters:

```python
        cell_counts_df = model.run_discrete_cell_assignment(
            nuclei_counts=nuclei_per_spot,
            max_em_iterations=args.max_em_iterations,
            lambda_sparse=args.lambda_sparse,
            lambda_prior=args.lambda_prior,
            prior_proportions=prior_props,
            global_solve=args.global_solve,
            global_time_limit=args.global_time_limit,
            global_mip_gap=args.global_mip_gap,
        )
```

### Step 6.3: Commit

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): add --global-solve CLI flags to discrete benchmark

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 7: Create SLURM Benchmark Scripts

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_global_mixed.sh`
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_global_high_seg.sh`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_global_xenium.sh`

### Step 7.1: Create mixed benchmark script

```bash
#!/bin/bash
#SBATCH --job-name=discrete_global_mixed
#SBATCH --output=logs/discrete_global_mixed_%A_%a.out
#SBATCH --error=logs/discrete_global_mixed_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id ${SLURM_ARRAY_TASK_ID} \
    --condition mixed \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global/mixed/dapi \
    --global-solve \
    --global-time-limit 600 \
    --global-mip-gap 0.05 \
    --max-em-iterations 10
```

### Step 7.2: Create high_seg benchmark script

```bash
#!/bin/bash
#SBATCH --job-name=discrete_global_highseg
#SBATCH --output=logs/discrete_global_highseg_%A_%a.out
#SBATCH --error=logs/discrete_global_highseg_%A_%a.err
#SBATCH --array=0-4
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id ${SLURM_ARRAY_TASK_ID} \
    --condition high_seg \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global/high_seg/dapi \
    --global-solve \
    --global-time-limit 600 \
    --global-mip-gap 0.05 \
    --max-em-iterations 10
```

### Step 7.3: Create xenium benchmark script

```bash
#!/bin/bash
#SBATCH --job-name=discrete_global_xenium
#SBATCH --output=logs/discrete_global_xenium_%A_%a.out
#SBATCH --error=logs/discrete_global_xenium_%A_%a.err
#SBATCH --array=0-13
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --partition=htc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3

eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/xenium_benchmarking/CITEgeist/src/benchmark_discrete_xenium.py \
    --region-id ${SLURM_ARRAY_TASK_ID} \
    --output-dir Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_global \
    --global-solve \
    --global-time-limit 600 \
    --global-mip-gap 0.05 \
    --max-em-iterations 10
```

### Step 7.4: Commit

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_global_mixed.sh
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_discrete_global_high_seg.sh
git add Benchmarking/xenium_benchmarking/CITEgeist/slurm/sbatch_discrete_global_xenium.sh
git commit -m "chore(benchmark): add SLURM scripts for global IQP benchmarks

Scripts for mixed, high_seg, and xenium datasets with global_solve=True.

Co-Authored-By: Claude Opus 4.5 <noreply@anthropic.com>"
```

---

## Task 8: Run Benchmarks and Validate

### Step 8.1: Create output directories

```bash
mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global/mixed/dapi
mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/output_discrete_global/high_seg/dapi
mkdir -p Benchmarking/simulation_benchmarking/CITEgeist/slurm/logs
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/output_discrete_global
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/slurm/logs
```

### Step 8.2: Submit jobs

```bash
cd Benchmarking/simulation_benchmarking/CITEgeist/slurm
sbatch sbatch_discrete_global_mixed.sh
sbatch sbatch_discrete_global_high_seg.sh

cd ../../xenium_benchmarking/CITEgeist/slurm
sbatch sbatch_discrete_global_xenium.sh
```

### Step 8.3: Check results after completion

```bash
# Check mixed results
for f in output_discrete_global/mixed/dapi/Wu_rep_*/benchmark_results.json; do
    python3 -c "import json; d=json.load(open('$f')); print(f'{f}: prop_corr={d[\"proportion_correlation\"]:.3f}')"
done

# Compare to baseline
echo "Baseline (per-spot):"
for f in output_discrete_geomfix_nolegacy/mixed/dapi/Wu_rep_*/benchmark_results.json; do
    python3 -c "import json; d=json.load(open('$f')); print(f'{f}: prop_corr={d[\"proportion_correlation\"]:.3f}')"
done
```

Expected results:
- mixed: prop_corr ≥0.65 (up from ~0.43)
- high_seg: prop_corr ≥0.76 (no regression)
- xenium: prop_corr ≥0.66 (no regression)

---

## Summary

| Task | Description | Files Changed |
|------|-------------|---------------|
| 1 | Add `solve_discrete_cell_counts_global()` | gurobi_impl.py, test_discrete_global_solve.py |
| 2 | Add time limit test | test_discrete_global_solve.py |
| 3 | Add warm-start test | test_discrete_global_solve.py |
| 4 | Integrate into EM | gurobi_impl.py, test_discrete_global_solve.py |
| 5 | Update CitegeistModel API | citegeist_model.py |
| 6 | Add CLI flags | benchmark_discrete_simulation.py |
| 7 | Create SLURM scripts | 3 new sbatch scripts |
| 8 | Run and validate | N/A |

**Total commits:** 7
**Total new lines:** ~250
