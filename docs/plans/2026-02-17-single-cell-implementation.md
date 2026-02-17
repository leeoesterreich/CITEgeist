# Single-Cell Output Mode Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Transform CITEgeist from spot-level layers to true single-cell AnnData output, with each nucleus as one observation.

**Architecture:** Four-phase pipeline: (1) IQP cell counts with sparsity, (2) nucleus-scaled GEX optimization, (3) GEX-informed count refinement, (4) single-cell assembly. Prerequisite Phase 0 fixes the IQP performance regression.

**Tech Stack:** Python 3.10, Gurobi (gurobipy), NumPy, Pandas, AnnData, scanpy

---

## Phase 0: Fix Discrete IQP Sparsity

### Task 0.1: Add Sparsity Parameter to solve_discrete_cell_counts

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:3102-3230`
- Test: `CITEgeist/tests/test_discrete_sparsity.py` (create)

**Step 1: Write the failing test**

Create `CITEgeist/tests/test_discrete_sparsity.py`:

```python
"""Tests for discrete cell assignment with sparsity regularization."""

import numpy as np
import pytest


def test_sparsity_parameter_exists():
    """Verify lambda_sparse parameter is accepted."""
    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts

    # Minimal test data: 2 spots, 3 markers, 2 cell types
    marker_data = np.array([[1.0, 0.5, 0.0], [0.0, 0.5, 1.0]])
    marker_names = ["M1", "M2", "M3"]
    # M1 -> Type A, M2 -> both, M3 -> Type B
    assignment_matrix = np.array([[1, 0], [1, 1], [0, 1]])
    cell_type_names = ["TypeA", "TypeB"]
    nuclei_counts = np.array([3, 3])
    beta_values = np.array([1.0, 1.0, 1.0])

    # Should accept lambda_sparse parameter without error
    result = solve_discrete_cell_counts(
        marker_level_data=marker_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        lambda_sparse=0.1,  # New parameter
    )

    assert result.shape == (2, 2)
    # Each row should sum to nuclei count
    assert result[0].sum() == 3
    assert result[1].sum() == 3


def test_sparsity_encourages_fewer_types():
    """Higher lambda_sparse should result in fewer non-zero cell types per spot."""
    from CITEgeist.model.gurobi_impl import solve_discrete_cell_counts

    # Create ambiguous data where multiple cell types could explain signal
    np.random.seed(42)
    marker_data = np.array([[0.5, 0.5, 0.5, 0.5]])  # Uniform signal
    marker_names = ["M1", "M2", "M3", "M4"]
    # Each marker maps to one of 4 cell types
    assignment_matrix = np.eye(4)
    cell_type_names = ["A", "B", "C", "D"]
    nuclei_counts = np.array([4])
    beta_values = np.ones(4)

    # Without sparsity: might spread cells across types
    result_no_sparse = solve_discrete_cell_counts(
        marker_level_data=marker_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        lambda_sparse=0.0,
    )

    # With high sparsity: should concentrate in fewer types
    result_high_sparse = solve_discrete_cell_counts(
        marker_level_data=marker_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        beta_values=beta_values,
        lambda_sparse=1.0,
    )

    n_types_no_sparse = (result_no_sparse[0] > 0).sum()
    n_types_high_sparse = (result_high_sparse[0] > 0).sum()

    # High sparsity should use fewer or equal cell types
    assert n_types_high_sparse <= n_types_no_sparse
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_sparsity.py -v`

Expected: FAIL with "unexpected keyword argument 'lambda_sparse'"

**Step 3: Add lambda_sparse parameter to function signature**

Modify `CITEgeist/model/gurobi_impl.py` at line 3102:

```python
def solve_discrete_cell_counts(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    assignment_matrix: np.ndarray,
    cell_type_names: List[str],
    nuclei_counts: np.ndarray,
    beta_values: np.ndarray,
    alpha_values: Optional[np.ndarray] = None,
    max_nuclei_cap: int = 30,
    timeout_per_spot: float = 60.0,
    lambda_sparse: float = 0.0,  # NEW: sparsity penalty weight
) -> np.ndarray:
    """
    Solve IQP for discrete cell counts given fixed beta (E-step).

    Mathematical formulation:
        minimize    Σᵢ Σₘ (X[i,m] - α[m] - Σₜ c[i,t] × profile[t,m] × β[m])²
                    + λ_sparse × Σₜ 1{c[i,t] > 0}
        subject to  Σₜ c[i,t] = N_i     ∀ spots i with N_i > 0
                    c[i,t] ∈ Z≥0        ∀ i, t

    Args:
        marker_level_data: (N, M) antibody data (preprocessed for discrete mode)
        marker_names: List of marker names (length M)
        assignment_matrix: (M, T) binary matrix where A[m,t]=1 if marker m belongs to type t
        cell_type_names: List of cell type names (length T)
        nuclei_counts: (N,) integer nuclei count per spot from Cellpose segmentation
        beta_values: (M,) per-marker scaling factors from previous EM iteration
        alpha_values: (M,) per-marker baselines (optional, defaults to zeros)
        max_nuclei_cap: Above this nuclei count, use continuous relaxation + rounding
        timeout_per_spot: Maximum seconds per spot optimization (default: 60)
        lambda_sparse: Weight for sparsity penalty (default: 0.0 = no sparsity).
            Higher values encourage fewer non-zero cell types per spot.

    Returns:
        c_values: (N, T) integer cell counts per type per spot
    """
```

**Step 4: Add indicator variables and sparsity penalty to Gurobi model**

Find the Gurobi model setup section (around line 3168) and modify:

```python
        try:
            model = gp.Model(f"discrete_cell_spot_{i}")
            model.setParam("OutputFlag", 0)
            model.setParam("TimeLimit", timeout_per_spot)

            # Decision variables: c[t] = count of cell type t at spot i
            if use_integer:
                c = model.addVars(T, lb=0, ub=N_i, vtype=GRB.INTEGER, name="c")
            else:
                c = model.addVars(T, lb=0, ub=N_i, vtype=GRB.CONTINUOUS, name="c")

            # Indicator variables for sparsity: y[t] = 1 if c[t] > 0
            if lambda_sparse > 0:
                y = model.addVars(T, vtype=GRB.BINARY, name="y")
                # Link indicators: c[t] > 0 => y[t] = 1
                # Equivalently: c[t] <= N_i * y[t]
                for t in range(T):
                    model.addConstr(c[t] <= N_i * y[t], name=f"indicator_{t}")

            # Constraint: sum of counts equals nuclei count
            model.addConstr(quicksum(c[t] for t in range(T)) == N_i, name="nuclei_sum")

            # Objective: minimize reconstruction error + sparsity penalty
            # For each marker m: error = X[i,m] - Σₜ c[t] × profile[t,m] × β[m]
            obj_terms = []
            for m in range(M):
                pred = quicksum(c[t] * profile_matrix[t, m] * beta_values[m] for t in range(T))
                residual = X_i[m] - pred
                obj_terms.append(residual * residual)

            # Add sparsity penalty: lambda_sparse * sum(y[t])
            if lambda_sparse > 0:
                sparsity_penalty = lambda_sparse * quicksum(y[t] for t in range(T))
                model.setObjective(quicksum(obj_terms) + sparsity_penalty, GRB.MINIMIZE)
            else:
                model.setObjective(quicksum(obj_terms), GRB.MINIMIZE)

            model.optimize()
```

**Step 5: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_sparsity.py -v`

Expected: PASS

**Step 6: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_discrete_sparsity.py
git commit -m "feat(discrete): add sparsity penalty to IQP solver

Add lambda_sparse parameter to solve_discrete_cell_counts() that
penalizes the number of non-zero cell types per spot. This encourages
sparser solutions matching biological reality (2-4 cell types per spot).

Uses Gurobi indicator constraints: c[t] <= N_i * y[t] where y[t] is binary."
```

---

### Task 0.2: Propagate lambda_sparse Through EM Algorithm

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:3237-3400` (optimize_discrete_cell_assignment_em)
- Test: `CITEgeist/tests/test_discrete_sparsity.py` (add test)

**Step 1: Write the failing test**

Add to `CITEgeist/tests/test_discrete_sparsity.py`:

```python
def test_em_accepts_lambda_sparse():
    """Verify EM algorithm propagates lambda_sparse to IQP solver."""
    from CITEgeist.model.gurobi_impl import optimize_discrete_cell_assignment_em

    # Minimal test data
    marker_data = np.array([[1.0, 0.5], [0.5, 1.0]])
    marker_names = ["M1", "M2"]
    assignment_matrix = np.array([[1, 0], [0, 1]])
    cell_type_names = ["TypeA", "TypeB"]
    nuclei_counts = np.array([2, 2])

    # Should accept lambda_sparse parameter
    c_values, beta_values, stats, alpha_values = optimize_discrete_cell_assignment_em(
        marker_level_data=marker_data,
        marker_names=marker_names,
        assignment_matrix=assignment_matrix,
        cell_type_names=cell_type_names,
        nuclei_counts=nuclei_counts,
        lambda_sparse=0.1,  # New parameter
        max_em_iterations=3,
    )

    assert c_values.shape == (2, 2)
    assert c_values[0].sum() == 2
    assert c_values[1].sum() == 2
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_sparsity.py::test_em_accepts_lambda_sparse -v`

Expected: FAIL with "unexpected keyword argument 'lambda_sparse'"

**Step 3: Add lambda_sparse to EM function signature**

Modify `CITEgeist/model/gurobi_impl.py` at the `optimize_discrete_cell_assignment_em` function (around line 3237):

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
    lambda_sparse: float = 0.0,  # NEW: sparsity penalty weight
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float], np.ndarray]:
```

**Step 4: Pass lambda_sparse to solve_discrete_cell_counts in E-step**

Find the E-step call (around line 3308-3320) and modify:

```python
        # ==================== E-Step ====================
        # Solve IQP for cell counts given current beta
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
            lambda_sparse=lambda_sparse,  # NEW: pass through
        )
```

**Step 5: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_sparsity.py -v`

Expected: PASS

**Step 6: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py CITEgeist/tests/test_discrete_sparsity.py
git commit -m "feat(discrete): propagate lambda_sparse through EM algorithm"
```

---

### Task 0.3: Add lambda_sparse to CitegeistModel.run_discrete_cell_assignment

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:912-1043`
- Test: `CITEgeist/tests/test_discrete_sparsity.py` (add test)

**Step 1: Write the failing test**

Add to `CITEgeist/tests/test_discrete_sparsity.py`:

```python
def test_model_accepts_lambda_sparse():
    """Verify CitegeistModel.run_discrete_cell_assignment accepts lambda_sparse."""
    import scanpy as sc
    from CITEgeist.model.citegeist_model import CitegeistModel

    # Create minimal AnnData objects
    n_spots, n_genes, n_markers = 10, 50, 4
    np.random.seed(42)

    gex_data = np.random.rand(n_spots, n_genes)
    adata_gex = sc.AnnData(gex_data)
    adata_gex.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata_gex.var_names = [f"gene_{i}" for i in range(n_genes)]

    cite_data = np.random.rand(n_spots, n_markers)
    adata_cite = sc.AnnData(cite_data)
    adata_cite.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata_cite.var_names = ["TypeA_Protein_1", "TypeA_Protein_2", "TypeB_Protein_1", "TypeB_Protein_2"]

    model = CitegeistModel(
        sample_name="test",
        output_folder="/tmp/citegeist_test",
        simulation=True,
        gene_expression_adata=adata_gex,
        antibody_capture_adata=adata_cite,
    )

    model.preprocess_antibody_discrete()
    model.load_cell_profile_dict({
        "TypeA": {"Major": ["TypeA_Protein_1", "TypeA_Protein_2"]},
        "TypeB": {"Major": ["TypeB_Protein_1", "TypeB_Protein_2"]},
    })

    # Create nuclei counts
    nuclei_counts = pd.Series([3] * n_spots, index=adata_gex.obs_names)

    # Should accept lambda_sparse parameter
    result = model.run_discrete_cell_assignment(
        nuclei_counts=nuclei_counts,
        lambda_sparse=0.1,  # New parameter
        max_em_iterations=2,
    )

    assert result is not None
    assert result.shape == (n_spots, 2)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_sparsity.py::test_model_accepts_lambda_sparse -v`

Expected: FAIL

**Step 3: Add lambda_sparse to run_discrete_cell_assignment**

Modify `CITEgeist/model/citegeist_model.py` at `run_discrete_cell_assignment` (around line 912):

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
        lambda_sparse: float = 0.0,  # NEW: sparsity penalty weight
    ) -> pd.DataFrame:
```

And update the call to `optimize_discrete_cell_assignment_em`:

```python
        c_values, beta_values, em_stats, alpha_values = optimize_discrete_cell_assignment_em(
            marker_level_data=marker_data,
            marker_names=marker_names,
            assignment_matrix=assignment_matrix,
            cell_type_names=cell_type_names,
            nuclei_counts=nuclei_arr,
            max_em_iterations=max_em_iterations,
            beta_convergence_tol=beta_convergence_tol,
            beta_min=beta_min,
            beta_max=beta_max,
            max_nuclei_cap=max_nuclei_cap,
            timeout_per_spot=timeout_per_spot,
            lambda_sparse=lambda_sparse,  # NEW
        )
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && module load gurobi/12.0.3 && python -m pytest CITEgeist/tests/test_discrete_sparsity.py::test_model_accepts_lambda_sparse -v`

Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/citegeist_model.py CITEgeist/tests/test_discrete_sparsity.py
git commit -m "feat(discrete): add lambda_sparse to CitegeistModel API"
```

---

### Task 0.4: Update Benchmark Script with lambda_sparse

**Files:**
- Modify: `Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py:366-373`

**Step 1: Add lambda_sparse parameter to benchmark script**

Modify the `run_discrete_cell_assignment` call around line 366:

```python
    # Run discrete assignment with sparsity
    start = time.time()
    cell_counts = model.run_discrete_cell_assignment(
        nuclei_counts=pred_aligned,
        max_em_iterations=max_em_iterations,
        beta_convergence_tol=1e-4,
        max_nuclei_cap=30,
        lambda_sparse=args.lambda_sparse,  # NEW: from CLI args
    )
    timings["discrete_assignment_sec"] = time.time() - start
```

**Step 2: Add CLI argument**

Add to the argument parser around line 496:

```python
    parser.add_argument("--lambda-sparse", type=float, default=0.1,
                        help="Sparsity penalty weight (default: 0.1)")
```

**Step 3: Update results to include lambda_sparse**

Add to results dict around line 271:

```python
    results = {
        "replicate_id": replicate_id,
        "condition": condition,
        "mode": mode,
        "lambda_sparse": args.lambda_sparse if hasattr(args, 'lambda_sparse') else 0.0,
    }
```

**Step 4: Commit**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py
git commit -m "feat(benchmark): add lambda_sparse parameter to discrete benchmark"
```

---

### Task 0.5: Run Benchmark Validation (SLURM Job)

**Files:**
- Create: `Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_sparsity_benchmark.sh`

**Step 1: Create SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=citegeist_sparsity
#SBATCH --output=logs/sparsity_%A_%a.out
#SBATCH --error=logs/sparsity_%A_%a.err
#SBATCH --array=0-9
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alc376@pitt.edu

module load gurobi/12.0.3
eval "$(conda shell.bash hook)"
conda activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

# Map array index to condition and replicate
CONDITIONS=("high_seg" "high_seg" "high_seg" "high_seg" "high_seg" "mixed" "mixed" "mixed" "mixed" "mixed")
REPLICATES=(0 1 2 3 4 0 1 2 3 4)

CONDITION=${CONDITIONS[$SLURM_ARRAY_TASK_ID]}
REPLICATE=${REPLICATES[$SLURM_ARRAY_TASK_ID]}

echo "Running: condition=$CONDITION, replicate=$REPLICATE, lambda_sparse=0.1"

cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist

python Benchmarking/simulation_benchmarking/CITEgeist/src/benchmark_discrete_simulation.py \
    --replicate-id $REPLICATE \
    --condition $CONDITION \
    --mode dapi \
    --output-dir Benchmarking/simulation_benchmarking/CITEgeist/output_sparsity \
    --lambda-sparse 0.1

echo "Done"
```

**Step 2: Submit and monitor**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/Benchmarking/simulation_benchmarking/CITEgeist/slurm
mkdir -p logs
sbatch sbatch_sparsity_benchmark.sh
```

**Step 3: Analyze results**

After jobs complete, compare proportion correlation:

```bash
python3 -c "
import json
from pathlib import Path

results_dir = Path('Benchmarking/simulation_benchmarking/CITEgeist/output_sparsity')
for cond in ['high_seg', 'mixed']:
    corrs = []
    for rep in range(5):
        f = results_dir / cond / 'dapi' / f'Wu_rep_{rep}' / 'benchmark_results.json'
        if f.exists():
            d = json.load(open(f))
            corrs.append(d['proportion_correlation'])
    if corrs:
        print(f'{cond}: mean={sum(corrs)/len(corrs):.4f}, n={len(corrs)}')
"
```

**Step 4: Commit results**

```bash
git add Benchmarking/simulation_benchmarking/CITEgeist/slurm/sbatch_sparsity_benchmark.sh
git commit -m "feat(benchmark): add SLURM script for sparsity benchmark validation"
```

---

## Phase 1-4: Single-Cell Pipeline (After Phase 0 Validation)

**Note:** Proceed with Phases 1-4 only after Phase 0 achieves proportion correlation ≥ 0.90.

### Task 1.1: Create Single-Cell Assembly Module

**Files:**
- Create: `CITEgeist/model/single_cell_assembly.py`
- Test: `CITEgeist/tests/test_single_cell_assembly.py`

(Detailed steps to be expanded after Phase 0 validation)

### Task 1.2: Implement Nucleus-Scaled GEX Optimization

**Files:**
- Create: `CITEgeist/model/nucleus_scaled_gex.py`
- Test: `CITEgeist/tests/test_nucleus_scaled_gex.py`

(Detailed steps to be expanded after Phase 0 validation)

### Task 1.3: Implement GEX-Informed Refinement Loop

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py`
- Test: `CITEgeist/tests/test_gex_refinement.py`

(Detailed steps to be expanded after Phase 0 validation)

### Task 1.4: Implement run_single_cell_deconvolution API

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py`
- Test: `CITEgeist/tests/test_single_cell_api.py`

(Detailed steps to be expanded after Phase 0 validation)

### Task 1.5: Add Deprecation Warnings

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py`
- Modify: `CITEgeist/model/analysis_functions.py`

(Detailed steps to be expanded after Phase 0 validation)

### Task 1.6: Final Integration Testing

**Files:**
- Create: `CITEgeist/tests/test_single_cell_integration.py`

(Detailed steps to be expanded after Phase 0 validation)

---

## Validation Checkpoint

After Phase 0 Tasks 0.1-0.5:

1. Run benchmark validation script
2. Compare proportion correlation to baseline:
   - high_seg: target ≥ 0.90 (baseline 0.78)
   - mixed: target ≥ 0.70 (baseline 0.44)
3. If targets met: proceed to Phase 1-4
4. If targets not met: investigate additional fixes before proceeding
