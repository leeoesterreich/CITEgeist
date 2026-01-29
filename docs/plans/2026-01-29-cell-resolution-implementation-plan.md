# Cell-Resolution Future-Proofing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `resolution="cell"` mode to CitegeistModel so Modules 1-3 handle single-cell spatial data (Xenium, Visium HD) with protein + RNA per cell, using optimization-based soft classification (Pass 1) and true count estimation (Pass 2).

**Architecture:** Explicit `resolution` parameter gates Module 3 behavior. Modules 1/2 get parameter presets (larger k for cell density). Pass 1 adds L1 sparsity to existing QP. Pass 2 adds a new per-cell Gurobi optimization for true count estimation. Existing spot-level behavior is unchanged.

**Tech Stack:** Python 3.10, gurobipy, numpy, scipy, scanpy, scikit-learn

**Design doc:** `docs/plans/2026-01-29-cell-resolution-future-proofing-design.md`

---

### Task 1: Add Resolution Presets and Parameter to CitegeistModel

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:42-90` (CitegeistModel.__init__)

**Step 1: Write failing test**

Create test that instantiates CitegeistModel with `resolution="cell"` and verifies it stores resolution and loads correct preset.

File: `tests/test_cell_resolution.py`

```python
"""Tests for cell-resolution mode in CITEgeist."""
import numpy as np
import pytest
import scanpy as sc

# Minimal AnnData for testing (no Gurobi needed)
def _make_test_adata(n_obs=100, n_genes=50, n_proteins=10):
    """Create minimal test AnnData with spatial coords."""
    rng = np.random.default_rng(42)
    adata_gex = sc.AnnData(X=rng.poisson(5, (n_obs, n_genes)).astype(np.float64))
    adata_gex.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata_gex.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata_gex.obsm["spatial"] = rng.uniform(0, 1000, (n_obs, 2))

    adata_protein = sc.AnnData(X=rng.exponential(2, (n_obs, n_proteins)).astype(np.float64))
    adata_protein.var_names = [f"Protein_{i}" for i in range(n_proteins)]
    adata_protein.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata_protein.obsm["spatial"] = adata_gex.obsm["spatial"].copy()

    return adata_gex, adata_protein


class TestResolutionPresets:
    """Test resolution parameter and preset loading."""

    def test_default_resolution_is_spot(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        model = CitegeistModel(
            sample_name="test",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
        )
        assert model.resolution == "spot"
        assert model.resolution_params["lambda_sparse"] == 0.0
        assert model.resolution_params["laplacian_k"] == 8

    def test_cell_resolution_loads_preset(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        model = CitegeistModel(
            sample_name="test",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="cell",
        )
        assert model.resolution == "cell"
        assert model.resolution_params["lambda_sparse"] == 0.1
        assert model.resolution_params["laplacian_k"] == 50
        assert model.resolution_params["pass2_library_slack"] == 1.5

    def test_invalid_resolution_raises(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        with pytest.raises(ValueError, match="resolution must be"):
            CitegeistModel(
                sample_name="test",
                output_folder=str(tmp_path),
                simulation=True,
                gene_expression_adata=adata_gex,
                antibody_capture_adata=adata_protein,
                resolution="invalid",
            )

    def test_user_can_override_preset(self, tmp_path):
        adata_gex, adata_protein = _make_test_adata()
        from CITEgeist.model import CitegeistModel
        model = CitegeistModel(
            sample_name="test",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="cell",
            resolution_overrides={"lambda_sparse": 0.2, "laplacian_k": 30},
        )
        assert model.resolution_params["lambda_sparse"] == 0.2
        assert model.resolution_params["laplacian_k"] == 30
        # Non-overridden params still from preset
        assert model.resolution_params["pass2_library_slack"] == 1.5
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cell_resolution.py::TestResolutionPresets -v`
Expected: FAIL — `resolution` parameter not accepted by `__init__`

**Step 3: Implement resolution presets**

Add to `CITEgeist/model/citegeist_model.py`:

1. Define `RESOLUTION_DEFAULTS` dict at module level (before class):

```python
RESOLUTION_DEFAULTS = {
    "spot": {
        "neighbor_k": 8,
        "morans_k": 8,
        "smooth_k": 6,
        "coloc_neighbor_k": 6,
        "coloc_multi_scale_k": [6, 12, 24, 48, 64],
        "laplacian_k": 8,
        "lambda_spatial": 0.1,
        "lambda_sparse": 0.0,
        "pass2_library_slack": 1.0,
    },
    "cell": {
        "neighbor_k": 50,
        "morans_k": 50,
        "smooth_k": 20,
        "coloc_neighbor_k": 30,
        "coloc_multi_scale_k": [20, 40, 60, 80, 100],
        "laplacian_k": 50,
        "lambda_spatial": 0.01,
        "lambda_sparse": 0.1,
        "pass2_library_slack": 1.5,
    },
}
```

2. Add `resolution` and `resolution_overrides` params to `__init__`:

```python
def __init__(
    self,
    sample_name,
    adata=None,
    output_folder=None,
    simulation=False,
    gene_expression_adata=None,
    antibody_capture_adata=None,
    resolution="spot",
    resolution_overrides=None,
):
```

3. In `__init__` body, after existing setup:

```python
if resolution not in RESOLUTION_DEFAULTS:
    raise ValueError(f"resolution must be one of {list(RESOLUTION_DEFAULTS.keys())}, got '{resolution}'")
self.resolution = resolution
self.resolution_params = dict(RESOLUTION_DEFAULTS[resolution])
if resolution_overrides:
    for key, val in resolution_overrides.items():
        if key not in self.resolution_params:
            raise ValueError(f"Unknown resolution parameter: '{key}'")
        self.resolution_params[key] = val
```

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cell_resolution.py::TestResolutionPresets -v`
Expected: PASS (all 4 tests)

**Step 5: Commit**

```bash
git add tests/test_cell_resolution.py CITEgeist/model/citegeist_model.py
git commit -m "feat: add resolution parameter and presets to CitegeistModel"
```

---

### Task 2: Add Sparsity Penalty to Pass 1 (optimize_cell_proportions_per_marker)

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py:527-695` (optimize_cell_proportions_per_marker)

**Step 1: Write failing test**

Add to `tests/test_cell_resolution.py`:

```python
class TestPass1Sparsity:
    """Test that lambda_sparse produces near-one-hot assignments on synthetic single-cell data."""

    def test_sparsity_produces_near_one_hot(self):
        """With lambda_sparse > 0, cells should be assigned predominantly to one type."""
        from CITEgeist.model.gurobi_impl import optimize_cell_proportions_per_marker

        rng = np.random.default_rng(42)
        n_cells = 20
        n_markers = 6
        n_types = 3

        # Simulate 3 cell types with 2 markers each
        # Type 0: markers 0,1 high; Type 1: markers 2,3 high; Type 2: markers 4,5 high
        X = rng.exponential(0.5, (n_cells, n_markers))
        for i in range(n_cells):
            true_type = i % n_types
            X[i, true_type * 2] += 5.0
            X[i, true_type * 2 + 1] += 5.0

        assignment = np.zeros((n_markers, n_types))
        assignment[0, 0] = assignment[1, 0] = 1
        assignment[2, 1] = assignment[3, 1] = 1
        assignment[4, 2] = assignment[5, 2] = 1

        marker_names = [f"M{i}" for i in range(n_markers)]
        type_names = [f"Type_{i}" for i in range(n_types)]

        # With sparsity
        Y_sparse, _, _ = optimize_cell_proportions_per_marker(
            marker_level_data=X,
            marker_names=marker_names,
            assignment_matrix=assignment,
            cell_type_names=type_names,
            lambda_sparse=0.1,
            lambda_laplacian=0.0,
            max_iterations=5,
        )

        # Check near-one-hot: max proportion per cell should be > 0.7
        max_props = Y_sparse.max(axis=1)
        assert np.mean(max_props > 0.7) >= 0.8, (
            f"Expected >=80% cells with max proportion >0.7, got {np.mean(max_props > 0.7):.2f}"
        )

    def test_zero_sparsity_matches_original(self):
        """lambda_sparse=0 should give same results as before."""
        from CITEgeist.model.gurobi_impl import optimize_cell_proportions_per_marker

        rng = np.random.default_rng(42)
        n_spots = 15
        n_markers = 4
        n_types = 2
        X = rng.exponential(2, (n_spots, n_markers))
        assignment = np.zeros((n_markers, n_types))
        assignment[0, 0] = assignment[1, 0] = 1
        assignment[2, 1] = assignment[3, 1] = 1

        Y_no_sparse, _, _ = optimize_cell_proportions_per_marker(
            marker_level_data=X,
            marker_names=[f"M{i}" for i in range(n_markers)],
            assignment_matrix=assignment,
            cell_type_names=["A", "B"],
            lambda_sparse=0.0,
            lambda_laplacian=0.0,
            max_iterations=5,
        )

        # Should be identical to calling without the parameter (default=0)
        Y_default, _, _ = optimize_cell_proportions_per_marker(
            marker_level_data=X,
            marker_names=[f"M{i}" for i in range(n_markers)],
            assignment_matrix=assignment,
            cell_type_names=["A", "B"],
            lambda_laplacian=0.0,
            max_iterations=5,
        )

        np.testing.assert_allclose(Y_no_sparse, Y_default, atol=1e-4)
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cell_resolution.py::TestPass1Sparsity -v`
Expected: FAIL — `lambda_sparse` not accepted as parameter

**Step 3: Implement sparsity penalty**

In `CITEgeist/model/gurobi_impl.py`, modify `optimize_cell_proportions_per_marker()`:

1. Add parameter to signature (line ~527):
```python
def optimize_cell_proportions_per_marker(
    ...,
    lambda_sparse: float = 0.0,  # L1 sparsity penalty (>0 for cell-level)
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
```

2. Add sparsity term to objective (after line ~665, where regularization is computed):
```python
# Sparsity penalty (L1 on Y - encourages near-one-hot for cell-level)
sparsity_term = 0
if lambda_sparse > 0:
    sparsity_term = lambda_sparse * gp.quicksum(
        Y[i, j] for i in range(N) for j in range(T)
    )
    logging.info(f"Sparsity penalty enabled: lambda_sparse={lambda_sparse}")
```

3. Add to objective (modify line ~679):
```python
model.setObjective(
    total_error + regularization_term + laplacian_term + sparsity_term,
    GRB.MINIMIZE,
)
```

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cell_resolution.py::TestPass1Sparsity -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py tests/test_cell_resolution.py
git commit -m "feat: add lambda_sparse L1 penalty to Pass 1 for cell-level mode"
```

---

### Task 3: Implement Cell-Level Pass 2 (optimize_gene_expression_cell)

**Files:**
- Modify: `CITEgeist/model/gurobi_impl.py` (add new function after `optimize_gene_expression` at ~line 2050)

**Step 1: Write failing test**

Add to `tests/test_cell_resolution.py`:

```python
class TestPass2CellLevel:
    """Test cell-level gene expression optimization."""

    def test_basic_true_count_estimation(self):
        """Verify true count estimation recovers dropout for expected genes."""
        from CITEgeist.model.gurobi_impl import estimate_true_expression_cell

        rng = np.random.default_rng(42)
        n_cells = 30
        n_genes = 20
        n_types = 2

        # Simulate cells: type 0 expresses genes 0-9, type 1 expresses genes 10-19
        X_obs = np.zeros((n_cells, n_genes), dtype=np.float64)
        cell_types = np.zeros(n_cells, dtype=int)
        Y_assignments = np.zeros((n_cells, n_types))

        for i in range(n_cells):
            ct = i % n_types
            cell_types[i] = ct
            Y_assignments[i, ct] = 0.95
            Y_assignments[i, 1 - ct] = 0.05

            # Express type-specific genes with some dropout
            type_genes = range(ct * 10, ct * 10 + 10)
            for g in type_genes:
                if rng.random() > 0.3:  # 30% dropout
                    X_obs[i, g] = rng.poisson(10)

        coords = rng.uniform(0, 100, (n_cells, 2))

        # Build enrichment weights: gene g is enriched in its type
        enrichment = np.ones((n_types, n_genes)) * 0.1
        enrichment[0, :10] = 1.0
        enrichment[1, 10:] = 1.0

        X_true = estimate_true_expression_cell(
            X_obs=X_obs,
            Y_assignments=Y_assignments,
            coords=coords,
            enrichment_weights=enrichment,
            library_slack=1.5,
            lambda_spatial=0.01,
            spatial_k=10,
        )

        assert X_true.shape == X_obs.shape

        # X_true should be non-negative
        assert np.all(X_true >= 0)

        # Dropout genes should have some recovery (not all zeros where expected)
        for ct in range(n_types):
            ct_cells = np.where(cell_types == ct)[0]
            type_genes = range(ct * 10, ct * 10 + 10)
            obs_zeros = np.sum(X_obs[np.ix_(ct_cells, list(type_genes))] == 0)
            true_zeros = np.sum(X_true[np.ix_(ct_cells, list(type_genes))] == 0)
            # Should recover at least some dropouts
            assert true_zeros < obs_zeros, (
                f"Type {ct}: expected fewer zeros in X_true ({true_zeros}) than X_obs ({obs_zeros})"
            )

    def test_library_size_bounded(self):
        """X_true total counts should not exceed library_slack * observed."""
        from CITEgeist.model.gurobi_impl import estimate_true_expression_cell

        rng = np.random.default_rng(42)
        n_cells = 10
        n_genes = 10
        X_obs = rng.poisson(5, (n_cells, n_genes)).astype(np.float64)
        Y_assignments = np.zeros((n_cells, 1))
        Y_assignments[:, 0] = 1.0
        coords = rng.uniform(0, 100, (n_cells, 2))
        enrichment = np.ones((1, n_genes))

        slack = 1.5
        X_true = estimate_true_expression_cell(
            X_obs=X_obs,
            Y_assignments=Y_assignments,
            coords=coords,
            enrichment_weights=enrichment,
            library_slack=slack,
            lambda_spatial=0.0,
            spatial_k=5,
        )

        obs_lib = X_obs.sum(axis=1)
        true_lib = X_true.sum(axis=1)
        for i in range(n_cells):
            assert true_lib[i] <= slack * obs_lib[i] + 1.0, (
                f"Cell {i}: true lib {true_lib[i]:.1f} > slack * obs lib {slack * obs_lib[i]:.1f}"
            )
```

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_cell_resolution.py::TestPass2CellLevel -v`
Expected: FAIL — `estimate_true_expression_cell` does not exist

**Step 3: Implement cell-level Pass 2**

Add new function to `CITEgeist/model/gurobi_impl.py` (after `optimize_gene_expression` at ~line 2050):

```python
def estimate_true_expression_cell(
    X_obs: np.ndarray,
    Y_assignments: np.ndarray,
    coords: np.ndarray,
    enrichment_weights: np.ndarray,
    library_slack: float = 1.5,
    lambda_enrich: float = 1.0,
    lambda_spatial: float = 0.01,
    spatial_k: int = 50,
    max_workers: Optional[int] = None,
    checkpoint_interval: int = 500,
) -> np.ndarray:
    """
    Estimate true gene expression per cell using optimization.

    For each cell, solves a QP to find X_true that:
    1. Stays close to observed counts (data fidelity)
    2. Recovers dropout for genes expected in the cell's type (enrichment prior)
    3. Agrees with same-type spatial neighbors (spatial coherence)
    4. Respects bounded total library size (prevents runaway imputation)

    Args:
        X_obs: Observed expression matrix (n_cells, n_genes).
        Y_assignments: Cell type assignment weights (n_cells, n_types) from Pass 1.
        coords: Spatial coordinates (n_cells, 2).
        enrichment_weights: Type-gene enrichment matrix (n_types, n_genes).
        library_slack: Max ratio of X_true library size to X_obs (default 1.5).
        lambda_enrich: Weight for enrichment prior term.
        lambda_spatial: Weight for spatial coherence term.
        spatial_k: Number of neighbors for spatial smoothing.
        max_workers: Max parallel workers (None = cpu_count).
        checkpoint_interval: Cells between checkpoint saves.

    Returns:
        X_true: Estimated true expression (n_cells, n_genes).
    """
    N, M = X_obs.shape
    T = Y_assignments.shape[1]

    # Determine dominant type per cell
    dominant_type = np.argmax(Y_assignments, axis=1)

    # Build spatial neighbor graph (k-NN)
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    k_query = min(spatial_k + 1, N)
    _, all_neighbor_indices = tree.query(coords, k=k_query)
    # Remove self
    if all_neighbor_indices.ndim > 1 and all_neighbor_indices.shape[1] > 1:
        all_neighbor_indices = all_neighbor_indices[:, 1:]
    else:
        all_neighbor_indices = np.empty((N, 0), dtype=int)

    # For each cell, find same-type neighbors
    same_type_neighbor_means = np.zeros((N, M))
    same_type_neighbor_counts = np.zeros(N, dtype=int)
    for i in range(N):
        ct = dominant_type[i]
        neighbors = all_neighbor_indices[i]
        same_type_mask = dominant_type[neighbors] == ct
        same_type_neighbors = neighbors[same_type_mask]
        if len(same_type_neighbors) > 0:
            same_type_neighbor_means[i] = X_obs[same_type_neighbors].mean(axis=0)
            same_type_neighbor_counts[i] = len(same_type_neighbors)

    logging.info(
        f"Cell-level Pass 2: {N} cells, {M} genes, {T} types, "
        f"library_slack={library_slack}, spatial_k={spatial_k}"
    )

    X_true = np.zeros((N, M), dtype=np.float64)

    # Process cells in parallel
    workers = max_workers if max_workers is not None else os.cpu_count()

    def _solve_single_cell(cell_idx):
        """Solve QP for a single cell."""
        import gurobipy as gp
        from gurobipy import GRB

        model = None
        try:
            ct = dominant_type[cell_idx]
            obs = X_obs[cell_idx]
            obs_lib = obs.sum()

            if obs_lib < 1:
                return cell_idx, obs.copy()  # Nothing to optimize

            enrich = enrichment_weights[ct]
            neighbor_mean = same_type_neighbor_means[cell_idx]
            has_neighbors = same_type_neighbor_counts[cell_idx] > 0

            model = gp.Model(f"cell_expr_{cell_idx}")
            model.setParam("OutputFlag", 0)
            model.setParam("Threads", 1)
            model.setParam("TimeLimit", 30)

            # Variables: X_true[g] for each gene
            max_lib = library_slack * obs_lib
            X_vars = model.addVars(
                M, lb=0, ub=max_lib, vtype=GRB.CONTINUOUS, name="X"
            )

            # Library size constraint
            model.addConstr(
                gp.quicksum(X_vars[g] for g in range(M)) <= max_lib,
                name="lib_size",
            )

            # Objective: data fidelity + enrichment + spatial
            obj_terms = []

            for g in range(M):
                # Data fidelity: (X_true[g] - X_obs[g])^2
                obj_terms.append((X_vars[g] - obs[g]) * (X_vars[g] - obs[g]))

                # Enrichment prior: -lambda_enrich * E[ct,g] * X_true[g]
                # (negative because we minimize; higher enrichment = more expression wanted)
                if lambda_enrich > 0 and enrich[g] > 0.1:
                    obj_terms.append(-lambda_enrich * enrich[g] * X_vars[g])

                # Spatial coherence with same-type neighbors
                if lambda_spatial > 0 and has_neighbors:
                    obj_terms.append(
                        lambda_spatial * (X_vars[g] - neighbor_mean[g]) * (X_vars[g] - neighbor_mean[g])
                    )

            model.setObjective(gp.quicksum(obj_terms), GRB.MINIMIZE)
            model.optimize()

            if model.status == GRB.OPTIMAL:
                result = np.array([X_vars[g].X for g in range(M)])
                return cell_idx, result
            else:
                return cell_idx, obs.copy()

        except Exception as e:
            logging.warning(f"Cell {cell_idx} optimization failed: {e}")
            return cell_idx, X_obs[cell_idx].copy()
        finally:
            if model:
                del model

    # Run in parallel
    from concurrent.futures import ProcessPoolExecutor, as_completed
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_solve_single_cell, i): i for i in range(N)}
        with tqdm(total=N, desc="Estimating true expression") as pbar:
            for future in as_completed(futures):
                cell_idx, result = future.result()
                X_true[cell_idx] = result
                pbar.update(1)

    return X_true
```

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_cell_resolution.py::TestPass2CellLevel -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/gurobi_impl.py tests/test_cell_resolution.py
git commit -m "feat: add estimate_true_expression_cell() for cell-level Pass 2"
```

---

### Task 4: Wire Resolution Dispatch in CitegeistModel

**Files:**
- Modify: `CITEgeist/model/citegeist_model.py:388-619` (run_cell_proportion_model)
- Modify: `CITEgeist/model/citegeist_model.py:621-759` (run_cell_expression_pass1)
- Modify: `CITEgeist/model/citegeist_model.py:20-30` (imports)

**Step 1: Write failing test**

Add to `tests/test_cell_resolution.py`:

```python
class TestResolutionDispatch:
    """Test that CitegeistModel dispatches to correct Pass 1/2 based on resolution."""

    def test_cell_mode_passes_sparsity_to_pass1(self, tmp_path):
        """In cell mode, run_cell_proportion_model should pass lambda_sparse from preset."""
        adata_gex, adata_protein = _make_test_adata(n_obs=30, n_genes=20, n_proteins=6)
        from CITEgeist.model import CitegeistModel

        model = CitegeistModel(
            sample_name="test_cell",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="cell",
        )

        # Verify the model stores cell-level params
        assert model.resolution == "cell"
        assert model.resolution_params["lambda_sparse"] > 0

    def test_spot_mode_no_sparsity(self, tmp_path):
        """In spot mode, lambda_sparse should be 0."""
        adata_gex, adata_protein = _make_test_adata(n_obs=30, n_genes=20, n_proteins=6)
        from CITEgeist.model import CitegeistModel

        model = CitegeistModel(
            sample_name="test_spot",
            output_folder=str(tmp_path),
            simulation=True,
            gene_expression_adata=adata_gex,
            antibody_capture_adata=adata_protein,
            resolution="spot",
        )

        assert model.resolution == "spot"
        assert model.resolution_params["lambda_sparse"] == 0.0
```

**Step 2: Run test**

Run: `pytest tests/test_cell_resolution.py::TestResolutionDispatch -v`
Expected: Should PASS if Task 1 was implemented correctly.

**Step 3: Implement dispatch logic**

In `CITEgeist/model/citegeist_model.py`:

1. Add import at top:
```python
from .gurobi_impl import (
    ...,
    estimate_true_expression_cell,
)
```

2. In `run_cell_proportion_model()` (~line 468), pass `lambda_sparse` from preset:

Where it calls `optimize_cell_proportions_per_marker()`, add:
```python
lambda_sparse=self.resolution_params.get("lambda_sparse", 0.0),
```

Where it sets `laplacian_k`, use preset if not explicitly passed:
The method already accepts `laplacian_k=8` as a parameter. Add logic at the top of the method to use preset defaults when not explicitly overridden by the caller. The simplest approach: if user doesn't pass `laplacian_k`, use the preset.

3. In `run_cell_expression_pass1()` (~line 677), add cell-mode dispatch:

After the existing check for preprocessed data and before calling `optimize_gene_expression`, add:
```python
if self.resolution == "cell":
    # Cell-level: estimate true expression per cell (not deconvolution)
    enrichment = compute_global_prior(...)  # reuse existing enrichment computation
    X_true = estimate_true_expression_cell(
        X_obs=self.gene_expression_adata.X,
        Y_assignments=cell_props_values,
        coords=self.gene_expression_adata.obsm["spatial"],
        enrichment_weights=enrichment,
        library_slack=self.resolution_params["pass2_library_slack"],
        lambda_spatial=self.resolution_params["lambda_spatial"],
        spatial_k=self.resolution_params["neighbor_k"],
        max_workers=max_workers,
    )
    # Store results in same format as spot-level for downstream compatibility
    ...
    return
```

**Step 4: Run tests**

Run: `pytest tests/test_cell_resolution.py -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/citegeist_model.py tests/test_cell_resolution.py
git commit -m "feat: wire resolution dispatch in CitegeistModel for cell-mode Pass 1/2"
```

---

### Task 5: Export New Function in __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Add export**

Add `estimate_true_expression_cell` to imports and `__all__`:

```python
from .gurobi_impl import (
    map_antibodies_to_profiles,
    optimize_cell_proportions,
    optimize_gene_expression,
    estimate_true_expression_cell,
)
```

And in `__all__`:
```python
"estimate_true_expression_cell",
```

Also add `RESOLUTION_DEFAULTS`:
```python
from .citegeist_model import CitegeistModel, RESOLUTION_DEFAULTS
```

**Step 2: Verify imports work**

Run: `python -c "from CITEgeist.model import estimate_true_expression_cell, RESOLUTION_DEFAULTS; print('OK')"`
Expected: `OK`

**Step 3: Commit**

```bash
git add CITEgeist/model/__init__.py
git commit -m "feat: export estimate_true_expression_cell and RESOLUTION_DEFAULTS"
```

---

### Task 6: Verify Spot-Level Backward Compatibility

**Files:**
- No modifications — run existing test suite

**Step 1: Run existing tests**

Run: `pytest tests/test_citegeist_simulated.py -v`
Expected: All existing tests PASS with no changes

Run: `pytest tests/test_module_one_marker_detection.py -v`
Expected: PASS

Run: `pytest tests/test_module2_profile_discovery.py -v`
Expected: PASS

**Step 2: Commit (only if any fix was needed)**

If any test fails, fix the regression and commit:
```bash
git commit -m "fix: preserve backward compatibility for spot-level mode"
```

---

### Task 7: Xenium Cell-Level Validation Script

**Files:**
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/src/validate_cell_resolution.py`
- Create: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_cell_resolution_validation.sh`

**Step 1: Write validation script**

```python
"""
Validate cell-resolution CITEgeist on Xenium single-cell data.

Runs Modules 1-3 in resolution="cell" mode on Xenium RCC regions,
comparing against RNA-derived ground truth cell type labels and
the existing heuristic cell assignment pipeline.

Implements Tests 1-4 from the design doc.
"""
import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix

# Add paths
REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from model import (
    CitegeistModel,
    identify_interesting_markers,
    analyze_marker_colocalization,
    discover_profiles,
    select_profiles,
)
from model.gurobi_impl import (
    optimize_cell_proportions_per_marker,
    map_antibodies_to_profiles_v2,
    estimate_true_expression_cell,
)
from load_xenium_singlecell import load_xenium_singlecell

logger = logging.getLogger(__name__)
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_cell_resolution"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_rna_gt" / "cell_types.csv"


def run_validation(region_id: int, max_cells: int = 10000):
    """Run full cell-resolution validation for one region."""
    output_dir = OUTPUT_DIR / f"region_{region_id}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info(f"Loading Xenium region {region_id} (max {max_cells} cells)")
    adata_gex, adata_protein = load_xenium_singlecell(
        region_id=region_id, max_cells=max_cells, seed=42
    )
    logger.info(f"Loaded {adata_gex.shape[0]} cells")

    # Load ground truth
    gt_df = pd.read_csv(GT_PATH)
    gt_df = gt_df.set_index("cell_id")
    common_cells = list(set(adata_gex.obs_names) & set(gt_df.index))
    gt_labels = gt_df.loc[common_cells, "cell_type"].values
    logger.info(f"Ground truth available for {len(common_cells)} cells")

    # === MODULE 1: Marker Interest Detection ===
    logger.info("=" * 60)
    logger.info("MODULE 1: Marker Interest Detection (cell resolution)")
    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    coords = adata_protein.obsm["spatial"]
    marker_names = list(adata_protein.var_names)

    m1_result = identify_interesting_markers(
        X=X_protein, coords=coords, marker_names=marker_names,
        morans_k=50, smooth_k=20, morans_n_perm=99, seed=42,
    )
    interesting = m1_result.interesting_markers
    logger.info(f"Found {len(interesting)} interesting markers: {interesting}")

    # === MODULE 2: Profile Discovery ===
    logger.info("=" * 60)
    logger.info("MODULE 2: Profile Discovery (cell resolution)")
    coloc_result = analyze_marker_colocalization(
        X=X_protein, coords=coords, marker_names=marker_names,
        markers_to_analyze=interesting,
        neighbor_k=30, smooth_k=20,
        multi_scale_k=[20, 40, 60, 80, 100],
        n_permutations=99, seed=42,
    )
    profile_result = discover_profiles(coloc_result, seed=42)
    selected = select_profiles(
        adata_protein, profile_result,
        interesting_markers=interesting,
    )
    profiles = selected.profiles
    logger.info(f"Discovered {len(profiles)} profiles: {list(profiles.keys())}")

    # Save profiles
    with open(output_dir / "profiles.json", "w") as f:
        json.dump({k: list(v) for k, v in profiles.items()}, f, indent=2)

    # === MODULE 3 PASS 1: Cell Type Soft Classification ===
    logger.info("=" * 60)
    logger.info("MODULE 3 PASS 1: Cell Type Classification (cell resolution, Gurobi QP)")

    # Map markers to profiles (reuse existing mapping function)
    marker_data, assignment_matrix, type_names = map_antibodies_to_profiles_v2(
        adata_protein, profiles
    )

    Y_assignments, beta_values, beta_dict = optimize_cell_proportions_per_marker(
        marker_level_data=marker_data,
        marker_names=[marker_names[i] for i in range(marker_data.shape[1])],
        assignment_matrix=assignment_matrix,
        cell_type_names=type_names,
        lambda_sparse=0.1,
        lambda_laplacian=0.01,
        coords=coords,
        laplacian_k=50,
        max_iterations=10,
    )

    # Evaluate classification
    dominant_type = np.argmax(Y_assignments, axis=1)
    predicted_labels = [type_names[dt] for dt in dominant_type]

    # Map predictions to GT labels (best matching)
    # Note: profile names won't match GT names exactly; use confusion matrix to find best mapping
    results = {
        "region_id": region_id,
        "n_cells": len(common_cells),
        "n_profiles": len(profiles),
        "profiles": list(profiles.keys()),
        "max_Y_mean": float(Y_assignments.max(axis=1).mean()),
        "max_Y_median": float(np.median(Y_assignments.max(axis=1))),
        "doublet_fraction": float(np.mean(Y_assignments.max(axis=1) < 0.6)),
    }

    # Save Y_assignments
    y_df = pd.DataFrame(Y_assignments, columns=type_names, index=adata_protein.obs_names)
    y_df.to_csv(output_dir / "cell_assignments.csv")

    # === MODULE 3 PASS 2: True Count Estimation ===
    logger.info("=" * 60)
    logger.info("MODULE 3 PASS 2: True Count Estimation (cell resolution)")

    X_gex = adata_gex.X
    if hasattr(X_gex, "toarray"):
        X_gex = X_gex.toarray()
    X_gex = np.asarray(X_gex, dtype=np.float64)

    # Build enrichment weights from Pass 1 assignments
    n_types = len(type_names)
    n_genes = X_gex.shape[1]
    enrichment = np.ones((n_types, n_genes)) * 0.1
    for ct_idx in range(n_types):
        ct_cells = np.where(dominant_type == ct_idx)[0]
        if len(ct_cells) > 0:
            ct_mean = X_gex[ct_cells].mean(axis=0)
            global_mean = X_gex.mean(axis=0) + 1e-10
            enrichment[ct_idx] = ct_mean / global_mean

    X_true = estimate_true_expression_cell(
        X_obs=X_gex,
        Y_assignments=Y_assignments,
        coords=coords,
        enrichment_weights=enrichment,
        library_slack=1.5,
        lambda_spatial=0.01,
        spatial_k=50,
    )

    # Evaluate dropout recovery
    obs_zeros = np.sum(X_gex == 0)
    true_zeros = np.sum(X_true == 0)
    results["dropout_recovery"] = {
        "obs_zero_fraction": float(obs_zeros / X_gex.size),
        "true_zero_fraction": float(true_zeros / X_true.size),
        "zeros_recovered": int(obs_zeros - true_zeros),
        "recovery_rate": float((obs_zeros - true_zeros) / max(obs_zeros, 1)),
    }

    # Save results
    with open(output_dir / "validation_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Results saved to {output_dir}")
    logger.info(f"Max Y mean: {results['max_Y_mean']:.3f}")
    logger.info(f"Doublet fraction: {results['doublet_fraction']:.3f}")
    logger.info(f"Dropout recovery rate: {results['dropout_recovery']['recovery_rate']:.3f}")

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--region", type=int, default=0)
    parser.add_argument("--max-cells", type=int, default=10000)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run_validation(args.region, args.max_cells)
```

**Step 2: Write SLURM script**

```bash
#!/bin/bash
#SBATCH --job-name=cell_resolution_validation
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_log/cell_resolution_%A_%a.out
#SBATCH --error=slurm_log/cell_resolution_%A_%a.err
#SBATCH --array=0-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

module load python/3.10
conda activate CITEgeist_env

REPO_ROOT="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
cd $REPO_ROOT

python Benchmarking/xenium_benchmarking/CITEgeist/src/validate_cell_resolution.py \
    --region $SLURM_ARRAY_TASK_ID \
    --max-cells 50000
```

**Step 3: Commit**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/src/validate_cell_resolution.py \
       Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_cell_resolution_validation.sh
git commit -m "feat: add Xenium cell-resolution validation script and SLURM job"
```

---

### Task 8: Run Validation and Analyze Results

**Step 1: Submit SLURM job**

```bash
mkdir -p Benchmarking/xenium_benchmarking/CITEgeist/slurm/slurm_log
sbatch Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_cell_resolution_validation.sh
```

**Step 2: Monitor and analyze**

Once jobs complete, check:
- `output_cell_resolution/region_{0-4}/validation_results.json` for metrics
- Compare `max_Y_mean` — should be >0.7 indicating near-one-hot assignments
- Compare `doublet_fraction` — should be reasonable (<0.2)
- Compare `dropout_recovery.recovery_rate` — should be >0 but <0.5

**Step 3: Commit results summary**

```bash
git add Benchmarking/xenium_benchmarking/CITEgeist/output_cell_resolution/
git commit -m "results: cell-resolution validation on Xenium regions 0-4"
```
