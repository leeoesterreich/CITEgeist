# Joint Program Discovery Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement `discover_joint_programs()` function that discovers spatial programs across all cell types simultaneously, enabling detection of cross-cell-type interaction programs.

**Architecture:** Stack deconvolved gene expression layers from all cell types, run single NMF on combined matrix, then assign cell type labels post-hoc based on correlation with cell type proportions. This complements the existing per-anchor discovery by finding programs that span multiple cell types.

**Tech Stack:** Python 3.10, numpy, sklearn.decomposition.NMF, scanpy, scipy

---

## Task 1: Create JointProgramResult Data Structure

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py:45-150` (add after existing dataclasses)

**Step 1: Write the test for the new dataclass**

Create test file:
```python
# tests/test_module4_joint.py
"""Test Module 4 Joint Program Discovery."""

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))


def test_joint_program_dataclass():
    """Test JointProgram dataclass initialization."""
    from CITEgeist.model.anchored_program_discovery import JointProgram

    program = JointProgram(
        program_id=0,
        top_genes=["GENE1", "GENE2", "GENE3"],
        gene_loadings={"GENE1": 0.8, "GENE2": 0.5, "GENE3": 0.3},
        variance_explained=0.15,
        spatial_moran_i=0.45,
        spatial_moran_pvalue=0.001,
        mean_activity=0.25,
        active_spots_fraction=0.6,
        cell_type_enrichments={"Cancer": 0.7, "Macrophage": 0.3},
        primary_cell_type="Cancer",
        secondary_cell_type="Macrophage",
        interaction_score=0.35,
        program_type="interaction",
    )

    assert program.program_id == 0
    assert program.primary_cell_type == "Cancer"
    assert program.interaction_score == 0.35
    assert program.program_type == "interaction"
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_joint_program_dataclass -v`
Expected: FAIL with "cannot import name 'JointProgram'"

**Step 3: Implement JointProgram dataclass**

Add to `CITEgeist/model/anchored_program_discovery.py` after line ~98 (after SpatialProgram):

```python
@dataclass
class JointProgram:
    """A spatial program discovered from joint analysis of all cell types."""

    program_id: int
    """Unique identifier for this program."""

    top_genes: List[str]
    """Top N genes by loading magnitude."""

    gene_loadings: Dict[str, float]
    """Gene name -> loading value for top genes."""

    variance_explained: float
    """Fraction of variance explained by this program."""

    spatial_moran_i: float
    """Moran's I spatial autocorrelation of program activity."""

    spatial_moran_pvalue: float
    """P-value for Moran's I test."""

    mean_activity: float
    """Mean program activity across all spots."""

    active_spots_fraction: float
    """Fraction of spots with above-median activity."""

    cell_type_enrichments: Dict[str, float]
    """Cell type -> enrichment score (correlation with proportions)."""

    primary_cell_type: str
    """Cell type with highest enrichment."""

    secondary_cell_type: Optional[str]
    """Second highest cell type if interaction program."""

    interaction_score: float
    """Score 0-1 indicating how multi-cell-type the program is (0=single, 1=balanced)."""

    program_type: str
    """'single_celltype', 'interaction', or 'microenvironment'."""
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_joint_program_dataclass -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_joint.py
git commit -m "feat(module4): add JointProgram dataclass for cross-cell-type programs"
```

---

## Task 2: Create JointDiscoveryResult Data Structure

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py` (add after JointProgram)
- Modify: `tests/test_module4_joint.py`

**Step 1: Write the test**

Add to `tests/test_module4_joint.py`:

```python
def test_joint_discovery_result_dataclass():
    """Test JointDiscoveryResult dataclass initialization."""
    from CITEgeist.model.anchored_program_discovery import JointProgram, JointDiscoveryResult

    program = JointProgram(
        program_id=0,
        top_genes=["GENE1"],
        gene_loadings={"GENE1": 0.8},
        variance_explained=0.15,
        spatial_moran_i=0.45,
        spatial_moran_pvalue=0.001,
        mean_activity=0.25,
        active_spots_fraction=0.6,
        cell_type_enrichments={"Cancer": 0.7},
        primary_cell_type="Cancer",
        secondary_cell_type=None,
        interaction_score=0.0,
        program_type="single_celltype",
    )

    result = JointDiscoveryResult(
        programs=[program],
        W=np.array([[0.8], [0.5]]),
        H=np.array([[0.3, 0.4, 0.5]]),
        gene_names=["GENE1", "GENE2"],
        cell_type_names=["Cancer", "Macrophage"],
        n_spots=3,
        reconstruction_error=0.05,
        parameters={"K_programs": 1},
    )

    assert len(result.programs) == 1
    assert result.n_spots == 3
    assert result.W.shape == (2, 1)

    # Test to_dataframe method
    df = result.to_dataframe()
    assert len(df) == 1
    assert "primary_cell_type" in df.columns
    assert "interaction_score" in df.columns
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_joint_discovery_result_dataclass -v`
Expected: FAIL with "cannot import name 'JointDiscoveryResult'"

**Step 3: Implement JointDiscoveryResult dataclass**

Add to `CITEgeist/model/anchored_program_discovery.py` after JointProgram:

```python
@dataclass
class JointDiscoveryResult:
    """Results from joint program discovery across all cell types."""

    programs: List[JointProgram]
    """List of discovered programs with cell type assignments."""

    W: NDArray[np.floating]
    """Gene loadings matrix (n_genes, K_programs)."""

    H: NDArray[np.floating]
    """Program activities matrix (K_programs, n_spots)."""

    gene_names: List[str]
    """Gene names corresponding to W rows."""

    cell_type_names: List[str]
    """Cell types included in analysis."""

    n_spots: int
    """Number of spots analyzed."""

    reconstruction_error: float
    """NMF reconstruction error."""

    H_by_celltype: Optional[Dict[str, NDArray]] = None
    """Program activities split by cell type (from unstacking)."""

    parameters: Dict[str, Any] = field(default_factory=dict)
    """Parameters used for discovery."""

    def summary(self) -> str:
        """Return summary string."""
        n_single = sum(1 for p in self.programs if p.program_type == "single_celltype")
        n_interaction = sum(1 for p in self.programs if p.program_type == "interaction")
        n_micro = sum(1 for p in self.programs if p.program_type == "microenvironment")

        lines = [
            "Joint Program Discovery Results",
            "=" * 40,
            f"Total programs: {len(self.programs)}",
            f"  Single cell-type: {n_single}",
            f"  Interaction: {n_interaction}",
            f"  Microenvironment: {n_micro}",
            f"Spots: {self.n_spots}",
            f"Genes: {len(self.gene_names)}",
            f"Cell types: {', '.join(self.cell_type_names)}",
            f"Reconstruction error: {self.reconstruction_error:.4f}",
        ]
        return "\n".join(lines)

    def to_dataframe(self) -> pd.DataFrame:
        """Convert programs to DataFrame."""
        records = []
        for prog in self.programs:
            records.append({
                "program_id": prog.program_id,
                "top_genes": ", ".join(prog.top_genes[:10]),
                "variance_explained": prog.variance_explained,
                "spatial_moran_i": prog.spatial_moran_i,
                "spatial_moran_pvalue": prog.spatial_moran_pvalue,
                "mean_activity": prog.mean_activity,
                "active_spots_fraction": prog.active_spots_fraction,
                "primary_cell_type": prog.primary_cell_type,
                "secondary_cell_type": prog.secondary_cell_type,
                "interaction_score": prog.interaction_score,
                "program_type": prog.program_type,
            })
        return pd.DataFrame(records)
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_joint_discovery_result_dataclass -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_joint.py
git commit -m "feat(module4): add JointDiscoveryResult dataclass"
```

---

## Task 3: Implement Cell Type Assignment Function

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py` (add helper function)
- Modify: `tests/test_module4_joint.py`

**Step 1: Write the test**

Add to `tests/test_module4_joint.py`:

```python
def test_assign_program_cell_types():
    """Test cell type assignment for joint programs."""
    from CITEgeist.model.anchored_program_discovery import _assign_program_cell_types

    # Mock H matrix: 2 programs, 100 spots
    np.random.seed(42)
    H = np.random.rand(2, 100)

    # Mock proportions: program 0 correlates with Cancer, program 1 with both
    proportions = pd.DataFrame({
        "Cancer": np.linspace(0, 1, 100) + np.random.normal(0, 0.1, 100),
        "Macrophage": np.linspace(1, 0, 100) + np.random.normal(0, 0.1, 100),
    })
    # Make program 0 correlate with Cancer
    H[0, :] = proportions["Cancer"].values + np.random.normal(0, 0.1, 100)
    H[0, H[0, :] < 0] = 0
    # Make program 1 correlate with both equally
    H[1, :] = 0.5 * proportions["Cancer"].values + 0.5 * proportions["Macrophage"].values

    result = _assign_program_cell_types(H, proportions)

    # Program 0 should be single cell type (Cancer)
    assert result[0]["primary_cell_type"] == "Cancer"
    assert result[0]["program_type"] == "single_celltype"
    assert result[0]["interaction_score"] < 0.3

    # Program 1 should be interaction or microenvironment
    assert result[1]["program_type"] in ["interaction", "microenvironment"]
    assert result[1]["interaction_score"] > 0.3
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_assign_program_cell_types -v`
Expected: FAIL with "cannot import name '_assign_program_cell_types'"

**Step 3: Implement _assign_program_cell_types function**

Add to `CITEgeist/model/anchored_program_discovery.py` in the helper functions section (~line 920):

```python
def _assign_program_cell_types(
    H: NDArray[np.floating],
    cell_type_proportions: pd.DataFrame,
    single_threshold: float = 0.7,
    interaction_threshold: float = 0.25,
) -> List[Dict[str, Any]]:
    """
    Assign cell type labels to joint programs based on correlation with proportions.

    Args:
        H: Program activities (K_programs, n_spots).
        cell_type_proportions: Cell type proportions per spot (n_spots, n_celltypes).
        single_threshold: Min enrichment for single cell-type classification.
        interaction_threshold: Min enrichment for secondary cell type in interaction.

    Returns:
        List of dicts with cell type assignment info per program.
    """
    K_programs = H.shape[0]
    cell_types = list(cell_type_proportions.columns)

    results = []

    for k in range(K_programs):
        h_k = H[k, :]

        # Compute correlation with each cell type's proportions
        enrichments = {}
        for ct in cell_types:
            if ct == "Unknown":
                continue
            props = cell_type_proportions[ct].values
            # Use Spearman correlation (rank-based, more robust)
            if np.std(h_k) > 1e-10 and np.std(props) > 1e-10:
                corr, _ = spearmanr(h_k, props)
                enrichments[ct] = max(0, corr)  # Only positive correlations
            else:
                enrichments[ct] = 0.0

        # Normalize to sum to 1
        total = sum(enrichments.values())
        if total > 0:
            enrichments = {ct: v / total for ct, v in enrichments.items()}

        # Sort by enrichment
        sorted_cts = sorted(enrichments.items(), key=lambda x: x[1], reverse=True)

        primary_ct = sorted_cts[0][0] if sorted_cts else "Unknown"
        primary_score = sorted_cts[0][1] if sorted_cts else 0.0

        secondary_ct = sorted_cts[1][0] if len(sorted_cts) > 1 else None
        secondary_score = sorted_cts[1][1] if len(sorted_cts) > 1 else 0.0

        # Compute interaction score: how balanced is the distribution?
        # 0 = all in one cell type, 1 = perfectly balanced
        if len(enrichments) > 1:
            max_enrich = max(enrichments.values())
            # Gini-like: 1 - (max / ideal_balanced)
            ideal_balanced = 1.0 / len(enrichments)
            interaction_score = 1.0 - (max_enrich - ideal_balanced) / (1.0 - ideal_balanced)
            interaction_score = max(0, min(1, interaction_score))
        else:
            interaction_score = 0.0

        # Classify program type
        if primary_score >= single_threshold:
            program_type = "single_celltype"
        elif secondary_score >= interaction_threshold:
            program_type = "interaction"
        else:
            program_type = "microenvironment"

        results.append({
            "cell_type_enrichments": enrichments,
            "primary_cell_type": primary_ct,
            "secondary_cell_type": secondary_ct if secondary_score >= interaction_threshold else None,
            "interaction_score": interaction_score,
            "program_type": program_type,
        })

    return results
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_assign_program_cell_types -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_joint.py
git commit -m "feat(module4): add cell type assignment for joint programs"
```

---

## Task 4: Implement discover_joint_programs Main Function

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py`
- Modify: `tests/test_module4_joint.py`

**Step 1: Write the test**

Add to `tests/test_module4_joint.py`:

```python
def test_discover_joint_programs_simulated():
    """Test joint program discovery on simulated data."""
    from CITEgeist.model.anchored_program_discovery import (
        discover_joint_programs,
        JointDiscoveryResult,
    )

    # Create simulated AnnData with deconvolved layers
    np.random.seed(42)
    n_spots = 100
    n_genes = 200

    # Base expression
    X = np.random.rand(n_spots, n_genes) * 10

    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Spot_{i}" for i in range(n_spots)]
    adata.obsm["spatial"] = np.random.rand(n_spots, 2) * 1000

    # Add deconvolved layers
    adata.layers["Cancer_genes_pass1"] = np.random.rand(n_spots, n_genes) * 5
    adata.layers["Macrophage_genes_pass1"] = np.random.rand(n_spots, n_genes) * 3
    adata.layers["Tcell_genes_pass1"] = np.random.rand(n_spots, n_genes) * 2

    # Create proportions
    proportions = pd.DataFrame({
        "Cancer": np.random.rand(n_spots),
        "Macrophage": np.random.rand(n_spots),
        "Tcell": np.random.rand(n_spots),
    }, index=adata.obs_names)
    # Normalize rows
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    # Run joint discovery
    result = discover_joint_programs(
        adata=adata,
        cell_type_proportions=proportions,
        K_programs=5,
        layer_pattern="_genes_pass1",
    )

    assert isinstance(result, JointDiscoveryResult)
    assert len(result.programs) == 5
    assert result.W.shape[0] == n_genes
    assert result.W.shape[1] == 5
    assert result.H.shape[0] == 5
    assert result.H.shape[1] == n_spots
    assert all(p.program_type in ["single_celltype", "interaction", "microenvironment"]
               for p in result.programs)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_discover_joint_programs_simulated -v`
Expected: FAIL with "cannot import name 'discover_joint_programs'"

**Step 3: Implement discover_joint_programs function**

Add to `CITEgeist/model/anchored_program_discovery.py` after the existing `discover_anchored_programs` function:

```python
def discover_joint_programs(
    adata: sc.AnnData,
    cell_type_proportions: pd.DataFrame,
    K_programs: int = 10,
    layer_pattern: str = "_genes_pass1",
    lambda_spatial: float = 0.1,
    lambda_sparsity: float = 0.01,
    top_n_genes: int = 50,
    random_state: int = 42,
) -> JointDiscoveryResult:
    """
    Discover spatial programs jointly across all cell types.

    Unlike discover_anchored_programs which runs NMF per cell type, this function
    stacks all deconvolved layers and runs a single NMF to find programs that
    may span multiple cell types (e.g., tumor-immune interface programs).

    Args:
        adata: AnnData with deconvolved layers from Module 3.
        cell_type_proportions: Module 3 output - cell type proportions per spot.
        K_programs: Number of programs to discover.
        layer_pattern: Pattern to identify deconvolved layers.
        lambda_spatial: Spatial smoothness regularization (not yet implemented).
        lambda_sparsity: Sparsity regularization via NMF alpha.
        top_n_genes: Number of top genes to report per program.
        random_state: Random seed.

    Returns:
        JointDiscoveryResult with programs and cell type assignments.
    """
    logger.info("Starting JOINT program discovery across all cell types")

    # Stack deconvolved layers
    stacked_adata = stack_deconvolved_layers(
        adata,
        layer_pattern=layer_pattern,
    )

    # Get stacked expression matrix
    if scipy.sparse.issparse(stacked_adata.X):
        X_stacked = stacked_adata.X.toarray()
    else:
        X_stacked = np.array(stacked_adata.X)

    # Ensure non-negative
    X_stacked = np.maximum(X_stacked, 0)

    gene_names = list(stacked_adata.var_names)
    cell_type_names = list(stacked_adata.obs["cell_type"].unique())
    n_spots = adata.shape[0]
    n_genes = len(gene_names)

    logger.info(f"Stacked data: {X_stacked.shape[0]} rows ({n_spots} spots × {len(cell_type_names)} cell types)")
    logger.info(f"Discovering {K_programs} joint programs")

    # Run NMF on stacked data
    nmf = NMF(
        n_components=K_programs,
        init="nndsvda",
        random_state=random_state,
        max_iter=500,
        alpha_W=lambda_sparsity,
        alpha_H=0.0,
        l1_ratio=0.5,
    )

    # W: (n_stacked_rows, K) -> transpose to (n_genes, K)
    # H: (K, n_stacked_rows)
    W_stacked = nmf.fit_transform(X_stacked)  # (n_stacked, K)
    H_stacked = nmf.components_  # (K, n_genes) - this is transposed from typical convention

    # NMF convention: X ≈ W @ H, so W is (n_samples, K), H is (K, n_features)
    # We want W as gene loadings, H as spot activities
    # So: X.T ≈ W_genes @ H_spots means we need to transpose

    # Actually for stacked: X_stacked is (n_spots*n_ct, n_genes)
    # W_stacked is (n_spots*n_ct, K) - activities per stacked row
    # H_stacked is (K, n_genes) - gene loadings

    # We want:
    # W = gene loadings (n_genes, K)
    # H = spot activities (K, n_spots) - aggregated across cell types

    W = H_stacked.T  # (n_genes, K)

    # Aggregate H across cell types by taking mean per spot
    H_by_celltype = unstack_program_results(W_stacked.T, stacked_adata, n_spots)

    # Average across cell types for overall spot activity
    H = np.zeros((K_programs, n_spots))
    for ct_H in H_by_celltype.values():
        H += ct_H
    H /= len(H_by_celltype)

    # Compute reconstruction error
    X_reconstructed = W_stacked @ H_stacked
    recon_error = np.mean((X_stacked - X_reconstructed) ** 2)

    # Assign cell types to programs
    cell_type_assignments = _assign_program_cell_types(H, cell_type_proportions)

    # Get spatial coordinates
    coords = adata.obsm.get("spatial", np.zeros((n_spots, 2)))

    # Build JointProgram objects
    programs = []
    total_var = np.var(X_stacked)

    for k in range(K_programs):
        # Top genes
        loadings = W[:, k]
        top_indices = np.argsort(loadings)[::-1][:top_n_genes]
        top_genes = [gene_names[i] for i in top_indices]
        gene_loadings = {gene_names[i]: float(loadings[i]) for i in top_indices}

        # Variance explained
        program_var = np.var(H[k, :]) * np.sum(loadings ** 2)
        var_explained = program_var / total_var if total_var > 0 else 0

        # Spatial Moran's I
        moran_i, moran_p = _compute_spatial_moran_i(H[k, :], coords, k=8, n_permutations=99)

        # Activity stats
        mean_activity = float(np.mean(H[k, :]))
        median_activity = float(np.median(H[k, :]))
        active_fraction = float(np.mean(H[k, :] > median_activity))

        # Cell type assignment
        ct_info = cell_type_assignments[k]

        programs.append(JointProgram(
            program_id=k,
            top_genes=top_genes,
            gene_loadings=gene_loadings,
            variance_explained=var_explained,
            spatial_moran_i=moran_i,
            spatial_moran_pvalue=moran_p,
            mean_activity=mean_activity,
            active_spots_fraction=active_fraction,
            cell_type_enrichments=ct_info["cell_type_enrichments"],
            primary_cell_type=ct_info["primary_cell_type"],
            secondary_cell_type=ct_info["secondary_cell_type"],
            interaction_score=ct_info["interaction_score"],
            program_type=ct_info["program_type"],
        ))

    logger.info(f"Discovered {len(programs)} joint programs")
    logger.info(f"  Single cell-type: {sum(1 for p in programs if p.program_type == 'single_celltype')}")
    logger.info(f"  Interaction: {sum(1 for p in programs if p.program_type == 'interaction')}")
    logger.info(f"  Microenvironment: {sum(1 for p in programs if p.program_type == 'microenvironment')}")

    return JointDiscoveryResult(
        programs=programs,
        W=W,
        H=H,
        gene_names=gene_names,
        cell_type_names=cell_type_names,
        n_spots=n_spots,
        reconstruction_error=recon_error,
        H_by_celltype=H_by_celltype,
        parameters={
            "K_programs": K_programs,
            "layer_pattern": layer_pattern,
            "lambda_spatial": lambda_spatial,
            "lambda_sparsity": lambda_sparsity,
            "random_state": random_state,
        },
    )
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_discover_joint_programs_simulated -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_joint.py
git commit -m "feat(module4): implement discover_joint_programs for cross-cell-type programs"
```

---

## Task 5: Export New Functions in __init__.py

**Files:**
- Modify: `CITEgeist/model/__init__.py`

**Step 1: Write the test**

Add to `tests/test_module4_joint.py`:

```python
def test_exports():
    """Test that new classes are exported from CITEgeist.model."""
    from CITEgeist.model import (
        JointProgram,
        JointDiscoveryResult,
        discover_joint_programs,
    )

    assert JointProgram is not None
    assert JointDiscoveryResult is not None
    assert callable(discover_joint_programs)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_exports -v`
Expected: FAIL with "cannot import name 'JointProgram'"

**Step 3: Update __init__.py**

Add to `CITEgeist/model/__init__.py` imports:

```python
from .anchored_program_discovery import (
    # ... existing imports ...
    JointProgram,
    JointDiscoveryResult,
    discover_joint_programs,
)
```

And add to `__all__`:

```python
__all__ = [
    # ... existing exports ...
    "JointProgram",
    "JointDiscoveryResult",
    "discover_joint_programs",
]
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py::test_exports -v`
Expected: PASS

**Step 5: Commit**

```bash
git add CITEgeist/model/__init__.py tests/test_module4_joint.py
git commit -m "feat(module4): export joint program discovery classes"
```

---

## Task 6: Run Full Test Suite and Integration Test

**Files:**
- Modify: `tests/test_module4_joint.py` (add integration test)

**Step 1: Add integration test with real-ish data**

Add to `tests/test_module4_joint.py`:

```python
def test_joint_vs_anchored_comparison():
    """Test that joint and anchored discovery produce comparable results."""
    from CITEgeist.model import (
        discover_joint_programs,
        discover_anchored_programs,
    )

    # Create richer simulated data
    np.random.seed(42)
    n_spots = 200
    n_genes = 500

    X = np.random.rand(n_spots, n_genes) * 10
    adata = sc.AnnData(X=X)
    adata.var_names = [f"Gene_{i}" for i in range(n_genes)]
    adata.obs_names = [f"Spot_{i}" for i in range(n_spots)]

    # Create spatial coords with some structure
    adata.obsm["spatial"] = np.column_stack([
        np.random.rand(n_spots) * 1000,
        np.random.rand(n_spots) * 1000,
    ])

    # Create structured deconvolved layers
    # Cancer high in left side
    cancer_weight = 1 - (adata.obsm["spatial"][:, 0] / 1000)
    adata.layers["Cancer_genes_pass1"] = X * cancer_weight[:, np.newaxis]

    # Macrophage high in right side
    macro_weight = adata.obsm["spatial"][:, 0] / 1000
    adata.layers["Macrophage_genes_pass1"] = X * macro_weight[:, np.newaxis]

    # Proportions matching the spatial pattern
    proportions = pd.DataFrame({
        "Cancer": cancer_weight,
        "Macrophage": macro_weight,
    }, index=adata.obs_names)
    proportions = proportions.div(proportions.sum(axis=1), axis=0)

    # Run joint discovery
    joint_result = discover_joint_programs(
        adata=adata,
        cell_type_proportions=proportions,
        K_programs=4,
    )

    # Verify structure
    assert len(joint_result.programs) == 4
    assert joint_result.H.shape == (4, n_spots)

    # At least one program should be spatial (Moran's I > 0.1)
    spatial_programs = [p for p in joint_result.programs if p.spatial_moran_i > 0.1]
    assert len(spatial_programs) > 0, "Expected at least one spatial program"

    # Programs should have different cell type assignments
    primary_types = [p.primary_cell_type for p in joint_result.programs]
    assert len(set(primary_types)) > 1 or any(
        p.program_type == "interaction" for p in joint_result.programs
    ), "Expected diverse cell type assignments"
```

**Step 2: Run full test suite**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_joint.py -v`
Expected: All tests PASS

**Step 3: Commit**

```bash
git add tests/test_module4_joint.py
git commit -m "test(module4): add integration test for joint program discovery"
```

---

## Summary

After completing all tasks, you will have:

1. **JointProgram** dataclass - represents a single joint program with cell type assignments
2. **JointDiscoveryResult** dataclass - holds all discovery results
3. **_assign_program_cell_types()** - helper to classify programs by cell type
4. **discover_joint_programs()** - main entry point for joint discovery
5. **Exports** in `__init__.py`
6. **Comprehensive tests** in `tests/test_module4_joint.py`

The implementation follows the existing patterns in the codebase and integrates with:
- `stack_deconvolved_layers()` for data preparation
- `unstack_program_results()` for splitting results
- `_compute_spatial_moran_i()` for spatial statistics
- Module 5 integration (uses same W/H format as anchored discovery)
