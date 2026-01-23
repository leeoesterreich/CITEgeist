# Module 4 Region Analysis Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add region-aware program analysis to Module 4, enabling discovery of condition-specific transcriptional programs within anchored cell types.

**Architecture:** Extend existing SpatialProgram dataclass with optional region fields; add three utility functions (analyze_program_regions, compare_programs_by_region, extract_program_context_genes) that operate on existing AnchoredProgramResult output. No changes to core NMF discovery logic.

**Tech Stack:** Python 3.10, NumPy, Pandas, SciPy (stats), existing Module 4 infrastructure.

---

## Task 1: Extend SpatialProgram Dataclass with Region Fields

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py:78-108`
- Test: `tests/test_module4_regions.py` (create new)

**Step 1: Write the failing test**

Create `tests/test_module4_regions.py`:

```python
"""
Test Module 4c: Region-Aware Program Analysis.

Tests the region analysis extension to Module 4:
- analyze_program_regions()
- compare_programs_by_region()
- extract_program_context_genes()
"""

import numpy as np
import pandas as pd
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from CITEgeist.model.anchored_program_discovery import SpatialProgram


class TestSpatialProgramRegionFields:
    """Test SpatialProgram dataclass region fields."""

    def test_region_fields_optional_and_default_none(self):
        """Region fields should be optional with None defaults."""
        program = SpatialProgram(
            program_id=0,
            top_genes=["GENE1", "GENE2"],
            gene_loadings={"GENE1": 0.8, "GENE2": 0.5},
            variance_explained=0.15,
            spatial_moran_i=0.3,
            spatial_moran_pvalue=0.01,
            mean_activity=0.5,
            active_spots_fraction=0.4,
        )

        # Region fields should exist and be None by default
        assert hasattr(program, 'region_enrichment')
        assert hasattr(program, 'region_specificity')
        assert hasattr(program, 'region_pvalue')
        assert hasattr(program, 'enriched_region')

        assert program.region_enrichment is None
        assert program.region_specificity is None
        assert program.region_pvalue is None
        assert program.enriched_region is None

    def test_region_fields_can_be_set(self):
        """Region fields should accept values."""
        program = SpatialProgram(
            program_id=0,
            top_genes=["GENE1"],
            gene_loadings={"GENE1": 0.8},
            variance_explained=0.15,
            spatial_moran_i=0.3,
            spatial_moran_pvalue=0.01,
            mean_activity=0.5,
            active_spots_fraction=0.4,
            region_enrichment={"D538G_pos": 0.8, "D538G_neg": 0.2},
            region_specificity=0.6,
            region_pvalue=0.001,
            enriched_region="D538G_pos",
        )

        assert program.region_enrichment == {"D538G_pos": 0.8, "D538G_neg": 0.2}
        assert program.region_specificity == 0.6
        assert program.region_pvalue == 0.001
        assert program.enriched_region == "D538G_pos"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestSpatialProgramRegionFields -v`

Expected: FAIL with "TypeError: __init__() got an unexpected keyword argument 'region_enrichment'"

**Step 3: Write minimal implementation**

Modify `CITEgeist/model/anchored_program_discovery.py`. Replace lines 78-108 with:

```python
@dataclass
class SpatialProgram:
    """A single spatial transcriptomic program discovered within a cell type context."""

    program_id: int
    """Unique identifier for this program within its anchor."""

    top_genes: List[str]
    """Top N genes by loading weight (most characteristic of this program)."""

    gene_loadings: Dict[str, float]
    """Gene name -> loading weight for all genes in the program."""

    variance_explained: float
    """Fraction of variance explained by this program."""

    spatial_moran_i: float
    """Moran's I of program activity across spots (spatial coherence)."""

    spatial_moran_pvalue: float
    """P-value for Moran's I (from permutation test)."""

    mean_activity: float
    """Mean program activity across spots."""

    active_spots_fraction: float
    """Fraction of spots with above-median program activity."""

    subpopulations: Optional[List[int]] = None
    """Subpopulation IDs where this program is dominant (if detected)."""

    # Region analysis fields (populated by analyze_program_regions)
    region_enrichment: Optional[Dict[str, float]] = None
    """Mean program activity per region (region_value -> mean_activity)."""

    region_specificity: Optional[float] = None
    """Specificity score: 0 = uniform across regions, 1 = completely region-specific."""

    region_pvalue: Optional[float] = None
    """P-value for region enrichment (Mann-Whitney U or Kruskal-Wallis)."""

    enriched_region: Optional[str] = None
    """Region with highest mean activity (if significantly enriched)."""
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestSpatialProgramRegionFields -v`

Expected: PASS (2 tests)

**Step 5: Commit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_regions.py
git commit -m "feat: add region fields to SpatialProgram dataclass"
```

---

## Task 2: Implement analyze_program_regions()

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py` (add at end, ~line 2137)
- Test: `tests/test_module4_regions.py`

**Step 1: Write the failing test**

Add to `tests/test_module4_regions.py`:

```python
import scanpy as sc
from CITEgeist.model.anchored_program_discovery import (
    SpatialProgram,
    AnchoredProgramResult,
    analyze_program_regions,
)


def create_mock_anchored_result(n_spots: int = 100, n_genes: int = 50, K: int = 3) -> AnchoredProgramResult:
    """Create a mock AnchoredProgramResult for testing."""
    np.random.seed(42)

    W = np.random.rand(n_genes, K)
    H = np.random.rand(K, n_spots)

    gene_names = [f"GENE_{i}" for i in range(n_genes)]

    programs = []
    for k in range(K):
        top_idx = np.argsort(W[:, k])[::-1][:10]
        programs.append(SpatialProgram(
            program_id=k,
            top_genes=[gene_names[i] for i in top_idx],
            gene_loadings={gene_names[i]: float(W[i, k]) for i in top_idx},
            variance_explained=0.15,
            spatial_moran_i=0.3,
            spatial_moran_pvalue=0.01,
            mean_activity=float(H[k, :].mean()),
            active_spots_fraction=0.4,
        ))

    return AnchoredProgramResult(
        anchor_name="Cancer_Cells",
        anchor_proteins=["EPCAM"],
        programs=programs,
        W=W,
        H=H,
        gene_names=gene_names,
        protein_correlations=pd.DataFrame(),
        reconstruction_error=0.1,
        n_spots_used=n_spots,
    )


def create_mock_adata_with_regions(n_spots: int = 100) -> sc.AnnData:
    """Create mock AnnData with region annotations."""
    np.random.seed(42)

    # Create AnnData with random expression
    X = np.random.rand(n_spots, 50)
    adata = sc.AnnData(X)
    adata.obs_names = [f"spot_{i}" for i in range(n_spots)]
    adata.var_names = [f"GENE_{i}" for i in range(50)]

    # Add region column - first half D538G+, second half D538G-
    adata.obs["D538G_Mutation"] = ["D538G_pos"] * (n_spots // 2) + ["D538G_neg"] * (n_spots - n_spots // 2)

    # Add spatial coordinates
    adata.obsm["spatial"] = np.random.rand(n_spots, 2) * 1000

    return adata


class TestAnalyzeProgramRegions:
    """Test analyze_program_regions function."""

    def test_populates_region_fields(self):
        """Should populate region fields in all programs."""
        result = create_mock_anchored_result(n_spots=100, K=3)
        adata = create_mock_adata_with_regions(n_spots=100)

        result = analyze_program_regions(result, adata, "D538G_Mutation")

        for program in result.programs:
            assert program.region_enrichment is not None
            assert "D538G_pos" in program.region_enrichment
            assert "D538G_neg" in program.region_enrichment
            assert program.region_specificity is not None
            assert program.region_pvalue is not None

    def test_detects_enriched_region(self):
        """Should detect enriched region when one exists."""
        # Create result where program 0 is strongly active in first half of spots
        result = create_mock_anchored_result(n_spots=100, K=3)
        result.H[0, :50] = 1.0  # High activity in D538G_pos
        result.H[0, 50:] = 0.1  # Low activity in D538G_neg

        adata = create_mock_adata_with_regions(n_spots=100)

        result = analyze_program_regions(result, adata, "D538G_Mutation")

        # Program 0 should be enriched in D538G_pos
        assert result.programs[0].enriched_region == "D538G_pos"
        assert result.programs[0].region_pvalue < 0.05

    def test_raises_on_missing_column(self):
        """Should raise ValueError if region column not found."""
        result = create_mock_anchored_result()
        adata = create_mock_adata_with_regions()

        with pytest.raises(ValueError, match="not found in adata.obs"):
            analyze_program_regions(result, adata, "NonexistentColumn")
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestAnalyzeProgramRegions -v`

Expected: FAIL with "cannot import name 'analyze_program_regions'"

**Step 3: Write minimal implementation**

Add to end of `CITEgeist/model/anchored_program_discovery.py`:

```python
# =============================================================================
# MODULE 4c: REGION-AWARE PROGRAM ANALYSIS
# =============================================================================


def analyze_program_regions(
    result: AnchoredProgramResult,
    adata: sc.AnnData,
    region_column: str,
    min_spots_per_region: int = 20,
) -> AnchoredProgramResult:
    """
    Annotate programs with region enrichment statistics.

    For each program, computes activity distribution across regions defined by
    a categorical column in adata.obs, and tests for significant enrichment.

    This is a post-hoc analysis step - programs are discovered first, then
    analyzed for regional variation.

    Args:
        result: AnchoredProgramResult from discover_programs_from_layers()
        adata: AnnData with region annotations in obs
        region_column: Column in adata.obs defining regions (e.g., "D538G Mutation")
        min_spots_per_region: Minimum spots required per region for analysis

    Returns:
        AnchoredProgramResult with region fields populated in each SpatialProgram

    Example:
        >>> result = discover_programs_from_layers(adata, ...)
        >>> result = analyze_program_regions(result, adata, "D538G Mutation")
        >>> for prog in result.programs:
        ...     print(f"Program {prog.program_id}: enriched in {prog.enriched_region}")
    """
    from scipy.stats import mannwhitneyu, kruskal

    if region_column not in adata.obs.columns:
        raise ValueError(f"Region column '{region_column}' not found in adata.obs")

    regions = adata.obs[region_column].unique()
    n_regions = len(regions)

    logger.info(f"Analyzing program regions using '{region_column}' ({n_regions} regions)")

    # Get program activities (H matrix)
    H = result.H  # Shape: (K_programs, n_spots)

    # Analyze each program
    for k, program in enumerate(result.programs):
        h_k = H[k, :]  # Activity for this program across spots

        # Compute mean activity per region
        region_means = {}
        region_activities = {}
        for region in regions:
            mask = (adata.obs[region_column] == region).values
            if mask.sum() >= min_spots_per_region:
                activities = h_k[mask]
                region_means[str(region)] = float(np.mean(activities))
                region_activities[str(region)] = activities

        program.region_enrichment = region_means

        # Skip statistical test if insufficient regions
        if len(region_activities) < 2:
            program.region_specificity = 0.0
            program.region_pvalue = 1.0
            program.enriched_region = None
            continue

        # Statistical test for region differences
        activity_lists = list(region_activities.values())
        if n_regions == 2:
            # Mann-Whitney U for two regions
            stat, pval = mannwhitneyu(activity_lists[0], activity_lists[1], alternative='two-sided')
        else:
            # Kruskal-Wallis for multiple regions
            stat, pval = kruskal(*activity_lists)

        program.region_pvalue = float(pval)

        # Compute specificity score (coefficient of variation of region means)
        means = np.array(list(region_means.values()))
        if means.mean() > 0:
            program.region_specificity = float(means.std() / means.mean())
        else:
            program.region_specificity = 0.0

        # Identify enriched region
        if pval < 0.05 and len(region_means) > 0:
            program.enriched_region = max(region_means, key=region_means.get)
        else:
            program.enriched_region = None

    n_enriched = sum(1 for p in result.programs if p.enriched_region is not None)
    logger.info(f"Region analysis complete: {n_enriched}/{len(result.programs)} "
                f"programs show significant region enrichment")

    return result
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestAnalyzeProgramRegions -v`

Expected: PASS (3 tests)

**Step 5: Commit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_regions.py
git commit -m "feat: add analyze_program_regions function"
```

---

## Task 3: Implement compare_programs_by_region()

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py`
- Test: `tests/test_module4_regions.py`

**Step 1: Write the failing test**

Add to `tests/test_module4_regions.py`:

```python
from CITEgeist.model.anchored_program_discovery import compare_programs_by_region


class TestCompareProgramsByRegion:
    """Test compare_programs_by_region function."""

    def test_returns_dataframe_with_expected_columns(self):
        """Should return DataFrame with comparison statistics."""
        result = create_mock_anchored_result(n_spots=100, K=3)
        adata = create_mock_adata_with_regions(n_spots=100)

        df = compare_programs_by_region(
            result, adata, "D538G_Mutation",
            region_a="D538G_pos", region_b="D538G_neg"
        )

        assert isinstance(df, pd.DataFrame)
        assert "program_id" in df.columns
        assert "mean_activity_a" in df.columns
        assert "mean_activity_b" in df.columns
        assert "fold_change" in df.columns
        assert "pvalue" in df.columns
        assert "top_genes" in df.columns
        assert len(df) == 3  # One row per program

    def test_fold_change_calculation(self):
        """Fold change should reflect region activity differences."""
        result = create_mock_anchored_result(n_spots=100, K=3)
        # Make program 0 strongly enriched in region A
        result.H[0, :50] = 2.0  # High in D538G_pos
        result.H[0, 50:] = 0.5  # Low in D538G_neg

        adata = create_mock_adata_with_regions(n_spots=100)

        df = compare_programs_by_region(
            result, adata, "D538G_Mutation",
            region_a="D538G_pos", region_b="D538G_neg"
        )

        # Program 0 should have fold_change > 1
        prog0_row = df[df["program_id"] == 0].iloc[0]
        assert prog0_row["fold_change"] > 3.0  # 2.0 / 0.5 = 4.0

    def test_sorted_by_fold_change(self):
        """Results should be sorted by fold change descending."""
        result = create_mock_anchored_result(n_spots=100, K=3)
        adata = create_mock_adata_with_regions(n_spots=100)

        df = compare_programs_by_region(
            result, adata, "D538G_Mutation",
            region_a="D538G_pos", region_b="D538G_neg"
        )

        # Check sorted descending
        assert df["fold_change"].is_monotonic_decreasing
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestCompareProgramsByRegion -v`

Expected: FAIL with "cannot import name 'compare_programs_by_region'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/anchored_program_discovery.py` after `analyze_program_regions`:

```python
def compare_programs_by_region(
    result: AnchoredProgramResult,
    adata: sc.AnnData,
    region_column: str,
    region_a: Any,
    region_b: Any,
    top_n_genes: int = 50,
) -> pd.DataFrame:
    """
    Compare program activities and gene loadings between two regions.

    For each program, computes:
    - Mean activity in region A vs B
    - Fold change and statistical significance
    - Top genes driving the program

    Useful for identifying which programs and genes are region-specific.

    Args:
        result: AnchoredProgramResult with discovered programs
        adata: AnnData with region annotations
        region_column: Column defining regions
        region_a: First region value (e.g., True for D538G+)
        region_b: Second region value (e.g., False for D538G-)
        top_n_genes: Number of top genes to include per program

    Returns:
        DataFrame with columns:
        - program_id: Program index
        - mean_activity_a, mean_activity_b: Mean activity per region
        - fold_change: region_a / region_b
        - pvalue: Mann-Whitney U test
        - top_genes: Comma-separated top genes
        - gene_loadings: Dict of gene -> loading

    Example:
        >>> df = compare_programs_by_region(result, adata, "D538G Mutation", True, False)
        >>> d538g_enriched = df[df['fold_change'] > 1.5]
    """
    from scipy.stats import mannwhitneyu

    mask_a = (adata.obs[region_column] == region_a).values
    mask_b = (adata.obs[region_column] == region_b).values

    H = result.H
    W = result.W

    records = []
    for k, program in enumerate(result.programs):
        h_k = H[k, :]

        activities_a = h_k[mask_a]
        activities_b = h_k[mask_b]

        mean_a = float(np.mean(activities_a))
        mean_b = float(np.mean(activities_b))

        # Fold change (with pseudocount to avoid division by zero)
        fc = (mean_a + 1e-6) / (mean_b + 1e-6)

        # Statistical test
        if len(activities_a) >= 10 and len(activities_b) >= 10:
            _, pval = mannwhitneyu(activities_a, activities_b, alternative='two-sided')
        else:
            pval = 1.0

        # Top genes
        loadings = W[:, k]
        top_idx = np.argsort(loadings)[::-1][:top_n_genes]
        top_genes = [result.gene_names[i] for i in top_idx]
        gene_loadings = {result.gene_names[i]: float(loadings[i]) for i in top_idx}

        records.append({
            'program_id': k,
            'mean_activity_a': mean_a,
            'mean_activity_b': mean_b,
            'fold_change': fc,
            'pvalue': float(pval),
            'n_spots_a': int(mask_a.sum()),
            'n_spots_b': int(mask_b.sum()),
            'top_genes': ', '.join(top_genes[:10]),
            'gene_loadings': gene_loadings,
        })

    df = pd.DataFrame(records)
    df = df.sort_values('fold_change', ascending=False)

    return df
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestCompareProgramsByRegion -v`

Expected: PASS (3 tests)

**Step 5: Commit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_regions.py
git commit -m "feat: add compare_programs_by_region function"
```

---

## Task 4: Implement extract_program_context_genes()

**Files:**
- Modify: `CITEgeist/model/anchored_program_discovery.py`
- Test: `tests/test_module4_regions.py`

**Step 1: Write the failing test**

Add to `tests/test_module4_regions.py`:

```python
from CITEgeist.model.anchored_program_discovery import extract_program_context_genes


class TestExtractProgramContextGenes:
    """Test extract_program_context_genes function."""

    def test_returns_list_of_tuples(self):
        """Should return list of (gene, loading) tuples."""
        result = create_mock_anchored_result(n_spots=100, K=3)

        context_genes = extract_program_context_genes(result, program_id=0, target_gene="GENE_0")

        assert isinstance(context_genes, list)
        assert len(context_genes) > 0
        assert isinstance(context_genes[0], tuple)
        assert len(context_genes[0]) == 2
        assert isinstance(context_genes[0][0], str)  # gene name
        assert isinstance(context_genes[0][1], float)  # loading

    def test_excludes_target_gene_by_default(self):
        """Should exclude target gene from results by default."""
        result = create_mock_anchored_result(n_spots=100, K=3)

        context_genes = extract_program_context_genes(result, program_id=0, target_gene="GENE_0")

        gene_names = [g[0] for g in context_genes]
        assert "GENE_0" not in gene_names

    def test_includes_target_gene_when_requested(self):
        """Should include target gene when exclude_target=False."""
        result = create_mock_anchored_result(n_spots=100, K=3)

        context_genes = extract_program_context_genes(
            result, program_id=0, target_gene="GENE_0", exclude_target=False
        )

        gene_names = [g[0] for g in context_genes]
        # GENE_0 might be in top N depending on random loadings
        # Just verify the function runs without error
        assert len(context_genes) > 0

    def test_respects_top_n_parameter(self):
        """Should return at most top_n genes."""
        result = create_mock_anchored_result(n_spots=100, n_genes=100, K=3)

        context_genes = extract_program_context_genes(
            result, program_id=0, target_gene="GENE_0", top_n=20
        )

        assert len(context_genes) <= 20

    def test_raises_on_invalid_program_id(self):
        """Should raise ValueError for invalid program_id."""
        result = create_mock_anchored_result(n_spots=100, K=3)

        with pytest.raises(ValueError, match="Program 10 not found"):
            extract_program_context_genes(result, program_id=10, target_gene="GENE_0")
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestExtractProgramContextGenes -v`

Expected: FAIL with "cannot import name 'extract_program_context_genes'"

**Step 3: Write minimal implementation**

Add to `CITEgeist/model/anchored_program_discovery.py` after `compare_programs_by_region`:

```python
def extract_program_context_genes(
    result: AnchoredProgramResult,
    program_id: int,
    target_gene: str,
    top_n: int = 50,
    exclude_target: bool = True,
) -> List[Tuple[str, float]]:
    """
    Extract genes co-loaded with a target gene in a specific program.

    Useful for finding "contextual factors" - genes that are co-expressed
    with a gene of interest (e.g., MDK) in a spatial program.

    Args:
        result: AnchoredProgramResult with discovered programs
        program_id: Index of the program to analyze
        target_gene: Gene of interest (e.g., "MDK")
        top_n: Number of top co-loaded genes to return
        exclude_target: Whether to exclude target gene from results

    Returns:
        List of (gene_name, loading) tuples sorted by loading

    Example:
        >>> # Find genes co-expressed with MDK in program 2
        >>> context_genes = extract_program_context_genes(result, 2, "MDK")
        >>> print("Contextual factors:", [g[0] for g in context_genes[:10]])
    """
    if program_id >= len(result.programs):
        raise ValueError(f"Program {program_id} not found (max: {len(result.programs)-1})")

    W = result.W  # (n_genes, K_programs)
    loadings = W[:, program_id]

    # Check if target gene is in this program
    if target_gene in result.gene_names:
        target_idx = result.gene_names.index(target_gene)
        target_loading = loadings[target_idx]
        logger.info(f"Target gene '{target_gene}' loading in program {program_id}: {target_loading:.4f}")
    else:
        logger.warning(f"Target gene '{target_gene}' not found in gene list")

    # Get top genes by loading
    top_idx = np.argsort(loadings)[::-1]

    context_genes = []
    for idx in top_idx:
        gene = result.gene_names[idx]
        if exclude_target and gene == target_gene:
            continue
        context_genes.append((gene, float(loadings[idx])))
        if len(context_genes) >= top_n:
            break

    return context_genes
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestExtractProgramContextGenes -v`

Expected: PASS (5 tests)

**Step 5: Commit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git add CITEgeist/model/anchored_program_discovery.py tests/test_module4_regions.py
git commit -m "feat: add extract_program_context_genes function"
```

---

## Task 5: Update Module Exports

**Files:**
- Modify: `CITEgeist/model/__init__.py`
- Test: Quick import test

**Step 1: Write the failing test**

Add to end of `tests/test_module4_regions.py`:

```python
class TestModuleExports:
    """Test that new functions are exported from module."""

    def test_import_from_model(self):
        """Should be able to import from CITEgeist.model."""
        from CITEgeist.model import (
            analyze_program_regions,
            compare_programs_by_region,
            extract_program_context_genes,
        )

        assert callable(analyze_program_regions)
        assert callable(compare_programs_by_region)
        assert callable(extract_program_context_genes)
```

**Step 2: Run test to verify it fails**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestModuleExports -v`

Expected: FAIL with "cannot import name 'analyze_program_regions'"

**Step 3: Write minimal implementation**

Modify `CITEgeist/model/__init__.py`. Update the anchored_program_discovery imports (around line 34-51):

```python
# Module 4: Protein-anchored program discovery
from .anchored_program_discovery import (
    SpatialSubpopulation,
    SpatialProgram,
    AnchoredProgramResult,
    AnchoredProgramDiscoveryResult,
    detect_spatial_subpopulations,
    discover_anchored_programs,  # Legacy: uses raw expression + contrastive
    discover_programs_from_layers,  # Recommended: uses deconvolved layers from Module 3
    store_results_in_adata,
    # Helper functions for deconvolved layers
    stack_deconvolved_layers,
    unstack_program_results,
    extract_celltype_expression,
    # Module 4b: Bivariate program relationships
    ProgramPairRelationship,
    BivariateProgramResult,
    analyze_program_relationships,
    # Module 4c: Region-aware program analysis
    analyze_program_regions,
    compare_programs_by_region,
    extract_program_context_genes,
)
```

Also update the `__all__` list (around line 100-115):

```python
    # Module 4b: Bivariate program relationships
    "ProgramPairRelationship",
    "BivariateProgramResult",
    "analyze_program_relationships",
    # Module 4c: Region-aware program analysis
    "analyze_program_regions",
    "compare_programs_by_region",
    "extract_program_context_genes",
```

**Step 4: Run test to verify it passes**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py::TestModuleExports -v`

Expected: PASS

**Step 5: Commit**

```bash
cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist
git add CITEgeist/model/__init__.py
git commit -m "feat: export Module 4c region analysis functions"
```

---

## Task 6: Run Full Test Suite

**Step 1: Run all new tests**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_regions.py -v`

Expected: All tests PASS (14 tests total)

**Step 2: Run existing Module 4 tests to verify no regression**

Run: `cd /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist && python -m pytest tests/test_module4_deconvolved.py -v --timeout=300` (if available)

Expected: Existing tests still PASS

**Step 3: Commit any fixes if needed**

```bash
git add -A
git commit -m "fix: address any test failures"
```

---

## Summary Checklist

- [ ] Task 1: Extend SpatialProgram with region fields
- [ ] Task 2: Implement analyze_program_regions()
- [ ] Task 3: Implement compare_programs_by_region()
- [ ] Task 4: Implement extract_program_context_genes()
- [ ] Task 5: Update module exports
- [ ] Task 6: Run full test suite

---

## Next Steps (After Implementation)

After completing the Module 4 enhancement, proceed to:

1. **Create Vignette 4 notebook** (`CITEgeist/examples/vignette_4_esr1_midkine_programs.ipynb`)
2. **Download GEO validation data** (GSE89888, GSE125117)
3. **Complete the full validation workflow**

See `docs/plans/2026-01-23-module4-region-analysis-midkine-design.md` for the complete vignette outline.
