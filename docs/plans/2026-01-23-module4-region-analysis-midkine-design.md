# Module 4 Region Analysis & Midkine Validation Design

**Goal**: Enhance Module 4 to analyze regional variation within anchored programs, enabling discovery of contextual factors driving phenotype-specific gene expression. Demonstrate utility through ESR1-D538G/midkine secretion story with orthogonal validation.

**Research Question**: What contextual factors in ESR1-mutant cancer cells drive MDK (midkine) secretion, and can these be validated in cell line models?

**Key Finding to Explain**: MCF7 cells show increased MDK secretion with D538G mutation, but T47D cells do not. Module 4 should help identify the "permissive context" genes that explain this difference.

---

## Part 1: Module 4 Enhancement

### Design Philosophy

- Programs discovered by Module 4 naturally capture spatial variation
- Region analysis is a **lens** applied after discovery, not a separate mode
- Keep core Module 4 unchanged; add lightweight utility functions
- Easily extensible to any categorical region definition

### 1.1 Extend SpatialProgram Dataclass

**File**: `CITEgeist/model/anchored_program_discovery.py`

Add optional fields to existing `SpatialProgram` dataclass (after line ~107):

```python
@dataclass
class SpatialProgram:
    """A single spatial transcriptomic program discovered within a cell type context."""

    program_id: int
    top_genes: List[str]
    gene_loadings: Dict[str, float]
    variance_explained: float
    spatial_moran_i: float
    spatial_moran_pvalue: float
    mean_activity: float
    active_spots_fraction: float
    subpopulations: Optional[List[int]] = None

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

### 1.2 Add analyze_program_regions() Function

**File**: `CITEgeist/model/anchored_program_discovery.py`

```python
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
    K = H.shape[0]

    # Analyze each program
    for k, program in enumerate(result.programs):
        h_k = H[k, :]  # Activity for this program across spots

        # Compute mean activity per region
        region_means = {}
        region_activities = {}
        for region in regions:
            mask = adata.obs[region_column] == region
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
        if pval < 0.05:
            program.enriched_region = max(region_means, key=region_means.get)
        else:
            program.enriched_region = None

    logger.info(f"Region analysis complete: {sum(1 for p in result.programs if p.enriched_region)} "
                f"programs show significant region enrichment")

    return result
```

### 1.3 Add compare_programs_by_region() Function

**File**: `CITEgeist/model/anchored_program_discovery.py`

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

    mask_a = adata.obs[region_column] == region_a
    mask_b = adata.obs[region_column] == region_b

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
            'pvalue': pval,
            'n_spots_a': int(mask_a.sum()),
            'n_spots_b': int(mask_b.sum()),
            'top_genes': ', '.join(top_genes[:10]),
            'gene_loadings': gene_loadings,
        })

    df = pd.DataFrame(records)
    df = df.sort_values('fold_change', ascending=False)

    return df
```

### 1.4 Add extract_program_context_genes() Function

**File**: `CITEgeist/model/anchored_program_discovery.py`

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

### 1.5 Update Module Exports

**File**: `CITEgeist/model/__init__.py`

Add new functions to exports:

```python
from .anchored_program_discovery import (
    # ... existing exports ...
    analyze_program_regions,
    compare_programs_by_region,
    extract_program_context_genes,
)
```

---

## Part 2: Vignette 4 - ESR1 Midkine Programs

**File**: `CITEgeist/examples/vignette_4_esr1_midkine_programs.ipynb`

### Vignette Outline

| Part | Title | Cells | Description |
|------|-------|-------|-------------|
| 1 | Setup & Data Loading | 3 | Load HCC22-088-P4-S2, D538G annotations |
| 2 | Module 4 Program Discovery | 3 | Discover programs on EPCAM+ cancer cells |
| 3 | Region Enrichment Analysis | 4 | Identify D538G-enriched programs |
| 4 | MDK Secretion Context | 4 | Extract contextual genes co-loaded with MDK |
| 5 | Load Bulk RNA-seq (GSE89888) | 4 | MCF7 & T47D, D538G vs WT samples |
| 6 | RNA-seq Interaction Validation | 5 | Test contextual genes for MCF7-specific effect |
| 7 | Load ChIP-seq (GSE125117) | 3 | ER binding data, D538G vs WT |
| 8 | ChIP-seq Mechanistic Validation | 4 | Differential ER binding at contextual genes |
| 9 | Integrated Model & Conclusions | 3 | High-confidence permissive factors, figures |

### Validation Data Sources

**GSE89888 (Bulk RNA-seq)**:
- Cell lines: MCF7 (responder), T47D (non-responder)
- ESR1 variants: WT, D538G (ignore Y537S for this analysis)
- Treatment: Vehicle, E2 (1nM, 24h)
- Replicates: 4 per condition
- Use case: Interaction effect test (ESR1_status × cell_line)

**GSE125117 (ChIP-seq)**:
- Same cell lines and variants
- ER binding ± E2 stimulation
- Use case: Validate contextual genes as direct ER targets

### Validation Logic

A high-confidence "permissive context" gene must satisfy:

1. **Spatial discovery**: Co-loaded with MDK in D538G-enriched program
2. **RNA-seq interaction**: Significant ESR1 × cell_line interaction effect
3. **Direction**: Upregulated in MCF7-D538G, unchanged in T47D-D538G
4. **Mechanistic support** (optional): Differential ER binding in D538G vs WT

---

## Part 3: Implementation Plan

### Task 1: Add SpatialProgram Region Fields

**Files**: `CITEgeist/model/anchored_program_discovery.py`

Add optional region fields to SpatialProgram dataclass.

### Task 2: Implement analyze_program_regions()

**Files**: `CITEgeist/model/anchored_program_discovery.py`

Implement function with Mann-Whitney/Kruskal-Wallis tests.

### Task 3: Implement compare_programs_by_region()

**Files**: `CITEgeist/model/anchored_program_discovery.py`

Implement comparative analysis function.

### Task 4: Implement extract_program_context_genes()

**Files**: `CITEgeist/model/anchored_program_discovery.py`

Implement context gene extraction function.

### Task 5: Update Module Exports

**Files**: `CITEgeist/model/__init__.py`

Add new functions to module exports.

### Task 6: Create Vignette 4 Notebook

**Files**: `CITEgeist/examples/vignette_4_esr1_midkine_programs.ipynb`

Create complete vignette with Parts 1-9.

### Task 7: Test with Spatial Data

Run vignette Parts 1-4 on HCC22-088-P4-S2 sample.

### Task 8: Integrate GEO Validation Data

Download and process GSE89888 (RNA-seq) and GSE125117 (ChIP-seq).

### Task 9: Complete Validation Analysis

Run vignette Parts 5-9 with full validation pipeline.

---

## Part 4: Expected Outputs

### From Spatial Analysis (Parts 1-4)

1. **Program summary table**: All programs with region enrichment statistics
2. **D538G-enriched programs**: Programs significantly active in D538G+ regions
3. **MDK secretion program**: Program with highest MDK loading among D538G-enriched
4. **Contextual gene list**: Top 50 genes co-loaded with MDK (validation candidates)

### From RNA-seq Validation (Parts 5-6)

1. **Interaction effect results**: For each contextual gene, p-value and effect size
2. **Validated contextual genes**: Genes with significant interaction AND correct direction
3. **Enrichment statistics**: Are contextual genes enriched for interaction effects vs background?

### From ChIP-seq Validation (Parts 7-8)

1. **ER binding at contextual genes**: Peak presence/absence at promoters
2. **Differential binding**: Genes with altered ER binding in D538G vs WT
3. **Mechanistic candidates**: Contextual genes that are direct ER targets

### Integrated Output (Part 9)

1. **High-confidence permissive factors**: Genes passing all validation criteria
2. **Summary figure**: Spatial program → RNA-seq validation → ChIP-seq mechanism
3. **Gene list for paper**: Exportable table with all statistics

---

## Part 5: Key Biological Interpretation

The vignette will address:

1. **Why does D538G mutation increase MDK secretion in some contexts?**
   - Contextual genes co-expressed with MDK in spatial data
   - Validated as MCF7-specific responders in bulk RNA-seq

2. **What makes MCF7 permissive but not T47D?**
   - Genes with interaction effect (ESR1 × cell_line)
   - Baseline expression differences in permissive factors

3. **Is this direct ER regulation?**
   - ChIP-seq evidence for ER binding at contextual genes
   - Differential binding with D538G mutation

4. **Autocrine signaling component?**
   - Check SDC4/NCL (MDK receptors) in same or correlated programs
   - Evidence for autocrine feedback loop in cancer cells

---

**Author**: Claude (Brainstorming Session)
**Date**: 2026-01-23
**Status**: Design Complete - Ready for Implementation
