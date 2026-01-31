# Module 1-2 Discovery Comparison: Spatial Co-expression vs Leiden Clustering

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Demonstrate that spatially-informed co-expression discovery (Module 1-2) finds more biologically coherent and resolution-robust protein marker programs than standard unsupervised clustering (Leiden). This comparison is for the final paper.

**Framing:** Module 1-2 is not a complete cell type profile generator. It discovers which protein co-expression patterns are spatially interesting and consistent. Users take this information to guide their Module 3 deconvolution profiles — e.g., knowing which macrophage subtypes actually exist in the tissue so they don't waste capacity on absent cell states.

---

## Experiment Structure

**Data:** 5 Xenium tissue regions, each at two resolutions:
- Single-cell (~85K-107K cells per region, 27 protein markers)
- Pseudo-Visium spots (~1,400 spots per region, same 27 markers)

**Comparison arms:**
- **Module 1-2:** Spatial colocalization analysis → co-expression module discovery (current pipeline)
- **Leiden clustering:** Standard Scanpy workflow — normalize, PCA, neighbors graph, Leiden at multiple resolutions (0.3, 0.5, 0.8, 1.0, 1.5), extract top enriched markers per cluster via `rank_genes_groups`

Multiple Leiden resolutions are shown because there is no principled way to choose one. Coarse merges distinct types, fine fragments them, and no single resolution matches Module 1-2's output.

Module 1-2 does not require a resolution parameter. It discovers co-expression modules from pairwise spatial statistics. This is a structural advantage, not a tuning advantage.

---

## Evaluation Metrics

### 1. Biological Coherence

Do discovered marker groups contain markers from a single lineage or merge across lineages?

- Score each module by checking whether all markers map to the same cell lineage (using the curated marker-to-lineage mapping below)
- Report fraction of modules that are "pure" (single-lineage) vs "mixed" (cross-lineage)

### 2. Rare Subtype Detection

Are biologically meaningful rare populations found?

Pre-defined subtypes detectable with this panel:
- Mast cells (KIT + CPA3)
- Smooth muscle (ACTA2 + DES + MYH11)
- Exhausted T cells (PD-1 + LAG-3 + HAVCR2)
- M2 macrophages (CD163 + MS4A4A)
- Plasma cells (MZB1 + SLAMF7)

Binary per method: did it find a module matching each subtype? At what resolution?

### 3. Spatial Coherence

Are the discovered programs spatially structured?

- For each module, compute mean Moran's I of member markers
- Compare distributions across methods
- Module 1-2 should produce higher spatial coherence by design

### 4. Cross-Resolution Consistency

Do the same programs appear at both single-cell and spot level?

- For each method, compute Jaccard overlap between modules found at single-cell vs spot resolution
- Module 1-2 should show higher overlap since pairwise spatial statistics are resolution-agnostic

### Supplementary: top_k Sensitivity Analysis

Sweep `top_k` = {2, 3, 4, 5} across all regions at both resolutions. Report module stability via pairwise Jaccard similarity of discovered modules across parameter values. Expectation: core modules (macrophage, B cell, T cell subtypes) are stable; only edge cases (singletons, weakly co-expressed pairs) change.

---

## Implementation

### Task 1: Marker-to-Lineage Mapping

Create curated dictionary for the 27 Xenium markers.

**File:** `Benchmarking/xenium_benchmarking/benchmark_constants.py` (add to existing)

```python
MARKER_LINEAGE_MAP: Dict[str, str] = {
    # T cell
    "CD3E": "T cell", "CD4": "T cell", "CD8A": "T cell",
    "CD45RO": "T cell", "GranzymeB": "T cell",
    # B cell
    "CD20": "B cell", "CD45RA": "B cell",
    # Myeloid
    "CD68": "Myeloid", "CD163": "Myeloid", "CD16": "Myeloid",
    "CD11c": "Myeloid", "HLA-DR": "Myeloid",
    # Stromal / Mesenchymal
    "alphaSMA": "Stromal", "Vimentin": "Stromal",
    # Epithelial
    "PanCK": "Epithelial", "E-Cadherin": "Epithelial", "Beta-catenin": "Epithelial",
    # Endothelial
    "CD31": "Endothelial",
    # Plasma cell
    "CD138": "Plasma cell",
    # Pan-immune
    "CD45": "Pan-immune",
    # Checkpoint / functional
    "PD-1": "Checkpoint", "PD-L1": "Checkpoint", "LAG-3": "Checkpoint",
    "VISTA": "Checkpoint",
    # Proliferation
    "Ki-67": "Proliferation", "PCNA": "Proliferation",
    # Tumor suppressor
    "PTEN": "Broadly expressed",
}
```

### Task 2: Leiden Baseline Script

**File:** `Benchmarking/xenium_benchmarking/evaluation/src/leiden_baseline_comparison.py`

- Load Xenium protein data (single-cell and pseudo-Visium) for each region
- Run Scanpy: `sc.pp.normalize_total()` → `sc.pp.log1p()` → `sc.pp.pca()` → `sc.pp.neighbors()` → `sc.tl.leiden()` at resolutions [0.3, 0.5, 0.8, 1.0, 1.5]
- Extract cluster signatures via `sc.tl.rank_genes_groups()` (Wilcoxon)
- For each cluster, collect markers with adjusted p-value < 0.05 and log2FC > 0.5
- Save as JSON: `{resolution: {cluster_id: [marker_list]}}` per region per resolution level

### Task 3: Evaluation Script

**File:** `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_discovery_methods.py`

Inputs:
- Module 1-2 output (existing profile discovery results)
- Leiden output (from Task 2)
- Marker-to-lineage mapping (from Task 1)

Computes:
1. Biological coherence: fraction of pure vs mixed modules
2. Rare subtype detection: binary presence/absence table
3. Spatial coherence: Moran's I distributions per method
4. Cross-resolution consistency: Jaccard overlap matrix

Outputs:
- `results/discovery_comparison.json` — all metrics
- Publication figures: grouped bar charts, violin plots, heatmaps

### Task 4: top_k Sensitivity Sweep

Extend Module 1-2 to run with top_k = {2, 3, 4, 5} across all regions at both resolutions. Compute pairwise Jaccard stability of discovered modules. Generate supplementary figure showing module stability across parameter values.

### Task 5: SLURM Submission and Validation

Submit all jobs, verify outputs, generate final comparison figures.

---

## Expected Results

- **Biological coherence:** Module 1-2 produces mostly pure single-lineage modules. Leiden at coarse resolution merges lineages; at fine resolution fragments them into noise clusters.
- **Rare subtypes:** Module 1-2 detects mast cells, smooth muscle, exhausted T cells at single-cell resolution (already demonstrated). Leiden may find some at fine resolution but not consistently.
- **Spatial coherence:** Module 1-2 modules have higher mean Moran's I since spatial statistics are baked into the discovery.
- **Cross-resolution:** Module 1-2 modules are more consistent between single-cell and spot level. Leiden clusters will differ substantially because cluster structure depends on resolution and cell count.
- **top_k sensitivity:** Core modules stable across top_k = 2-5. Only singletons and weak pairs change.

---

## Files

| File | Purpose |
|------|---------|
| `benchmark_constants.py` | Add MARKER_LINEAGE_MAP |
| `evaluation/src/leiden_baseline_comparison.py` | Leiden baseline (new) |
| `evaluation/src/evaluate_discovery_methods.py` | Head-to-head evaluation (new) |
| `evaluation/slurm/run_discovery_comparison.sh` | SLURM job script (new) |
| `evaluation/results/discovery_comparison/` | Output directory |

---

**Last updated:** 2026-01-31
