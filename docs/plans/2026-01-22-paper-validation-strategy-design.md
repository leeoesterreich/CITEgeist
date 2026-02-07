# Paper Validation Strategy: Autodiscovery vs Manual Profiles

**Date**: 2026-01-22
**Status**: Design Discussion

## The Core Question

How do we frame CITEgeist's autodiscovery advantage in the paper while maintaining fair comparisons to other methods?

## Current Benchmarking Approach

**Data**: Xenium pseudo-Visium (465,534 cells aggregated into 7,054 spots)

**Ground Truth**: Cell counts per spot, typed using predefined marker rules

**Comparison Methods**:
- Cell2Location (scRNA-seq reference-based)
- Tangram (scRNA-seq reference-based)
- RCTD (scRNA-seq reference-based)
- Seurat (anchor-based transfer)

**Problem**: All comparisons use artificial cell type definitions that don't match spatial biology.

## Proposed Paper Framing

### Narrative Arc

1. **Standard benchmark shows CITEgeist works** (Figures 3-4)
   - Manual 7-type profiles → comparable performance to other methods
   - Proves the core algorithm is sound

2. **But the standard benchmark is flawed** (Figure 5a)
   - Show PanCK + E-Cadherin don't co-localize (score = 0.42)
   - Show Vimentin + E-Cadherin DO co-localize (score = 0.85)
   - The "ground truth" itself is wrong

3. **Autodiscovery finds what's actually there** (Figure 5b)
   - Cross-region consistency of discovered profiles
   - Biologically interpretable: EMT cells, exhausted T cells, M2 macs
   - These are cell STATES, not just cell TYPES

4. **Single-cell validation proves autodiscovery is right** (Figure 6)
   - This is the key novel analysis

## Single-Cell Validation Strategy

### The Unique Opportunity

We have **single-cell protein expression** for all 465,534 Xenium cells:
- 27 protein markers per cell
- Exact spatial coordinates
- Can validate profile assignments at cellular resolution

### Validation Analyses

#### Analysis 1: Profile Marker Co-expression at Cell Level

For each autodiscovered profile, validate that marker co-expression holds at single-cell level:

```
For EMT profile (Vimentin + E-Cadherin):
  - Count cells where both markers are high
  - Compare to cells where only one marker is high
  - Validate: Do these cells exist? Are they spatially clustered?
```

**Expected result**: EMT cells are real - Vimentin+E-Cadherin double-positive cells exist and cluster spatially.

#### Analysis 2: Soft Assignment Accuracy

For each spot, we know:
- CITEgeist's predicted proportions (continuous)
- Actual cell counts by type

**Key metric**: Which profile set explains more variance?

```python
# For manual 7-type profiles:
variance_explained_manual = compare_proportions(predicted_manual, actual_counts)

# For autodiscovered profiles:
variance_explained_auto = compare_proportions(predicted_auto, actual_counts_retyped)

# Hypothesis: autodiscovery explains MORE variance
```

#### Analysis 3: Cell-Level Assignment Comparison

For each cell, compare what different approaches would assign:

| Approach | Cell Type | Based On |
|----------|-----------|----------|
| Manual protein typing | "Fibroblast" (Vimentin+) | Arbitrary thresholds |
| scRNA-seq transfer | "Stromal" | Reference mapping |
| CITEgeist autodiscovery | "EMT cell" (Vim+ECad+) | Spatial colocalization |

**Key question**: Which assignment captures the cell's full phenotype?

#### Analysis 4: Spatial Coherence Test

For cells assigned to each profile:
- Compute Moran's I of profile scores
- Compare: Are autodiscovered profiles more spatially coherent than manual profiles?

**Expected result**: Autodiscovered profiles have higher spatial coherence because they're derived from spatial patterns.

## Comparison to Other Methods: The Fundamental Limitation

### Why scRNA-seq reference methods miss cell states

1. **Cell2Location, Tangram, RCTD** require scRNA-seq reference with predefined cell types
2. The reference may not contain:
   - EMT transitional cells (captured mid-transition)
   - Spatially-induced states (exhausted T cells at tumor margin)
   - Tissue-specific macrophage polarization

3. **These methods can only find what's in the reference**

### CITEgeist's Advantage

- Uses spatial protein data directly
- Discovers what's actually present in the tissue
- Captures cell states that wouldn't survive dissociation for scRNA-seq

## Proposed Figure Layout

### Figure 5: The Ground Truth Problem

**5a**: Colocalization heatmap showing:
- PanCK vs E-Cadherin: weak (0.42)
- Vimentin vs E-Cadherin: strong (0.85)
- "Epithelial" definition is wrong

**5b**: Autodiscovered profiles across 5 regions
- Dendrogram showing consistent clustering
- Profile composition per region

**5c**: Biological interpretation
- EMT profile → metastatic potential
- Exhausted T → immunotherapy response
- M2 Mac → tumor-promoting environment

### Figure 6: Single-Cell Validation

**6a**: Single-cell marker co-expression
- Scatter: Vimentin vs E-Cadherin at cell level
- Highlight EMT population (double positive)

**6b**: Spatial map of EMT cells
- Show they cluster at specific tissue locations
- Compare to "Fibroblast" assignment (scattered)

**6c**: Proportion accuracy comparison
- Manual profiles vs autodiscovered
- Variance explained metric

### Figure 7: Method Comparison with Revised Ground Truth

**7a**: Re-evaluate all methods with spatially-validated GT
- CITEgeist + autodiscovery
- CITEgeist + manual profiles
- Cell2Location
- Tangram
- RCTD/Seurat

**7b**: Key insight: Reference-free method finds more biology

## Implementation Plan

### Step 1: Single-Cell Profile Scoring

```python
def score_cells_by_profile(cell_protein_matrix, profiles):
    """
    Score each cell for each autodiscovered profile.

    Args:
        cell_protein_matrix: (n_cells, 27) protein expression
        profiles: Dict[profile_name, List[marker_names]]

    Returns:
        scores: (n_cells, n_profiles) profile scores
    """
    # For each profile, compute mean expression of profile markers
    # Normalize across profiles for soft assignment
```

### Step 2: Compare to Spot-Level Predictions

```python
def validate_spot_proportions(cell_scores, cell_to_spot_mapping, spot_predictions):
    """
    Compare cell-level scores aggregated to spots vs CITEgeist predictions.

    Metric: Pearson correlation between:
    - Aggregated cell scores per spot
    - CITEgeist proportion predictions
    """
```

### Step 3: Cross-Method Comparison

```python
def compare_cell_assignments(cell_scores_citegeist,
                             cell_types_manual,
                             cell_types_scrnaseq_transfer):
    """
    For cells assigned differently, examine which captures true phenotype.

    Focus on: Cells in EMT profile vs "Fibroblast" vs "Epithelial"
    """
```

## Key Messages for Paper

1. **Standard benchmarks are necessary but insufficient**
   - CITEgeist performs comparably on artificial ground truth
   - This proves the algorithm works

2. **The ground truth itself can be wrong**
   - Marker colocalization reveals biological reality
   - Some predefined cell types don't exist as defined

3. **Autodiscovery captures cell states, not just types**
   - EMT cells are transitional states
   - Reference-based methods can't find what's not in the reference

4. **Single-cell validation proves this isn't overfitting**
   - Autodiscovered profiles exist at the cell level
   - They're spatially coherent
   - They capture more biological variance

## Next Steps

1. [ ] Implement cell-level profile scoring
2. [ ] Generate single-cell co-expression plots
3. [ ] Compute spatial coherence metrics
4. [ ] Create revised ground truth based on spatial profiles
5. [ ] Re-run all method comparisons with revised GT
6. [ ] Generate Figure 5-7 visualizations
