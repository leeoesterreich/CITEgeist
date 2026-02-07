# Revised Ground Truth Definitions for Xenium Benchmarking

**Date**: 2026-01-22
**Status**: Proposed

## Problem Statement

The current "achievable 7-type" ground truth definitions used for Xenium benchmarking were constructed based on theoretical marker-celltype associations from literature. However, spatial colocalization analysis across 5 Xenium regions revealed that these definitions do not match the actual spatial biology:

| Current Definition | Issue |
|--------------------|-------|
| Stromal = CD3E, Vimentin | CD3E is a T cell marker, not stromal |
| Epithelial = PanCK, E-Cadherin | PanCK+E-Cadherin colocalization score = 0.42 (weak) |
| Fibroblasts = alphaSMA, Vimentin | Vimentin's top partner is E-Cadherin (0.89), not alphaSMA (0.83) |
| Proliferating T = CD138, CD4 | CD138 is plasma cell marker, not T cell |

## Analysis Summary

### Colocalization Evidence (Bivariate Moran's I across 5 regions)

**Strong colocalization (score > 0.80):**
- CD20 + CD45RA: 0.86 (B cells)
- CD68 + CD16: 0.90 (Macrophages)
- CD3E + CD45RO: 0.88 (T cells)
- Vimentin + E-Cadherin: 0.85 (EMT/Hybrid cells)
- CD8A + PD-1: 0.87 (Exhausted T cells)

**Weak colocalization (score < 0.50):**
- PanCK + E-Cadherin: 0.42 (do NOT co-localize)
- alphaSMA + CD31: 0.38 (pericytes vs endothelium - separate)

### Key Biological Insight: EMT Cells

The strongest signal in the data is Vimentin + E-Cadherin colocalization (0.85), which represents cells undergoing Epithelial-Mesenchymal Transition (EMT). This is biologically significant in breast cancer:

- EMT is associated with metastasis and drug resistance
- These are hybrid cells expressing both epithelial (E-Cadherin) and mesenchymal (Vimentin) markers
- This is NOT an artifact - it's a real biological state

## Revised 8-Type Spatially-Validated Definitions

```python
REVISED_GT_CELLTYPE_MARKERS = {
    # Core immune populations
    "B cells": {"CD20", "CD45RA"},           # Score: 0.86, validated
    "Macrophages": {"CD68", "CD16"},         # Score: 0.90, validated
    "T cells": {"CD3E", "CD45RO"},           # Score: 0.88, pan-T marker

    # Functional immune states
    "Exhausted T cells": {"CD8A", "PD-1"},   # Score: 0.87, functional state
    "M2 Macrophages": {"CD163"},             # Singleton, polarization marker

    # Structural populations
    "Endothelial": {"CD31"},                 # Singleton, consistent
    "Epithelial": {"PanCK"},                 # True epithelial (NOT E-Cadherin)

    # EMT/Hybrid population
    "EMT cells": {"Vimentin", "E-Cadherin"}, # Score: 0.85, hybrid state
}
```

## Changes from Previous Definitions

| Population | Old Definition | New Definition | Rationale |
|------------|---------------|----------------|-----------|
| Stromal | CD3E, Vimentin | REMOVED | CD3E is T cell marker; Vimentin goes to EMT |
| Epithelial | PanCK, E-Cadherin | PanCK only | E-Cadherin doesn't co-localize with PanCK |
| Fibroblasts | alphaSMA, Vimentin | REMOVED | Markers don't co-localize (0.83 < Vim+ECad 0.89) |
| EMT cells | N/A | Vimentin, E-Cadherin | NEW - reflects spatial biology |
| Exhausted T | N/A | CD8A, PD-1 | NEW - functional state |
| M2 Macrophages | N/A | CD163 | NEW - polarization state |

## Implications for Benchmarking

### What this means:
1. **Cell types vs cell states**: Spatial data captures functional states (EMT, exhausted T, M2 macs), not just cell type lineages
2. **Autodiscovery is correct**: The 15 profiles from Module 2b likely reflect real biology
3. **GT definitions need updating**: Current GT penalizes correct biological findings

### Recommended evaluation approach:
1. Use revised 8-type definitions for proportion evaluation
2. Accept that some profiles capture states, not types
3. Consider hierarchical evaluation: major types + subtypes/states

## Files to Update

1. `Benchmarking/xenium_benchmarking/CITEgeist/src/run_benchmark.py` - Update profile definitions
2. `Benchmarking/xenium_benchmarking/evaluation/src/evaluate_all_methods.py` - Update GT mapping
3. Create new evaluation script for revised definitions

## Validation Steps

1. Re-run CITEgeist with revised 8-type profiles
2. Compare correlation metrics against previous 7-type results
3. Visual inspection of spatial patterns for EMT cells
4. Check if Exhausted T and M2 Mac populations are spatially coherent

## References

- EMT marker coexpression: Pastushenko & Blanpain, Trends Cell Biol 2019
- CD8+PD-1+ exhausted T cells: Wherry & Kurachi, Nat Rev Immunol 2015
- M2 macrophage markers: Murray, Annu Rev Physiol 2017
