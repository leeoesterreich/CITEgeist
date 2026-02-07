# Pipeline Overview Figure Design

**Date**: 2026-01-28
**Status**: Complete

## Summary

High-level pipeline overview SVG for the CITEgeist paper. Horizontal left-to-right flow showing Modules 1-5 as colored boxes with arrows indicating data flow between stages.

## Specifications

- **Format**: SVG vector graphic (editable in Illustrator/Inkscape)
- **Layout**: Horizontal, 5 module boxes, left-to-right
- **Canvas**: 1400 x 320 px
- **Style**: Graphical abstract / methods figure

## Module Content

| Module | Title | Description | Color |
|--------|-------|-------------|-------|
| 1 | Marker Interest Detection | Filter antibody markers by spatial signal | `#5B7FA5` (slate blue) |
| 2 | Profile Assembly | Discover cell type marker profiles | `#7BA58A` (sage green) |
| 3 | Deconvolution | Estimate proportions & deconvolve gene expression | `#C4965A` (warm amber) |
| 4 | Program Discovery | NMF-based spatial gene programs per cell type | `#B07A8F` (dusty rose) |
| 5 | Cross-Sample Integration | Align programs across patients | `#8A7EB5` (steel purple) |

## Arrow Labels (data flow)

1. Module 1 -> 2: "Interesting markers"
2. Module 2 -> 3: "Cell type profiles"
3. Module 3 -> 4: "Deconvolved layers"
4. Module 4 -> 5: "Programs + relationships"

## Files

- Generator script: `figures/generate_pipeline_overview.py`
- Output SVG: `figures/pipeline_overview.svg`

## Usage

```bash
python figures/generate_pipeline_overview.py
```

The SVG can be opened in Inkscape or Adobe Illustrator for further refinement before submission.
