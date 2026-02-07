# Figure Audit and Gap Analysis

**Date**: 2026-02-06
**Auditor**: Claude (overnight task)
**Purpose**: Identify gaps and quality issues in Figures 1-6 for Cell/Cell Reports publication

---

## Summary

| Figure | Current Panels | Quality Issues | Critical Gaps |
|--------|---------------|----------------|---------------|
| 1 | A-C (Pipeline) | Font sizes small (9pt), loose layout | Data thumbnails missing |
| 2 | A-D (Profile Discovery) | Font sizes 8pt, some overlap | Cell resolution data sparse |
| 3 | A-B (Benchmarking) | Improved layout, 11pt fonts | No scatter plot (Panel C removed), no simulated benchmark |
| 4 | A-D (Programs) | Good layout, 11pt fonts | Bivariate heatmap missing |
| 5 | A-D (Integration) | Good layout | UMAP legend crowded |
| 6 | A-B (Interoperability) | Sparse content | Missing vignettes, wet lab validation panels |

---

## Figure 1: Pipeline Overview

### Current Panels
- **A**: Module flow schematic (5 boxes with arrows)
- **B**: Spatial operations (Moran's I formula, colocalization, spatial graph)
- **C**: Resolution flexibility (Visium vs single-cell)

### Quality Issues
1. **Font sizes**: Base 9pt is too small for print
2. **Whitespace**: Excessive padding between panels
3. **Colors**: Default matplotlib colors, not uniform palette
4. **Legend**: Cell type legend in Panel C uses small 6pt font

### Missing Content
- No actual data thumbnails (all schematic)
- Could add mini spatial heatmaps showing real output
- Input/output labels very small (7pt)

### Recommended Changes
- [ ] Increase base font to 10pt minimum
- [ ] Apply consistent PALETTE from figure_style.py
- [ ] Tighten layout with constrained_layout
- [ ] Add data preview thumbnails in appropriate panels

---

## Figure 2: Modules 1-2 Profile Discovery

### Current Panels
- **A**: Marker interest detection (kurtosis, GMM, Moran's I gates)
- **B**: Colocalization -> clustering -> profiles workflow
- **C**: Xenium single-cell demonstration summary
- **D**: Discovered profiles vs known markers table

### Quality Issues
1. **Font sizes**: Base 8pt is below minimum
2. **Panel A**: Distribution sketches are synthetic, not real data
3. **Panel C**: Summary box is text-heavy, low data density
4. **Panel D**: Table layout has some crowding

### Missing Content
- Real marker histograms with kurtosis values
- Actual colocalization heatmap (not just schematic)
- Dendrogram with real marker names
- P-values on colocalization comparisons

### Recommended Changes
- [ ] Increase fonts to 10pt minimum
- [ ] Replace schematic distributions with real data
- [ ] Add actual colocalization matrix visualization
- [ ] Include p-value annotations
- [ ] Tighten Panel C text into bullets

---

## Figure 3: Benchmarking

### Current Panels
- **A**: Deconvolution schematic (simplified, clean)
- **B**: Three bar charts (Correlation, JSD, RMSE)

### Quality Issues
1. **Panel C removed**: Scatter plot no longer in figure
2. **Bar charts**: Method labels rotated 45 degrees (acceptable)
3. **Error bars**: Present and properly styled

### Missing Content (from v3)
- **Scatter plot**: Predicted vs ground truth visualization
- **Simulated benchmark**: Results from `Benchmarking/simulation_benchmarking/Figures/`
- **Spatial maps**: Per-method proportion visualization
- Simulation benchmark data available in CSV files

### Recommended Changes
- [ ] Add Panel C: Predicted vs GT scatter plot
- [ ] Add Panel D: Simulated benchmark summary (supporting role)
- [ ] Or: Replace with spatial proportion maps comparison
- [ ] Ensure consistent method colors across panels

---

## Figure 4: Module 4 Spatial Programs

### Current Panels
- **A**: NMF schematic (clean, good size)
- **B**: Top programs horizontal bar chart with genes
- **C**: Moran's I boxplot by cell type
- **D**: Summary statistics table

### Quality Issues
1. **Panel B**: Gene labels may overlap at small widths
2. **Panel D**: Table could be more compact
3. Good font sizes (11pt base)

### Missing Content
- **Bivariate relationships heatmap**: Module 4b output
- **Spatial plots**: Program activity maps with Moran's I annotation
- Real expression heatmaps showing top genes per program

### Recommended Changes
- [ ] Add bivariate relationship heatmap (replace or supplement Panel D)
- [ ] Include spatial activity maps for top programs
- [ ] Annotate Moran's I values on spatial plots
- [ ] Check gene label legibility at 89mm width

---

## Figure 5: Module 5 Integration

### Current Panels
- **A**: Integration workflow schematic
- **B**: UMAP of programs colored by cell type
- **C**: Response-associated programs bar chart
- **D**: Volcano plot (PyDESeq2)

### Quality Issues
1. **Panel B legend**: 6+ cell types in small legend, may crowd
2. **Schematic**: Uses hard-coded sample colors
3. Good font sizes (11pt base)

### Missing Content
- Conserved relationship network (edge weights visible)
- More detail on aligned program count (590 programs mentioned in plan)
- 127 DE genes mentioned but volcano may not show all labels

### Recommended Changes
- [ ] UMAP: Annotate program count clearly
- [ ] Volcano: Ensure top gene labels don't overlap
- [ ] Add p-value threshold markers on volcano
- [ ] Consider adding relationship network panel

---

## Figure 6: Interoperability & Validation

### Current Panels
- **A**: Workflow diagram (CITEgeist outputs -> tools)
- **B**: Output formats table

### Quality Issues
1. **Panel A**: Only 2 columns of content, sparse
2. **Panel B**: Table is basic, could be more informative
3. Missing actual demonstration data

### Missing Content (Critical)
- **Vignette outputs**: Real patient examples from notebooks
  - Files exist: `CITEgeist/examples/vignette_*.ipynb`
- **Wet lab validation**: Midkine discovery
  - Prism files: `midkine/CITEgeistImageAnalysis.prism`, `midkine/ELISA Results.prism`
  - Pipeline outputs: `midkine/mdk_saturation_pipeline/outputs/figures/`
- **ELISA quantification**: p-values available
- **Image analysis**: Spatial patterns from CITEgeist

### Available Validation Data
```
midkine/mdk_saturation_pipeline/outputs/
├── figures/
│   ├── fig10_unsupervised_hsp90b1.png
│   ├── fig11_venn_expression_binding.png
│   └── ... (16 figures total)
└── tables/
    ├── evidence_summary.csv
    ├── mechanism_discrimination_scorecard.csv
    └── ... (10 CSV files)
```

### Recommended Changes
- [ ] Add Panel C: Vignette example outputs (real spatial patterns)
- [ ] Add Panel D: Wet lab validation summary
  - Use figures from `midkine/mdk_saturation_pipeline/outputs/figures/`
  - Include ELISA results summary
- [ ] Rework Panel A to be denser
- [ ] Tighten Panel B table

---

## Cross-Figure Issues

### Style Consistency
- **Font families**: All use 'sans-serif' (OK)
- **Font sizes**: Vary between 8-11pt base (need standardization to 10pt minimum)
- **Cell type colors**: Different palettes across figures
  - Figure 1: `#e74c3c`, `#3498db`, `#2ecc71`, `#9b59b6`, `#f39c12`
  - Figure 3: `#E41A1C`, `#377EB8`, `#4DAF4A`, `#984EA3`, `#FF7F00`
  - Need CELL_TYPE_COLORS standardization

### Color Accessibility
- Some color combinations may not be color-blind safe
- Need to verify with accessibility checkers

### Data Density
- Figures 1, 2, 6 have more schematic content than data
- Figures 3, 4, 5 have good data density

---

## Priority Actions

### High Priority (Publication Blocking)
1. Create `figure_style.py` with unified PALETTE and CELL_TYPE_COLORS
2. Fix font sizes to minimum 10pt across all figures
3. Add Panel C (scatter) and Panel D (simulated) to Figure 3
4. Add wet lab validation panels to Figure 6

### Medium Priority (Quality Enhancement)
5. Replace schematic distributions with real data in Figure 2
6. Add bivariate heatmap to Figure 4
7. Tighten whitespace across all figures
8. Add Moran's I annotations to spatial plots

### Low Priority (Polish)
9. Add data thumbnails to Figure 1
10. Improve legend readability in Figure 5
11. Color accessibility review

---

## Data Sources Identified

| Figure | Data Path | Status |
|--------|-----------|--------|
| 2 | `Benchmarking/xenium_benchmarking/CITEgeist/output_singlecell/` | Available |
| 3 | `Benchmarking/xenium_benchmarking/evaluation/results/method_comparison/` | Available |
| 3 | `Benchmarking/simulation_benchmarking/Figures/*.csv` | Available |
| 4 | `Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/` | Available |
| 5 | `examples/output_module5_pydeseq/` | Available |
| 6 | `midkine/mdk_saturation_pipeline/outputs/` | Available |
| 6 | `CITEgeist/examples/vignette_*.ipynb` | Available |

---

*Generated by overnight task execution - 2026-02-06*
