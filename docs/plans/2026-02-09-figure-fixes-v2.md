# Figure Fixes Plan v2 — 2026-02-09

## Issues Reported

### Figure 1
- [ ] **1B**: Bottom text overlapping gray box - need more vertical space

### Figure 2
- [ ] **2A**: Right edge cutoff - SVG schematic content extends past canvas

### Figure 3
- [ ] **3A**: Overlapping error between deconvolution schematic and profile discovery table
- [ ] **Missing**: GEX benchmark for Xenium data (currently only simulated GEX)
- [ ] **3F**: May need more comprehensive cell type comparison

### Figure 4
- [ ] **4D**: Pathway dot plot extending too far right - needs horizontal shrink

### Figure 5 (Major Redesign)
- [ ] **5B**: Not clear what it shows - needs redesign for better information conveyance
- [ ] **5D**: Network needs expanded co-localized groups with labels explaining connections
- [ ] **5E/F**: Need clearer text and colors showing co-localization/exclusion phenotype
- [ ] **All spatial panels (E, F, H)**: Currently using direct scatter - should use sc.pl.spatial() with H&E backdrop
- [ ] **Colors**: Switch to viridis or RdBu diverging colormap throughout

---

## Fixes by Figure

### Figure 1 Fixes

**Issue**: Panel B bottom text overlapping gray box

**Fix**:
1. Increase figure height from (12, 11) to (12, 12)
2. Increase hspace from 0.25 to 0.30
3. Or: Modify the SVG schematic to have more bottom padding

**File**: `generate_figure1.py` line 46-47

---

### Figure 2 Fixes

**Issue**: Panel A right edge cutoff

**Fix**:
1. The script already increased figsize to (14, 10) and wspace to 0.18
2. Check if SVG schematic `figure2_panel_a_marker_interest.png` itself is clipped
3. May need to regenerate the SVG with more canvas width or adjust viewBox

**File**: `svg_schematics.py` or the source SVG file

---

### Figure 3 Fixes

**Issue 1**: Panel A overlap between table and deconvolution schematic

**Current code** (lines 232-281):
- Table bbox: `[0.0, 0.05, 0.58, 0.88]` (60% width)
- Schematic inset: `[0.62, 0.15, 0.36, 0.70]` (38% width, 4% gap)

**Fix**:
- Reduce table width to 55%: `[0.0, 0.05, 0.55, 0.88]`
- Move schematic right: `[0.58, 0.15, 0.40, 0.70]`
- Or shrink schematic and add padding

**Issue 2**: Missing Xenium GEX benchmark

**Reality check**: Xenium benchmark only has proportion ground truth (from cell segmentation), NOT GEX ground truth. The GEX benchmark can only be done on simulated data where we know true cell-type-specific expression.

**Decision**: This is not a bug - Xenium doesn't have GEX ground truth. Add clarifying text to panel D.

**Issue 3**: Panel F comprehensive comparison

**Current**: Already shows all 6 cell types in 2x6 grid

**Potential improvement**: Increase subplot sizes for better visibility

---

### Figure 4 Fixes

**Issue**: Panel D pathway dot plot extending too far right

**Fix**:
- Tighten x-axis limits in `panel_d_pathway_dotplot()`
- Reduce right padding
- Possibly shrink dot sizes or move legend inside

**File**: `generate_figure4.py`

---

### Figure 5 Fixes (Major)

**Issue 1**: Panel B unclear

**Current**: Horizontal bar chart showing # aligned programs per cell type
**Problem**: Doesn't convey conservation patterns well

**Redesign options**:
1. Heatmap: programs (rows) x samples (cols), colored by presence/intensity
2. Stacked bar: show breakdown by conservation level
3. Better annotations explaining what "aligned programs" means

**Issue 2**: Panel D network needs expansion

**Current**: NetworkX spring layout with colored edges
**Problem**: Hard to understand which nodes are connected and why

**Fix**:
1. Increase canvas size for network panel
2. Add node labels for ALL key relationship participants
3. Draw co-localized groups as clusters with background shading
4. Add arrows or annotations explaining "Fibroblast-CD4 T co-localization in responders"

**Issue 3**: Panels E/F need clearer phenotype illustration

**Current**: Two side-by-side spatial maps + scatter inset
**Problems**:
- Colors not distinct enough
- Text annotation at bottom may be missed
- Co-localization vs exclusion not visually obvious

**Fixes**:
1. Use diverging colormap (RdBu) for proportion difference
2. Add arrows or annotations pointing to key regions
3. Make the scatter inset larger with clearer regression line
4. Add explicit "High overlap" vs "Non-overlapping" annotations

**Issue 4**: Use sc.pl.spatial() with H&E

**Current status**:
- Panels E/F: Already use `sc.pl.spatial()` with `alpha_img=0.7`
- Panel H: Already uses `sc.pl.spatial()` with `alpha_img=0.6`

**Verification needed**: Check if H&E is actually loading. The images showed mostly white backgrounds which suggests either:
1. `load_images=True` failing silently
2. `alpha_img` too low
3. Image not in expected location

**Fix**:
1. Debug image loading in `load_visium_with_image()`
2. Increase `alpha_img` to 0.9
3. Verify Visium outs folder has `spatial/tissue_hires_image.png`

**Issue 5**: Color scheme

**Current**: Using viridis in many places, some custom colormaps
**Fix**: Standardize on viridis for sequential data, RdBu for diverging (proportion differences)

---

## Implementation Order

1. **Quick fixes** (can do immediately):
   - Figure 1B: Increase hspace
   - Figure 3A: Adjust table/schematic positioning
   - Figure 4D: Tighten x-axis

2. **Medium fixes** (require careful testing):
   - Figure 5 H&E loading verification
   - Figure 5 color scheme standardization
   - Figure 5E/F annotation improvements

3. **Major redesign** (need discussion):
   - Figure 5B: What should it actually show?
   - Figure 5D: Network expansion - how much detail?

---

## Questions for User

1. **Figure 5B**: The current bar chart shows # programs per cell type. What information would be more useful? Options:
   - Heatmap of program presence across samples
   - Stacked bar showing conservation levels
   - Something else?

2. **Figure 5D network**: How much annotation is desired?
   - Just the key 4 relationships labeled?
   - All co-localized relationships expanded?
   - Background shading for clusters?

3. **H&E visibility**: The current `alpha_img=0.7` may be too transparent. Should we increase to 0.9 or keep some transparency to see spots better?
