# Overnight Figure Enhancement - Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.
> **For each figure task:** Use ralph-loop with the specified completion promise.

**Goal:** Elevate Figures 1-6 to Cell/Cell Reports publication quality through systematic review and iterative refinement.

**Design Doc:** `docs/plans/2026-02-06-figure-enhancement-design.md`

**Branch:** overnight/2026-02-06-figure-enhancement

---

## Shared Style Configuration

Before starting figure tasks, create this shared style module:

```python
# manuscript/figures/figure_style.py
import matplotlib.pyplot as plt

PALETTE = {
    'primary': '#2171b5',
    'secondary': '#6baed6',
    'accent1': '#fe9929',
    'accent2': '#41ab5d',
    'neutral': '#636363',
    'highlight': '#c51b8a',
}

CELL_TYPE_COLORS = {
    'Epithelial': '#2171b5',
    'Fibroblasts': '#fe9929',
    'Macrophages': '#41ab5d',
    'T cells': '#c51b8a',
    'CD4+ T': '#c51b8a',
    'CD8+ T': '#ae017e',
    'B cells': '#6baed6',
    'Endothelial': '#636363',
}

def apply_style():
    plt.rcParams.update({
        'font.family': 'Arial',
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'figure.constrained_layout.use': True,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })
```

---

## Task 1: Audit & Gap Analysis

**Files:**
- Read: `manuscript/figures/generate_figure*.py` (all 6)
- Read: `manuscript/CITEgeistManuscript_v3.docx` or related v3 files
- Read: `midkine/*.prism` exports (check for CSV exports)
- Create: `FIGURE_AUDIT.md`

**Steps:**

1. Read all 6 figure generation scripts, note current panel structure
2. Cross-reference v3 manuscript to identify:
   - Simulated benchmark content that should appear (supporting role)
   - Vignette figures that got dropped
   - Wet lab validation panels missing
3. Check `midkine/` for Prism CSV exports or analysis outputs
4. For each figure, document in `FIGURE_AUDIT.md`:
   - Current panels
   - Missing content from v3
   - Specific quality issues (text overlap, sparse panels, etc.)
   - Recommended additions

**Verification:**
```bash
cat FIGURE_AUDIT.md | head -100
```
Expected: Shows structured gap analysis for all 6 figures

**Commit:** `docs: create figure audit with gap analysis`

---

## Task 2: Create Shared Style Module

**Files:**
- Create: `manuscript/figures/figure_style.py`

**Steps:**

1. Create the style module with PALETTE and apply_style() function (see above)
2. Test import works

**Verification:**
```bash
python -c "from manuscript.figures.figure_style import apply_style, PALETTE; print('OK')"
```
Expected: Prints "OK"

**Commit:** `fig: add shared style module for Cell/Cell Reports aesthetic`

---

## Task 3: Figure 1 Enhancement (Pipeline Overview)

**Files:**
- Modify: `manuscript/figures/generate_figure1.py`
- Output: `manuscript/figures/output/figure1_pipeline_overview.pdf`

**Use Ralph Loop:**
```
/ralph-loop "Review and enhance Figure 1 (Pipeline Overview).

Current state: Read generate_figure1.py, understand panel structure.

For EACH panel, check against these criteria:
1. TEXT OVERLAP: Any labels/annotations overlapping? Fix with constrained_layout, adjust positions
2. FONT SIZE: All text ≥10pt? Increase if smaller
3. WHITESPACE: Excessive empty space? Tighten with pad=0.3, reduce margins
4. COLORS: Using gaudy/default matplotlib colors? Replace with PALETTE from figure_style.py
5. DATA DENSITY: Does panel show real data or just schematic? Add data thumbnails if schematic-only
6. LEGIBILITY: Readable at 89mm width? Test by viewing at small size

After checking all panels:
- Apply fixes to generate_figure1.py
- Import and use figure_style.apply_style() at top
- Regenerate: python manuscript/figures/generate_figure1.py
- Re-check all criteria
- Log what was fixed in OVERNIGHT_LOG.md

Iterate until ALL panels pass ALL criteria." \
--completion-promise "Figure 1 all panels pass publication criteria" \
--max-iterations 10
```

**On loop completion:**
**Commit:** `fig(1): enhance Figure 1 pipeline overview to publication quality`

---

## Task 4: Figure 2 Enhancement (Modules 1-2 Profile Discovery)

**Files:**
- Modify: `manuscript/figures/generate_figure2.py`
- Output: `manuscript/figures/output/figure2_profile_discovery.pdf`

**Use Ralph Loop:**
```
/ralph-loop "Review and enhance Figure 2 (Modules 1-2 Profile Discovery).

Current state: Read generate_figure2.py, understand panel structure.

For EACH panel, check criteria:
1. TEXT OVERLAP: Fix overlapping labels
2. FONT SIZE: Minimum 10pt
3. WHITESPACE: Aggressive minimization
4. COLORS: Use PALETTE, ensure color-blind safe
5. DATA DENSITY: Show real colocalization data, dendrograms with actual labels
6. STATISTICAL: Add p-values to marker comparisons if applicable

Key content to ensure present:
- Marker interest detection with real spatial patterns
- Colocalization network/heatmap with actual marker names
- Dendrogram showing profile discovery
- Comparison of discovered vs known markers

Apply figure_style.apply_style(), regenerate, re-check, iterate." \
--completion-promise "Figure 2 all panels pass publication criteria" \
--max-iterations 10
```

**On loop completion:**
**Commit:** `fig(2): enhance Figure 2 profile discovery to publication quality`

---

## Task 5: Figure 3 Enhancement (Benchmarking)

**Files:**
- Modify: `manuscript/figures/generate_figure3.py`
- Output: `manuscript/figures/output/figure3_benchmarking.pdf`
- Reference: `Benchmarking/simulation_benchmarking/` (simulated benchmark - supporting role)
- Reference: `Benchmarking/xenium_benchmarking/evaluation/results/`

**Use Ralph Loop:**
```
/ralph-loop "Review and enhance Figure 3 (Benchmarking).

Current state: Read generate_figure3.py, check what benchmark data is shown.

CRITICAL ADDITIONS from v3:
- Simulated benchmark results should appear in supporting role (check Benchmarking/simulation_benchmarking/)
- Xenium benchmark is primary focus
- Per-method spatial proportion maps (not just bar charts)

For EACH panel, check criteria:
1. TEXT OVERLAP: Fix overlapping method names, axis labels
2. FONT SIZE: Minimum 10pt
3. WHITESPACE: Pack more data, reduce margins
4. COLORS: Consistent method colors across panels
5. DATA DENSITY: Show spatial maps alongside metrics
6. STATISTICAL: Error bars, significance markers on comparisons

Content structure:
- Panel A: Deconvolution schematic (keep but tighten)
- Panel B: Method comparison metrics (JSD, RMSE, correlation) with error bars
- Panel C: Spatial proportion visualization per method
- Panel D (if space): Simulated benchmark summary (supporting)

Apply style, regenerate, re-check, iterate." \
--completion-promise "Figure 3 all panels pass publication criteria" \
--max-iterations 10
```

**On loop completion:**
**Commit:** `fig(3): enhance Figure 3 benchmarking to publication quality`

---

## Task 6: Figure 4 Enhancement (Module 4 Programs)

**Files:**
- Modify: `manuscript/figures/generate_figure4.py`
- Output: `manuscript/figures/output/figure4_module4_programs.pdf`
- Reference: `Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/`

**Use Ralph Loop:**
```
/ralph-loop "Review and enhance Figure 4 (Module 4 Spatial Programs).

Current state: Read generate_figure4.py, check panel structure.

For EACH panel, check criteria:
1. TEXT OVERLAP: Gene names in heatmaps must be readable
2. FONT SIZE: Minimum 10pt for all labels
3. WHITESPACE: Tight layout
4. COLORS: Consistent, accessible colormap for expression
5. DATA DENSITY: Show multiple programs per cell type
6. STATISTICAL: Moran's I values annotated on spatial plots

Key content:
- NMF schematic (keep but tighten)
- Program heatmap showing top genes per program
- Spatial plots of program activity with Moran's I annotation
- Bivariate relationship heatmap

Apply style, regenerate, re-check, iterate." \
--completion-promise "Figure 4 all panels pass publication criteria" \
--max-iterations 10
```

**On loop completion:**
**Commit:** `fig(4): enhance Figure 4 spatial programs to publication quality`

---

## Task 7: Figure 5 Enhancement (Module 5 Integration)

**Files:**
- Modify: `manuscript/figures/generate_figure5.py`
- Output: `manuscript/figures/output/figure5_integration.pdf`
- Reference: `examples/output_module5_pydeseq/`

**Use Ralph Loop:**
```
/ralph-loop "Review and enhance Figure 5 (Module 5 Cross-Sample Integration).

Current state: Read generate_figure5.py, check panel structure.

For EACH panel, check criteria:
1. TEXT OVERLAP: UMAP labels, volcano gene names
2. FONT SIZE: Minimum 10pt
3. WHITESPACE: Aggressive minimization
4. COLORS: Cell type colors consistent with other figures
5. DATA DENSITY: Pack information into UMAP (590 programs), volcano (127 DE genes)
6. STATISTICAL: p-value thresholds clearly marked on volcano

Key content:
- Integration schematic (tighten)
- UMAP of aligned programs (colored by cell type, annotated count)
- Responder vs progressor comparison
- Volcano plot with top genes labeled
- Conserved relationship network (edge weights visible)

Apply style, regenerate, re-check, iterate." \
--completion-promise "Figure 5 all panels pass publication criteria" \
--max-iterations 10
```

**On loop completion:**
**Commit:** `fig(5): enhance Figure 5 integration to publication quality`

---

## Task 8: Figure 6 Enhancement (Interoperability)

**Files:**
- Modify: `manuscript/figures/generate_figure6.py`
- Output: `manuscript/figures/output/figure6_interoperability.pdf`
- Reference: `midkine/CITEgeistImageAnalysis.prism` (check for CSV)
- Reference: `midkine/ELISA Results.prism` (check for CSV)
- Reference: `midkine/mdk_saturation_pipeline/outputs/`
- Reference: `CITEgeist/examples/vignette_*.ipynb`

**Use Ralph Loop:**
```
/ralph-loop "Review and enhance Figure 6 (Interoperability & Validation).

Current state: Read generate_figure6.py, check panel structure.

CRITICAL ADDITIONS from v3:
- Vignette real outputs (patient examples from notebooks)
- Wet lab validation panels (check midkine/ for Prism exports or processed data)
- Midkine discovery pathway visualization

For EACH panel, check criteria:
1. TEXT OVERLAP: Fix all overlapping annotations
2. FONT SIZE: Minimum 10pt
3. WHITESPACE: Pack validation data densely
4. COLORS: Professional, consistent
5. DATA DENSITY: Show actual experimental data, not placeholders
6. STATISTICAL: Include p-values from wet lab validation

Content structure:
- Panel A: Workflow diagram (CITEgeist → external tools)
- Panel B: Vignette example outputs (real spatial patterns)
- Panel C: Midkine discovery summary
- Panel D: Wet lab validation (ELISA, image quantification)

Check midkine/ folder for:
- CSV exports from Prism files
- Python analysis outputs in mdk_saturation_pipeline/outputs/
- Any processed validation data

Apply style, regenerate, re-check, iterate." \
--completion-promise "Figure 6 all panels pass publication criteria" \
--max-iterations 10
```

**On loop completion:**
**Commit:** `fig(6): enhance Figure 6 interoperability to publication quality`

---

## Task 9: Cross-Figure Style Consistency

**Files:**
- Review: All 6 figure scripts
- Check: All output PDFs

**Steps:**

1. Verify all figures import and use `figure_style.apply_style()`
2. Check color consistency:
   - Same cell type → same color across figures
   - Method colors consistent in benchmarking
3. Verify font sizes consistent
4. Check legend styles match (inside panel, no border)
5. Verify all outputs regenerated with today's date

**Verification:**
```bash
ls -la manuscript/figures/output/*.pdf
```
Expected: All 6 PDFs with today's date

**Commit:** `fig: verify cross-figure style consistency`

---

## Task 10: Final Verification

**Use:** `/superpowers:verification-before-completion`

**Checklist:**
- [ ] All 6 figures regenerated
- [ ] No text overlap in any panel
- [ ] All fonts ≥10pt
- [ ] Whitespace minimized
- [ ] Colors consistent and accessible
- [ ] V3 content migrated (simulated benchmark, vignettes, wet lab)
- [ ] FIGURE_AUDIT.md documents changes made
- [ ] All commits have descriptive messages

**Final Commit:** `fig: complete overnight figure enhancement`

---

## Summary

| Task | Description | Ralph Loop | Max Iter |
|------|-------------|------------|----------|
| 1 | Audit & gap analysis | No | - |
| 2 | Create style module | No | - |
| 3 | Figure 1 enhancement | Yes | 10 |
| 4 | Figure 2 enhancement | Yes | 10 |
| 5 | Figure 3 enhancement | Yes | 10 |
| 6 | Figure 4 enhancement | Yes | 10 |
| 7 | Figure 5 enhancement | Yes | 10 |
| 8 | Figure 6 enhancement | Yes | 10 |
| 9 | Cross-figure consistency | No | - |
| 10 | Final verification | No | - |

**Expected Duration:** 8-10 hours (with iteration time)
**Commits:** ~10 (1 per task)
