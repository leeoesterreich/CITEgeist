# Figure Enhancement Design - Cell/Cell Reports Quality

**Date:** 2026-02-06
**Status:** Approved for overnight execution

---

## Problem Statement

Current manuscript figures (1-6) are underwhelming for a Nature/Cell-level publication:
- Too sparse - not enough data packed per panel
- Gaudy colors that "scream AI-generated"
- Text overlap and layout issues
- Missing statistical annotations in places
- V3 content (simulated benchmarks, vignettes, wet lab validation) not properly migrated

## Goal

Elevate all 6 main figures to Cell/Cell Reports publication quality through systematic review and iterative refinement.

## Success Criteria (Per Panel)

| Criterion | Check |
|-----------|-------|
| **Visual** | No text overlap, readable at 89mm column width, proper spacing |
| **Legibility** | Minimum 10pt font, adequate contrast, clear labels |
| **Information Density** | ≥2 data dimensions shown, no wasted whitespace |
| **Statistical Rigor** | Comparisons show p-values/CIs where applicable |
| **Style** | Cell/Cell Reports aesthetic - professional, saturated but not gaudy |
| **Color-blind Safe** | Palette passes accessibility check |

## Style Guide

**Cell/Cell Reports Reference Style:**

| Element | Specification |
|---------|---------------|
| **Font** | Arial/Helvetica, 9-10pt labels, 11pt panel letters |
| **Line weight** | 0.5-1pt axes, 1-1.5pt data lines |
| **Panel letters** | Bold, uppercase, top-left outside panel |
| **Legends** | Inside panel when space allows, no border |
| **Whitespace** | Aggressive minimization - tight margins |

**Color Palette:**
```python
PALETTE = {
    'primary': '#2171b5',    # Deep blue
    'secondary': '#6baed6',  # Medium blue
    'accent1': '#fe9929',    # Warm orange
    'accent2': '#41ab5d',    # Forest green
    'neutral': '#636363',    # Gray
    'highlight': '#c51b8a',  # Magenta (sparingly)
}
```

**Matplotlib Defaults:**
```python
plt.rcParams.update({
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
})
# Whitespace: tight_layout(pad=0.3), wspace=0.15, hspace=0.2
```

## Content Sources

**Current Figure Scripts:**
- `manuscript/figures/generate_figure[1-6].py`

**V3 Content to Migrate:**
- `manuscript/CITEgeistManuscript_v3.docx` - reference for dropped content
- `Benchmarking/simulation_benchmarking/` - simulated benchmark results
- `Benchmarking/xenium_benchmarking/evaluation/results/` - Xenium benchmark data

**Vignette Outputs:**
- `CITEgeist/examples/vignette_*.ipynb`
- `examples/output_module5_pydeseq/`

**Wet Lab Validation:**
- `midkine/CITEgeistImageAnalysis.prism`
- `midkine/ELISA Results.prism`
- `midkine/mdk_saturation_pipeline/outputs/`

## Architecture: Three-Layer Execution

```
┌─────────────────────────────────────────────────────────┐
│  hpc-isolate:overnight                                  │
│  (Isolated Claude on HPC, autonomous execution)         │
│                                                         │
│  ┌───────────────────────────────────────────────────┐  │
│  │  superpowers:executing-plans                      │  │
│  │  (Batch execution, commits after each task)       │  │
│  │                                                   │  │
│  │  For each figure (Tasks 2-7):                     │  │
│  │  ┌─────────────────────────────────────────────┐  │  │
│  │  │  ralph-loop                                 │  │  │
│  │  │  --completion-promise "Fig N passes all"    │  │  │
│  │  │  --max-iterations 10                        │  │  │
│  │  │                                             │  │  │
│  │  │  Iterate: check → fix → regen → re-check   │  │  │
│  │  │  Until: all panels genuinely pass          │  │  │
│  │  └─────────────────────────────────────────────┘  │  │
│  │  ► Commit after each figure                       │  │
│  └───────────────────────────────────────────────────┘  │
│                                                         │
│  ► OVERNIGHT_SUMMARY.md + notification                  │
└─────────────────────────────────────────────────────────┘
```

## Constraints

- **Don't delete existing content** - additive improvements only
- **Can run code** - regenerate figures as needed
- **Avoid regenerating large datasets** - use existing benchmark outputs
- **New scripts OK** - especially for vignette integration
- **Full autonomy** - no human interaction expected

## Expected Outputs

1. Enhanced figure PDFs in `manuscript/figures/output/`
2. Updated figure scripts with style improvements
3. `FIGURE_AUDIT.md` documenting gaps found
4. Commits per figure with descriptive messages
5. `OVERNIGHT_SUMMARY.md` with morning review instructions

---

**Implementation Plan:** `docs/plans/2026-02-06-figure-enhancement-plan.md`
