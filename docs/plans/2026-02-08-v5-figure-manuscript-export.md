# V5 Manuscript Figure Alignment & Export Plan

**Date**: 2026-02-08
**Goal**: Align figures with v5 manuscript text, export panels for Illustrator, create Word docs

## Context

v5 manuscript (`CITEgeist_CRM_v5.md`) restructured results significantly from v4:
- Section 2.4 = midkine/ESR1 case study (NEW Figure 4, 8 panels A-H)
- Section 2.5 = full pipeline application with Module 4+5 results (NEW Figure 5)
- Figure 6 (interoperability) retired - message woven into text
- ESR1 D538G patient = HCC22-088-P4-S2

Reference formatting: `CITEgeistFigures.docx` (unversioned, manually polished)

## Deliverables

| # | Deliverable | Format |
|---|------------|--------|
| 1 | v5 manuscript as Word doc | `manuscript/CITEgeist_CRM_v5.docx` |
| 2 | v5 figure legends + embedded images as Word doc | `manuscript/CITEgeist_Patterns_v5_Figures.docx` |
| 3 | v5 figure legends markdown | `manuscript/CITEgeist_Patterns_v5_Figures.md` |
| 4 | Rewritten Figure 4 script (midkine) | `manuscript/figures/generate_figure4.py` |
| 5 | Restructured Figure 5 script (Mod4+5) | `manuscript/figures/generate_figure5.py` |
| 6 | Individual panel exports (PNG 300DPI + SVG) | `manuscript/figures/export/Figure{1-5}/` |
| 7 | Text-figure cross-reference audit | `docs/review_2026-02-07/v5_figure_audit.md` |

## Figure Content Map (v5)

### Figure 1: Pipeline Overview (UNCHANGED)
- 1A: Modular pipeline schematic
- 1B: Spatial statistics at each stage
- 1C: Resolution flexibility

### Figure 2: Profile Discovery (UNCHANGED)
- 2A: Marker interest detection
- 2B: Profile discovery workflow
- 2C: Xenium single-cell demo
- 2D: Profile validation table

### Figure 3: Benchmarking (UNCHANGED)
- 3A: Two-pass deconvolution schematic
- 3B: Proportion benchmarking bars
- 3C: GEX benchmarking by cell type
- 3D: Predicted vs ground truth scatter

### Figure 4: Midkine/ESR1 Case Study (REWRITE)
- 4A: Basal cytokeratin signatures in cancer layer (spatial plot, P4-S2)
- 4B: Spatial distribution of mutation signal (not clustered)
- 4C: EstroGene mutation signatures
- 4D: MSigDB Hallmark pathways upregulated in D538G spots
- 4E: KEGG pathways upregulated in D538G spots
- 4F: COMMOT signaling bar plot (MDK, PTN, MIF)
- 4G: ELISA validation (PLACEHOLDER - data in Prism)
- 4H: IF validation (PLACEHOLDER - images manual)

### Figure 5: Full Pipeline - Module 4+5 (RESTRUCTURE)
- 5A: Module 4 NMF program discovery + Moran's I validation
- 5B: UMAP of 590 programs colored by cell type
- 5C: Response-associated programs bar chart
- 5D: PyDESeq2 volcano plot

## Execution Phases

### Phase 1: Parallel script work + text conversion
- [ ] Convert v5 md to docx (pandoc)
- [ ] Rewrite generate_figure4.py for midkine
- [ ] Restructure generate_figure5.py for Module 4+5
- [ ] Update figure legends markdown to v5

### Phase 2: Generate figures (sbatch)
- [ ] Run updated figure scripts
- [ ] Export individual panels to per-figure folders

### Phase 3: Assembly + review
- [ ] Create figures Word doc with embedded images
- [ ] Deploy review agents for cross-checking

## Review Agent Assignments

### Agent R1: Text-to-Figure Coherence
- For every `(Figure X)` reference in v5 text, verify the panel exists and content matches
- Check numbers, gene names, statistics in legends vs text body

### Agent R2: Figure-to-Text Coherence
- For every panel in each figure, verify it's referenced in the text
- Check no orphaned panels or missing references
- Verify panel letter sequence (A, B, C...) is continuous

### Agent R3: v3-to-v5 Comparison
- Compare old unversioned figures doc against v5 figure plan
- Identify what changed, what's missing, what's new
- Flag any content from v3 figures that should be preserved in v5

### Agent R4: Internal Consistency
- Cross-check statistics between text, legends, and figure data
- Verify sample IDs, patient counts, gene names are consistent
- Check that responder/progressor classification is consistent throughout
