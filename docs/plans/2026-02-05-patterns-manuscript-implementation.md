# CITEgeist Patterns Manuscript - Implementation Plan

**Design Document**: `docs/plans/2026-02-05-patterns-manuscript-rewrite-design.md`
**Date**: 2026-02-05

---

## Phase 1: Computational Analyses

### Task 1.1: Xenium Module 4 Full Demonstration
**Priority**: HIGH (blocks Figure 4)
**Status**: ✅ COMPLETE

**Existing outputs** in `Benchmarking/xenium_benchmarking/CITEgeist/output_module4_validation/singlecell/`:
- 5 regions processed
- 7 cell types: CD8+ T, CD4+ T, Fibroblasts, Macrophages, Endothelial, B cells, Epithelial
- 5 programs per cell type (175 total)
- Moran's I validation complete (range: 0.04-0.37, programs with I>0.2 indicate spatial coherence)
- Module 5 integration: 27 aligned programs, 207 conserved relationships

**Key results**:
| Cell Type | Best Moran's I | Example Genes |
|-----------|----------------|---------------|
| Fibroblasts | 0.367 | HLA-DRA, FCGR3A, HLA-DRB1 |
| Endothelial | 0.277 | HLA-DRA, HLA-DRB1, PTPRC |
| CD4+ T | 0.263 | VIM, CD3E, PTEN, ACTA2 |
| CD8+ T | 0.260 | VIM, CXCL6, PTEN, VCAN |

**No further computation needed** - ready for Figure 4 generation.

---

### Task 1.2: Module 5 PyDESeq2 Analysis
**Priority**: HIGH (blocks Figure 5)
**Status**: ✅ COMPLETE
**Executed**: 2026-02-06

**Approach** (pseudo-bulk aggregation):
1. Aggregate deconvolved gene expression within each of 14 samples (pseudo-bulk)
2. Compare responder samples (P1, P5: 5 samples) vs progressor samples (P2-P4, P6: 9 samples)
3. Run PyDESeq2 differential expression with 8 CPUs
4. Pathway enrichment with GSEApy (KEGG, GO, MSigDB Hallmark)

**Key Results**:
| Metric | Value |
|--------|-------|
| Genes tested | 13,371 |
| Significant (padj<0.05) | 127 |
| Responder-enriched | 5 |
| Progressor-enriched | 122 |

**Top Progressor-Enriched Genes** (negative log2FC = higher in progressors):
- SLC37A2, MMP13, BNC2, ARRDC4, ACAA2, ADAMTS4, MLKL

**Pathway Enrichment**:
- 6 GO BP pathways for responder-up genes
- 7 GO BP pathways for progressor-up genes
- 1 MSigDB Hallmark pathway for progressor-up genes

**Output** (`examples/output_module5_pydeseq/`):
- `pseudobulk_de_results.csv` - All DE genes
- `pseudobulk_*_enrichment.csv` - Pathway enrichment files
- `module5_pydeseq_summary.json` - Summary statistics

**Scripts**:
- `examples/run_module5_pydeseq_pseudobulk.py` ✅
- `examples/sbatch_module5_pydeseq_pseudobulk.sh` ✅

---

## Phase 2: Figure Generation

### Task 2.1: Figure 1 - Pipeline Overview
**Priority**: MEDIUM
**Status**: ✅ COMPLETE
**Dependencies**: None

**Panels**:
- A: Module flow schematic
- B: Spatial operations highlight
- C: Resolution flexibility diagram

**Script**: `manuscript/figures/generate_figure1.py` ✅
**Output**: `manuscript/figures/output/figure1_pipeline_overview.pdf`

---

### Task 2.2: Figure 2 - Modules 1-2 Profile Discovery
**Priority**: MEDIUM
**Status**: ✅ COMPLETE
**Dependencies**: Xenium Module 1-2 outputs (already exist)

**Panels**:
- A: Marker interest detection examples
- B: Colocalization network → dendrogram → profiles
- C: Xenium single-cell workflow
- D: Discovered vs known markers comparison

**Script**: `manuscript/figures/generate_figure2.py` ✅
**Output**: `manuscript/figures/output/figure2_profile_discovery.pdf`

---

### Task 2.3: Figure 3 - Benchmarking
**Priority**: MEDIUM
**Status**: ✅ COMPLETE
**Dependencies**: Existing benchmark results

**Panels**:
- A: Deconvolution schematic
- B: Method comparison (bar plots: JSD, RMSE, correlation)
- C: Spatial proportion visualization

**Data source**: `Benchmarking/xenium_benchmarking/evaluation/results/`
**Script**: `manuscript/figures/generate_figure3.py` ✅
**Output**: `manuscript/figures/output/figure3_benchmarking.pdf`

---

### Task 2.4: Figure 4 - Module 4 Programs
**Priority**: HIGH
**Status**: ✅ COMPLETE
**Dependencies**: Task 1.1 (Xenium Module 4)

**Panels**:
- A: NMF schematic
- B: Program examples (heatmap of top genes)
- C: Moran's I validation (scatter: program activity vs spatial coherence)
- D: Xenium single-cell programs (spatial plots)
- E: Bivariate relationships (heatmap)

**Script**: `manuscript/figures/generate_figure4.py` ✅
**Output**: `manuscript/figures/output/figure4_module4_programs.pdf`

---

### Task 2.5: Figure 5 - Module 5 Integration
**Priority**: HIGH
**Status**: ✅ COMPLETE
**Dependencies**: Task 1.2 (PyDESeq2)

**Panels**:
- A: Integration schematic
- B: UMAP of aligned programs (590 programs colored by cell type)
- C: Responder vs progressor bar plot (3 resp-enriched, 4 prog-enriched)
- D: PyDESeq2 volcano plot (127 sig genes, 5 resp-up, 122 prog-up)
- E: Conserved relationship network (191 relationships)

**Script**: `manuscript/figures/generate_figure5.py` ✅
**Output**: `manuscript/figures/output/figure5_integration.pdf`

---

### Task 2.6: Figure 6 - Interoperability
**Priority**: MEDIUM
**Status**: ✅ COMPLETE
**Dependencies**: Existing vignette outputs

**Panels**:
- A: Workflow diagram (CITEgeist → external tools)
- B: Midkine discovery summary
- C: Wet lab validation (condensed)

**Script**: `manuscript/figures/generate_figure6.py` ✅
**Output**: `manuscript/figures/output/figure6_interoperability.pdf`

---

## Phase 3: Manuscript Writing

### Task 3.1: Introduction Rewrite
**Priority**: HIGH
**Dependencies**: None

**Key points**:
- Spatial transcriptomics limitations
- Single-cell methods applied to spatial data problem
- Same-slide proteomics solution
- CITEgeist three principles

---

### Task 3.2: Results Sections
**Priority**: HIGH
**Dependencies**: Figures

**Sections**:
- 2.1 Framework Overview
- 2.2 Modules 1-2
- 2.3 Module 3 + Benchmarking
- 2.4 Module 4
- 2.5 Module 5
- 2.6 Interoperability

---

### Task 3.3: Discussion
**Priority**: MEDIUM
**Dependencies**: Results

**Key points**:
- Spatial-native distinction
- Future technologies
- Limitations

---

### Task 3.4: Supplement
**Priority**: LOW
**Dependencies**: Main text complete

**Content**:
- Extended midkine mechanism
- Additional benchmarks
- Parameter sensitivity
- Per-sample details

---

## Execution Order

```
Week 1: Computational
├── Day 1-2: Task 1.1 (Xenium Module 4)
└── Day 2-3: Task 1.2 (PyDESeq2)

Week 2: Figures
├── Day 1: Tasks 2.1, 2.2, 2.3 (pipeline, modules 1-2, benchmarks)
└── Day 2-3: Tasks 2.4, 2.5, 2.6 (modules 4-5, interop)

Week 3: Writing
├── Day 1-2: Task 3.1 (Introduction)
├── Day 3-4: Task 3.2 (Results)
└── Day 5: Task 3.3 (Discussion)

Week 4: Polish
├── Day 1-2: Task 3.4 (Supplement)
├── Day 3: Internal review
└── Day 4-5: Revisions
```

---

## File Structure

```
manuscript/
├── figures/
│   ├── generate_figure1.py
│   ├── generate_figure2.py
│   ├── generate_figure3.py
│   ├── generate_figure4.py
│   ├── generate_figure5.py
│   ├── generate_figure6.py
│   └── output/
│       ├── figure1_pipeline.pdf
│       ├── figure2_modules12.pdf
│       ├── ...
├── CITEgeistManuscript_v4_patterns.docx  (new draft)
└── supplement/
    ├── extended_midkine.docx
    └── supplementary_figures/
```

---

## Checkpoints

- [x] Task 1.1: Xenium Module 4 complete ✅ (outputs exist in output_module4_validation/singlecell/)
- [x] Task 1.2: PyDESeq2 analysis complete ✅ (127 DE genes, pseudo-bulk approach)
- [x] Task 2.1-2.6: All figures generated ✅ (2026-02-06)
      - Figure 1: Pipeline Overview
      - Figure 2: Profile Discovery (Modules 1-2)
      - Figure 3: Benchmarking
      - Figure 4: Spatial Programs (Module 4)
      - Figure 5: Cross-Sample Integration (Module 5 + PyDESeq2 volcano)
      - Figure 6: Interoperability
- [ ] Task 3.1: Introduction draft
- [ ] Task 3.2: Results draft
- [ ] Task 3.3: Discussion draft
- [ ] Internal review
- [ ] Final submission draft
