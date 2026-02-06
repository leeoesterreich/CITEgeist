# CITEgeist Patterns Manuscript - Implementation Plan

**Design Document**: `docs/plans/2026-02-05-patterns-manuscript-rewrite-design.md`
**Date**: 2026-02-05

---

## Phase 1: Computational Analyses

### Task 1.1: Xenium Module 4 Full Demonstration
**Priority**: HIGH (blocks Figure 4)
**Estimated compute**: ~2-4 hours on GPU node

**Steps**:
1. Load Xenium single-cell data with Module 1-2 outputs
2. Run Module 4 `discover_programs_from_layers()` for each cell type anchor
3. Compute Moran's I for spatial coherence validation
4. Run Module 4b bivariate relationships
5. Save outputs to `Benchmarking/xenium_benchmarking/CITEgeist/output_module4_singlecell/`

**Input**:
- Xenium data: `Benchmarking/xenium_pseudovisium/data/`
- Module 1-2 outputs: Already run on Xenium

**Output**:
- Program W matrices (gene loadings)
- Program H matrices (spot/cell loadings)
- Moran's I statistics per program
- Top genes per program
- Bivariate relationship matrix

**Script location**: `Benchmarking/xenium_benchmarking/CITEgeist/src/run_module4_xenium.py`
**SLURM script**: `Benchmarking/xenium_benchmarking/CITEgeist/slurm/run_module4_xenium.sh`

---

### Task 1.2: Module 5 PyDESeq2 Analysis
**Priority**: HIGH (blocks Figure 5)
**Estimated compute**: ~1-2 hours

**Steps**:
1. Load Module 5 outputs (`output/module5_integration/`)
2. Identify responder-enriched vs progressor-enriched programs
3. For each patient sample:
   - Stratify spots by program enrichment (top 25% vs bottom 25%)
   - Extract raw counts from deconvolved layers
4. Run PyDESeq2 comparing responder-high vs progressor-high spots
5. Run pathway analysis (GSEApy) on DE genes
6. Save results

**Input**:
- Module 5: `output/module5_integration/module5_response_analysis.json`
- Patient data: Deconvolved layers from Module 3

**Output**:
- DE gene lists (responder vs progressor program spots)
- Pathway enrichment results
- Summary statistics

**Script location**: `examples/run_module5_pydeseq.py`
**SLURM script**: `examples/sbatch_module5_pydeseq.sh`

---

## Phase 2: Figure Generation

### Task 2.1: Figure 1 - Pipeline Overview
**Priority**: MEDIUM
**Dependencies**: None

**Panels**:
- A: Module flow schematic (draw in Python/matplotlib or external tool)
- B: Spatial operations highlight
- C: Resolution flexibility diagram

**Script**: `manuscript/figures/generate_figure1.py`

---

### Task 2.2: Figure 2 - Modules 1-2 Profile Discovery
**Priority**: MEDIUM
**Dependencies**: Xenium Module 1-2 outputs (already exist)

**Panels**:
- A: Marker interest detection examples
- B: Colocalization network → dendrogram → profiles
- C: Xenium single-cell workflow
- D: Discovered vs known markers comparison

**Script**: `manuscript/figures/generate_figure2.py`

---

### Task 2.3: Figure 3 - Benchmarking
**Priority**: MEDIUM
**Dependencies**: Existing benchmark results

**Panels**:
- A: Deconvolution schematic
- B: Method comparison (bar plots: JSD, RMSE, correlation)
- C: Spatial proportion visualization

**Data source**: `Benchmarking/xenium_benchmarking/evaluation/results/`
**Script**: `manuscript/figures/generate_figure3.py`

---

### Task 2.4: Figure 4 - Module 4 Programs
**Priority**: HIGH
**Dependencies**: Task 1.1 (Xenium Module 4)

**Panels**:
- A: NMF schematic
- B: Program examples (heatmap of top genes)
- C: Moran's I validation (scatter: program activity vs spatial coherence)
- D: Xenium single-cell programs (spatial plots)
- E: Bivariate relationships (heatmap)

**Script**: `manuscript/figures/generate_figure4.py`

---

### Task 2.5: Figure 5 - Module 5 Integration
**Priority**: HIGH
**Dependencies**: Task 1.2 (PyDESeq2)

**Panels**:
- A: Integration schematic
- B: UMAP of aligned programs
- C: Responder vs progressor bar plot
- D: PyDESeq2 volcano plot / pathway results
- E: Conserved relationship network

**Script**: `manuscript/figures/generate_figure5.py`

---

### Task 2.6: Figure 6 - Interoperability
**Priority**: MEDIUM
**Dependencies**: Existing vignette outputs

**Panels**:
- A: Workflow diagram (CITEgeist → external tools)
- B: Midkine discovery summary
- C: Wet lab validation (condensed)

**Script**: `manuscript/figures/generate_figure6.py`

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

- [ ] Task 1.1: Xenium Module 4 complete
- [ ] Task 1.2: PyDESeq2 analysis complete
- [ ] Task 2.1-2.6: All figures generated
- [ ] Task 3.1: Introduction draft
- [ ] Task 3.2: Results draft
- [ ] Task 3.3: Discussion draft
- [ ] Internal review
- [ ] Final submission draft
