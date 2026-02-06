# CITEgeist Patterns/Cell Reports Methods Manuscript Rewrite

**Date**: 2026-02-05
**Target Journals**: Patterns, Cell Reports Methods (Cell Press multi-journal submission)
**Word Limit**: 5,000-7,000 words (research article)

---

## Core Framing

### Title Direction
"CITEgeist: A spatially-native framework for multi-modal integration and program discovery in spatial transcriptomics"

### Three Principles

1. **Spatially-native** - Every module uses spatial information (Moran's I, spatial colocalization, neighborhood optimization, spatial coherence). Contrasts with methods that treat coordinates as metadata.

2. **Integration-first** - Designed for same-slide proteomics+transcriptomics from the start. As technologies mature (Visium HD + antibodies), CITEgeist is ready.

3. **Resolution-agnostic** - Works at spot-level (Visium) and single-cell (Xenium/CosMx). Demonstrated on both.

### Key Differentiator
"Unlike reference-dependent methods that transfer information from external atlases, CITEgeist leverages the protein layer from the same tissue section to anchor all downstream discovery - no hallucination of spatial patterns from non-spatial data."

---

## Paper Structure (~6,000 words)

### 1. Introduction (~800 words)
- Spatial transcriptomics revolution + limitations of current methods
- The "single-cell methods applied to spatial data" problem
- Same-slide proteomics as the solution (CITE-seq, antibody capture)
- CITEgeist: spatially-native, integration-first, resolution-agnostic

### 2. Results

#### 2.1 Framework Overview (~500 words)
- Module architecture diagram (Figure 1)
- Data flow: raw data → markers → profiles → proportions → programs → integration

#### 2.2 Modules 1-2: Automated Profile Discovery (~600 words)
- Spatial marker detection (Moran's I, GMM)
- Colocalization-based profile assembly
- Xenium single-cell demonstration

#### 2.3 Module 3: Spatial Deconvolution (~500 words)
- Protein-anchored proportion estimation
- Neighborhood-aware refinement
- Benchmarking results (vs Cell2Location, RCTD, Tangram, Seurat)

#### 2.4 Module 4: Spatial Program Discovery (~700 words)
- Per-cell-type NMF with spatial coherence
- Xenium full demonstration (resolution-agnostic proof)
- Program examples with Moran's I validation

#### 2.5 Module 5: Cross-Sample Integration (~700 words)
- Harmony-based program alignment across 14 patients
- Responder vs progressor enriched programs **(TOUR DE FORCE)**
- PyDESeq2 differential expression on program-stratified spots
- Conserved relationships network

#### 2.6 Interpretable Outputs Enable Downstream Analysis (~500 words)
- Midkine discovery via cell-cell signaling tools (COMMOT)
- Compatibility with external methods (PyDESeq2, GSEApy, Harmony)
- Wet lab validation summary

### 3. Discussion (~700 words)
- Spatial-native vs spatial-aware distinction
- Future: Visium HD, emerging technologies
- Limitations and future work

---

## Figure Strategy

### Main Figures (6)

**Figure 1: Pipeline Overview (Full page)**
- Panel A: Schematic of 5 modules with data flow
- Panel B: Spatial-native operations highlighted (Moran's I, colocalization graph, neighborhood optimization)
- Panel C: Resolution flexibility diagram (Visium spot → single cell)

**Figure 2: Modules 1-2 Profile Discovery**
- Panel A: Marker interest detection (kurtosis, GMM, Moran's I gates)
- Panel B: Colocalization network → hierarchical clustering → profiles
- Panel C: Xenium single-cell demonstration (same workflow, cell resolution)
- Panel D: Discovered profiles vs known markers (validation)

**Figure 3: Module 3 Deconvolution Benchmarking**
- Panel A: Proportion estimation schematic
- Panel B: Benchmark comparison (JSD, RMSE, correlation) - CITEgeist vs Cell2Location, RCTD, Tangram, Seurat
- Panel C: Spatial visualization of proportions (example region)

**Figure 4: Module 4 Spatial Programs (TOUR DE FORCE #1)**
- Panel A: NMF program discovery schematic
- Panel B: Example programs per cell type with top genes
- Panel C: Moran's I validation (programs are spatially coherent)
- Panel D: Xenium single-cell programs (resolution-agnostic proof)
- Panel E: Bivariate relationships (co-localized vs exclusive programs)

**Figure 5: Module 5 Cross-Sample Integration (TOUR DE FORCE #2)**
- Panel A: 14-sample integration schematic
- Panel B: UMAP of aligned programs colored by sample/cell type
- Panel C: Responder vs Progressor enriched programs
- Panel D: PyDESeq2 differential expression results
- Panel E: Conserved relationship network

**Figure 6: Interpretability & External Tool Compatibility**
- Panel A: CITEgeist outputs → cell-cell signaling analysis workflow
- Panel B: Midkine discovery pathway (condensed)
- Panel C: Wet lab validation summary

### Supplementary Figures
- Extended midkine mechanism (ER saturation model, full wet lab)
- Additional benchmarking metrics
- Parameter sensitivity analyses
- Per-sample Module 4 results
- Xenium Module 1-2 detailed results

---

## Key Claims & Evidence

### Claim 1: "CITEgeist is spatially-native, not just spatially-aware"
**Evidence:**
- Module 1: Moran's I for marker detection
- Module 2: Spatial colocalization (adjacent-spot enrichment, bivariate Moran's I)
- Module 3: Laplacian smoothing, neighborhood finetuning
- Module 4: Moran's I validation of programs
- Module 5: Spatial relationship conservation
- Comparison table contrasting with methods that just use coordinates as features

### Claim 2: "Same-slide proteomics eliminates reference dependency and hallucination"
**Evidence:**
- No external scRNA-seq atlas required
- Profiles derived from same tissue section
- Benchmarks show competitive/superior accuracy without reference

### Claim 3: "Resolution-agnostic - works spot-level to single-cell"
**Evidence:**
- Visium patient cohort (14 samples, Modules 1-5)
- Xenium single-cell (Modules 1-2-4 full demonstration)
- Discussion: Architectural readiness for Visium HD + proteomics

### Claim 4: "Discovers clinically relevant biology"
**Evidence:**
- Responder vs progressor programs (n=14, Module 5)
- PyDESeq2 DE on program-stratified spots
- Spatially coherent programs per cell type (Module 4)
- Midkine pathway discovered and validated (supplement)

### Claim 5: "Outputs are interpretable and interoperable"
**Evidence:**
- Standard scanpy/anndata workflows (Leiden, rank_genes_groups, Harmony)
- PyDESeq2 differential expression (Vignette 3, Module 5)
- GSEApy/Enrichr pathway analysis
- COMMOT cell-cell communication
- Human Protein Atlas validation
- Bulk RNA-seq interaction effects
- ChIP-seq binding validation
- Published gene signatures

---

## Pre-Submission Work

### Computational (To Run)

1. **Xenium Module 4 Full Demonstration**
   - Run Module 4 on Xenium single-cell data
   - Multiple cell type anchors
   - Report programs, Moran's I, top genes
   - **Status**: NOT DONE

2. **Module 5 PyDESeq2 Analysis**
   - Stratify spots by responder vs progressor program enrichment
   - Run PyDESeq2 differential expression
   - Pathway analysis on DE genes
   - **Status**: NOT DONE

3. **Figure Generation**
   - Pipeline overview (Figure 1)
   - Module 1-2 discovery (Figure 2)
   - Benchmarking comparison (Figure 3)
   - Module 4 programs (Figure 4)
   - Module 5 integration (Figure 5)
   - Interoperability/midkine (Figure 6)
   - **Status**: PARTIAL - some exist, need consolidation

### Writing

4. **Manuscript Rewrite**
   - New intro emphasizing spatial-native + integration-first
   - Results sections per design above
   - Condensed midkine to interoperability demonstration
   - Discussion on future technologies
   - **Status**: NOT DONE - major rewrite from current draft

---

## Midkine Story Handling

**Current status**: Extensive wet lab validation exists, but only 1/6 patients has ESR1 mutation, limiting generalizability.

**Decision**: Move to supplement as interoperability demonstration.

**Main text role**:
- Brief mention in Section 2.6 as example of what pipeline enables
- Emphasize COMMOT compatibility and interpretable outputs
- Condensed to ~1 figure panel

**Supplement**:
- Full ER saturation mechanism
- MCF7 vs T47D differential response explanation
- Complete wet lab validation data

---

## Existing Resources

### Vignette Demonstrations (can adapt for paper)
- **Vignette 1**: COMMOT, GSEApy, HPA validation
- **Vignette 2**: D538G mapping, mutation-stratified signaling, GSEA
- **Vignette 3**: PyDESeq2 on macrophage-filtered spots, Harmony integration
- **Vignette 4**: Module 4 region analysis, contextual gene extraction

### Module 4/5 Outputs Available
- `output/module4_discovery/` - 7 patient samples with programs
- `output/module5_integration/` - 14-sample integration, 73 aligned programs, 191 conserved relationships
- Response analysis: 3 responder-enriched, 4 progressor-enriched programs

### Benchmarking Data
- `Benchmarking/xenium_benchmarking/` - Xenium pseudo-Visium with 10 granular cell types
- Comparison methods: Cell2Location, RCTD, Tangram, Seurat

---

## Timeline Considerations

1. Run Xenium Module 4 (computational)
2. Run Module 5 PyDESeq2 (computational)
3. Generate/consolidate figures
4. Write manuscript draft
5. Internal review
6. Submit to Patterns (or Cell Reports Methods via multi-journal)

---

## Notes

- Cell Press multi-journal submission allows targeting both Patterns and Cell Reports Methods
- Patterns emphasizes data science utility; Cell Reports Methods emphasizes biological methods
- Both have similar formatting requirements
- No strict figure limit for either journal
