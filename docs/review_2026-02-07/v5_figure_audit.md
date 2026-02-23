# V5 Manuscript vs V4 Figure Legends: Cross-Reference Audit

**Date**: 2026-02-08
**Manuscript**: `manuscript/CITEgeist_CRM_v5.md`
**Figure Legends**: `manuscript/CITEgeist_Patterns_v4_Figures.md`

---

## A. Text to Figure References (Every Reference in V5 Manuscript)

### Figure 1 References

| # | Manuscript Line | Text Context | What Text Says Panel Contains | Legend Match? | Status |
|---|---|---|---|---|---|
| 1 | L63 | "...passes interpretable outputs to subsequent modules **(Figure 1A)**" | Schematic of five-module pipeline with data flow | Legend 1A: "Schematic of the five-module pipeline...Arrows indicate data flow between modules" | MATCH |
| 2 | L67 | "...interoperable with tools such as PyDESeq2, GSEApy, and COMMOT **(Figure 1B-C)**" | Interoperability with downstream tools + output formats | Legend 1B: Spatial-native operations at each stage; Legend 1C: Resolution flexibility across platforms | **WARNING** |

**Issue #1 (WARNING)**: The v5 text references Figure 1B-C in the context of tool interoperability and standard output formats. However, the v4 legends describe:
- **1B** = Spatial-native operations (Moran's I, colocalization networks, Laplacian regularization)
- **1C** = Resolution flexibility (Visium vs Xenium)

Neither 1B nor 1C is about interoperability with PyDESeq2/GSEApy/COMMOT. The interoperability content actually appears in **Figure 6A** of the v4 legends. Either the text description needs to match what 1B-C actually shows, or 1B-C needs to be redesigned to include interoperability content.

---

### Figure 2 References

| # | Manuscript Line | Text Context | What Text Says Panel Contains | Legend Match? | Status |
|---|---|---|---|---|---|
| 3 | L71 | "...automated discovery pipeline...biologically coherent cell type profiles **(Figure 2A-B)**" | Module 1 marker detection + Module 2 profile assembly | Legend 2A: Module 1 three statistical gates; Legend 2B: Module 2 colocalization network + clustering | MATCH |
| 4 | L77 | "...applied Modules 1-2 to Xenium single-cell data where ground truth cell types were available **(Figure 2C)**" | Xenium single-cell validation of Modules 1-2 | Legend 2C: "Xenium single-cell demonstration...discovers biologically coherent profiles without reference data" | MATCH |
| 5 | L77 | "...spatial colocalization provides sufficient signal to recover biologically meaningful cell type definitions **(Figure 2D)**" | Marker clustering validation (CD3E/CD8A, CD68/CD163, EPCAM/KRT) | Legend 2D: "discovered profiles correctly recover known cell type markers. CD3E and CD8A cluster together (T cells)..." | MATCH |

---

### Figure 3 References

| # | Manuscript Line | Text Context | What Text Says Panel Contains | Legend Match? | Status |
|---|---|---|---|---|---|
| 6 | L81 | "...estimates cell type proportions...deconvolves bulk gene expression **(Figure 3A)**" | Two-pass deconvolution schematic | Legend 3A: "Two-pass deconvolution schematic...Pass 1...Pass 2..." | MATCH |
| 7 | L83-85 | "...Tangram...reporting average JSD of 0.56... **(Figure 3B)**" | Simulated and Xenium benchmarking metrics across methods | Legend 3B: "Benchmarking on Xenium pseudo-Visium spots...Pearson correlation...JSD...RMSE" | **WARNING** |

**Issue #2 (WARNING)**: The v5 text paragraph referencing Figure 3B covers BOTH simulated benchmarking (scCube, mixed vs segmented) AND Xenium pseudo-Visium benchmarking. The v4 Figure 3B legend only describes Xenium pseudo-Visium benchmarking. The simulated benchmarking content (CITEgeist r=0.95, JSD 0.16; Cell2Location RMSE degradation from 0.08 to 0.167; RCTD reference dependency 0.05 to 0.21; Tangram JSD 0.56) has no corresponding panel in the figure legends. Either:
- A separate panel should be added for simulated benchmarking, or
- Panel 3B legend should be expanded to cover both simulated and real benchmarks, or
- The simulated results live in a supplementary figure not yet defined.

**Note**: The v5 text also describes gene expression deconvolution benchmarking (NRMSE 0.04 vs 0.16, scResolve r=0.43 comparison) without a specific figure panel reference. This content likely corresponds to Figure 3C (GEX benchmarking) but is not explicitly cited with a panel letter in the text.

---

### Figure 4 References

| # | Manuscript Line | Text Context | What Text Says Panel Contains | Legend Match? | Status |
|---|---|---|---|---|---|
| 8 | L93 | "...confirmed the presence of basal gene signatures in the cancer cell layer **(Figure 4A)**" | ESR1-mutant basal cytokeratin spatial expression | Legend 4A (v4): NMF-based program discovery schematic | **CRITICAL** |
| 9 | L93 | "...signals were not clustered in a single location **(Figure 4B)**" | Spatial distribution of basal signals (dispersed, not clonal) | Legend 4B (v4): Program examples across cell types with Moran's I | **CRITICAL** |
| 10 | L93 | "...expected gene mutation signatures from EstroGene **(Figure 4C)**" | EstroGene mutation signature spatial map | Legend 4C (v4): Moran's I validation box plot by cell type | **CRITICAL** |
| 11 | L93 | "...known pathway alterations from prior work **(Figure 4D-E)**" | Known ESR1-mutant pathway analysis spatial maps | Legend 4D (v4): Summary statistics table, bivariate Moran's I | **CRITICAL** |
| 12 | L95 | "...COMMOT analysis then revealed increased MDK, PTN, and MIF signaling **(Figure 4F)**" | COMMOT signaling analysis results | No panel 4F in v4 legends | **CRITICAL** |
| 13 | L97 | "...mutant cells secreted approximately double the midkine...ELISA **(Figure 4G)**" | ELISA bar chart for MCF7 WT vs D538G MDK | No panel 4G in v4 legends | **CRITICAL** |
| 14 | L97 | "...double the midkine at the cell membrane...immunofluorescence **(Figure 4H)**" | IF quantification, membrane + intracellular MDK | No panel 4H in v4 legends | **CRITICAL** |

**Issue #3 (CRITICAL)**: Figure 4 in v5 text describes the ESR1-mutant midkine case study (8 panels, A through H), but the v4 figure legends describe Figure 4 as "Module 4 - Spatial Program Discovery" (4 panels, A through D). There is a complete mismatch:

| Panel | V5 Text Says | V4 Legend Says |
|-------|-------------|----------------|
| 4A | Basal gene signatures in cancer layer | NMF program discovery schematic |
| 4B | Spatial distribution of basal signals | Program examples with Moran's I |
| 4C | EstroGene mutation signatures | Moran's I validation box plot |
| 4D | Pathway alterations | Summary statistics table |
| 4E | Pathway alterations (cont.) | *(does not exist)* |
| 4F | COMMOT MDK/PTN/MIF signaling | *(does not exist)* |
| 4G | ELISA results | *(does not exist)* |
| 4H | Immunofluorescence results | *(does not exist)* |

This indicates the v4 figure legends were written for a different figure numbering scheme where Module 4 results are in Figure 4 and the midkine case study content was elsewhere. The v5 text reorganized figures so that the midkine/ESR1 case study IS Figure 4, and the Module 4+5 pipeline results are in Figure 5.

---

### Figure 5 References

| # | Manuscript Line | Text Context | What Text Says Panel Contains | Legend Match? | Status |
|---|---|---|---|---|---|
| 15 | L109 | "...discovering gene expression programs within each cell type population **(Figure 5A)**" | Module 4 NMF program discovery + Moran's I validation | Legend 5A (v4): Integration schematic (PCA, Harmony, clustering) | **CRITICAL** |
| 16 | L111 | "...73 aligned programs representing conserved transcriptional states from 590 individual programs **(Figure 5B)**" | UMAP of 590 programs, 73 aligned clusters | Legend 5B (v4): UMAP of 590 programs, 73 aligned | PARTIAL MATCH |
| 17 | L113 | "...Response analysis compared program enrichment **(Figure 5C)**" | Response-associated programs bar plot (3 responder, 4 progressor) | Legend 5C (v4): Response-associated programs bar plot | MATCH |
| 18 | L113 | "...127 significantly differentially expressed genes...volcano plot **(Figure 5D)**" | PyDESeq2 volcano plot, 127 DE genes | Legend 5D (v4): PyDESeq2 volcano plot, 127 significant | MATCH |

**Issue #4 (CRITICAL)**: Figure 5A mismatch:
- **V5 text**: Figure 5A = Module 4 NMF program discovery within cell types, Moran's I validation, fibroblast coherence (I=0.28)
- **V4 legend**: Figure 5A = Cross-sample integration schematic (PCA + Harmony + hierarchical clustering)

The v5 text puts Module 4 content (program discovery) in 5A and Module 5 content (integration) starting at 5B. The v4 legends put Module 5 content starting at 5A. The v5 text effectively adds a new panel 5A for Module 4 results and shifts the old 5A-D to 5B-E, but only uses 5A-D.

**Issue #5 (WARNING)**: Figure 5B partial mismatch. Both text and legend agree on UMAP of 590 programs / 73 aligned. However:
- V5 text describes 5B as Harmony integration output with hierarchical clustering
- V4 legend describes 5B as UMAP colored by cell type with the "5 appearing in >50% of samples" detail
- The integration schematic content (PCA, Harmony, patient breakdown) from v4 legend 5A is what v5 references in 5B. So the content overlaps but the schematic vs. UMAP split may be off.

---

### Supplementary Figure References

| # | Manuscript Line | Text Context | What Text Says Panel Contains | Legend Match? | Status |
|---|---|---|---|---|---|
| 19 | L101 | "...T47D...did not recapitulate the MCF7 results **(Supplemental Figure S3)**" | T47D ELISA and IF results | V4 legends S3: "Parameter sensitivity analysis (lambda regularization)" | **CRITICAL** |
| 20 | L254 | "Supplementary Figure S3: T47D ELISA and immunofluorescence results." | T47D wet lab validation | V4 legends S3: Parameter sensitivity analysis | **CRITICAL** |

**Issue #6 (CRITICAL)**: Supplementary Figure S3 mismatch:
- **V5 text + v5 supplementary list**: S3 = T47D ELISA and immunofluorescence results
- **V4 figure legends**: S3 = Parameter sensitivity analysis (lambda regularization)

The v4 legends have a completely different supplementary figure numbering:
| V4 Legend | Content |
|-----------|---------|
| S1 | Extended Module 1-2 results for all Xenium regions |
| S2 | Per-cell-type benchmarking breakdown |
| S3 | Parameter sensitivity analysis (lambda regularization) |
| S4 | Per-sample Module 4 programs |
| S5 | Extended differential expression and pathway analysis |

The v5 manuscript only references S3 (T47D), and lists it in supplementary information. The v5 text does not reference S1, S2, S4, or S5 from the v4 legends. The v5 text also mentions old figures S1 (ddPCR) and S2 (runtime comparison) from the original manuscript, which are completely absent from v4 legends.

---

## B. Figure to Text Coverage (Every Panel in V4 Legends)

### Figure 1 Panels

| Panel | Legend Description | Referenced in V5 Text? | Status |
|-------|-------------------|----------------------|--------|
| 1A | Five-module pipeline schematic | Yes (L63) | OK |
| 1B | Spatial-native operations at each stage | Referenced as part of 1B-C (L67) but text describes interoperability, not spatial operations | **WARNING** |
| 1C | Resolution flexibility (Visium vs Xenium) | Referenced as part of 1B-C (L67) but text describes interoperability, not resolution | **WARNING** |

### Figure 2 Panels

| Panel | Legend Description | Referenced in V5 Text? | Status |
|-------|-------------------|----------------------|--------|
| 2A | Module 1 statistical gates | Yes (L71) | OK |
| 2B | Module 2 colocalization network | Yes (L71) | OK |
| 2C | Xenium single-cell demonstration | Yes (L77) | OK |
| 2D | Known marker clustering validation | Yes (L77) | OK |

### Figure 3 Panels

| Panel | Legend Description | Referenced in V5 Text? | Status |
|-------|-------------------|----------------------|--------|
| 3A | Two-pass deconvolution schematic | Yes (L81) | OK |
| 3B | Xenium benchmarking (Pearson, JSD, RMSE) | Yes (L83-85), but text includes simulated data not in legend | **WARNING** |
| 3C | GEX deconvolution benchmarking bar chart | Not explicitly referenced by panel letter | **WARNING** |
| 3D | Predicted vs ground truth scatter (r=0.39) | Not referenced in v5 text at all | **INFO** - Orphaned panel |

**Issue #7 (WARNING)**: Figure 3C is described in the legend as "Gene expression deconvolution (Pass 2) benchmarking...Pearson correlation between predicted and ground truth cell type-specific expression. Coverage percentage...CITEgeist achieves 100% spatial coverage across all cell types." The v5 text discusses GEX benchmarking (NRMSE 0.04 vs 0.16, scResolve comparison) but never cites a specific panel letter for this content.

**Issue #8 (INFO)**: Figure 3D (predicted vs ground truth scatter, r=0.39) is described in the v4 legends but never referenced in the v5 manuscript text. This panel is orphaned.

### Figure 4 Panels (V4 = Module 4 Program Discovery)

| Panel | Legend Description | Referenced in V5 Text? | Status |
|-------|-------------------|----------------------|--------|
| 4A | NMF program discovery schematic | **No** - v5 moved Module 4 content to Figure 5A | **CRITICAL** - Entire figure redefined |
| 4B | Program examples with Moran's I | **No** - v5 Figure 4 is now the midkine case study | **CRITICAL** |
| 4C | Moran's I validation box plot | **No** | **CRITICAL** |
| 4D | Summary statistics, bivariate Moran's I | **No** | **CRITICAL** |

**Issue #9 (CRITICAL)**: The entire v4 Figure 4 (Module 4 Program Discovery) content has been reassigned. In v5, Module 4 content moved to Figure 5A. The v4 Figure 4 legends (4A-D) describe content that no longer exists as "Figure 4" in v5. The v4 legends must be completely rewritten for Figure 4 to match the v5 midkine case study (panels A-H).

### Figure 5 Panels (V4 = Module 5 Integration)

| Panel | Legend Description | Referenced in V5 Text? | Status |
|-------|-------------------|----------------------|--------|
| 5A | Integration schematic (PCA, Harmony) | Not as 5A; this content is described in v5 text before the 5B reference | **CRITICAL** - Panel shifted |
| 5B | UMAP of 590 programs | Yes (L111) | OK (content matches) |
| 5C | Response-associated programs | Yes (L113) | OK |
| 5D | PyDESeq2 volcano plot | Yes (L113) | OK |

### Figure 6 Panels (V4 = Interoperability)

| Panel | Legend Description | Referenced in V5 Text? | Status |
|-------|-------------------|----------------------|--------|
| 6A | Workflow diagram: CITEgeist -> PyDESeq2/GSEApy/COMMOT | **No** - Figure 6 retired in v5 | **INFO** - Orphaned/Retired |
| 6B | Scanpy UMAP of deconvolved layers | **No** - Figure 6 retired in v5 | **INFO** - Orphaned/Retired |
| 6C | PyDESeq2 volcano from CITEgeist output | **No** - Figure 6 retired in v5 (DE content moved to 5D) | **INFO** - Orphaned/Retired |
| 6D | GSEApy pathway enrichment bars | **No** - Figure 6 retired in v5 | **INFO** - Orphaned/Retired |
| 6E | MDK discovery-to-validation workflow table | **No** - Figure 6 retired in v5 (MDK content in Fig 4) | **INFO** - Orphaned/Retired |

**Issue #10 (INFO)**: Figure 6 is entirely retired in v5. All 5 panels (6A-6E) are orphaned in the v4 legends. Some content was absorbed into other figures:
- 6A interoperability content -> partially in v5 text at L117 but not assigned to a panel
- 6C PyDESeq2 volcano -> became Figure 5D
- 6E MDK workflow -> absorbed into Figure 4 narrative
- 6B and 6D appear to be dropped entirely

### Supplementary Figure Panels

| Panel | V4 Legend Description | Referenced in V5 Text? | Status |
|-------|----------------------|----------------------|--------|
| S1 | Extended Module 1-2 for all Xenium regions | Not referenced | **INFO** - Orphaned or TBD |
| S2 | Per-cell-type benchmarking breakdown | Not referenced | **INFO** - Orphaned or TBD |
| S3 | Parameter sensitivity analysis | Not referenced; v5 S3 = T47D ELISA/IF | **CRITICAL** - Reassigned content |
| S4 | Per-sample Module 4 programs | Not referenced | **INFO** - Orphaned or TBD |
| S5 | Extended DE and pathway analysis | Not referenced | **INFO** - Orphaned or TBD |

---

## C. Number/Statistic Consistency

### Pearson Correlation Values

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| CITEgeist Xenium Pearson r | 0.60 +/- 0.05 (L85) | r = 0.60 (Fig 3B legend) | MATCH |
| Cell2Location Pearson r | 0.61 +/- 0.04 (L85) | r = 0.61 (Fig 3B legend) | MATCH |
| RCTD Pearson r | 0.62 +/- 0.03 (L85) | r = 0.62 (Fig 3B legend) | MATCH |
| Tangram Pearson r | 0.14 +/- 0.08 (L85) | Not in legend | INFO - Text only |
| Seurat Pearson r | 0.17 +/- 0.07 (L85) | Not in legend | INFO - Text only |
| GEX CITEgeist Pearson r | 0.44 across 7 cell types (L87) | Not in legend | INFO - Text only |
| scResolve Pearson r | 0.43 (L87) | Not in legend | INFO - Text only |
| Aggregate scatter r | Not mentioned in text | r = 0.39 (Fig 3D legend) | **WARNING** - Legend only |

**Issue #11 (WARNING)**: The v4 legend for Figure 3D reports an aggregate correlation of r = 0.39. This value never appears in the v5 text.

### Spot Counts

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| Total spots | "7,054 spots" (L85) | "n = 7,054 spots" (Fig 3B legend) | MATCH |
| Tissue regions | "5 tissue regions" (L85) | "5 tissue regions" (Fig 3B legend) | MATCH |

### Program Counts

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| Individual programs | 590 (L111) | 590 (Fig 5B legend) | MATCH |
| Aligned programs | 73 (L111) | 73 (Fig 5B legend) | MATCH |
| Conserved (>50% samples) | 5 (L111) | 5 (Fig 5B legend) | MATCH |
| Specimen-specific | 36 (L111) | Not in legend | INFO |

### Differential Expression Numbers

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| Genes tested | 13,371 (L113) | 13,371 (Fig 5D legend) | MATCH |
| Significant DE genes | 127 (L113) | 127 (Fig 5D legend) | MATCH |
| Up in progressors | 122 (L113) | 122 (Fig 5D legend) | MATCH |
| Up in responders | 5 (L113) | 5 (Fig 5D legend) | MATCH |

### Sample Counts

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| Total specimens | 14 (L107) | 14 (Fig 5A legend) | MATCH |
| Total patients | 6 (L107) | 6 (Fig 5A legend) | MATCH |
| Responder specimens | 5 (L107) | 5 (Fig 5A, 5D legends) | MATCH |
| Progressor specimens | 9 (L107) | 9 (Fig 5A legend) | MATCH |
| Responder patients | P1, P5 (implied L107) | P1, P5 (Fig 5A legend) | MATCH |
| Progressor patients | P2-P4, P6 (implied) | P2-P4, P6 (Fig 5A legend) | MATCH |

### Gene Names Referenced

| Gene Set | V5 Text | V4 Legend | Match? |
|----------|---------|-----------|--------|
| Progressor-upregulated MMP genes | MMP13, MMP3, ADAMTS4, ADAMTS15 (L113) | MMP13, MMP3, ADAMTS4, ADAMTS15 (5D legend) | MATCH |
| Other progressor genes | MLKL, ALOX5AP, CLEC5A (L113) | MLKL, ALOX5AP, CLEC5A (5D legend) | MATCH |
| Responder-upregulated genes | TMEM38B, TRIM72 (L113) | TMEM38B, TRIM72 (5D legend) | MATCH |
| Responder macrophage program | FABP4, HBA2, TNXB (L113) | FABP4, HBA2, TNXB (5C legend) | MATCH |
| Responder cancer luminal | MGP, MT-CO3, FOS (L113) | MGP, MT-CO3, FOS (5C legend) | MATCH |
| Progressor DC program | FTL, FN1, TIMP1 (L113) | FTL, FN1, TIMP1 (5C legend) | MATCH |
| MDK/PTN/MIF signaling | MDK, PTN, MIF (L95) | Not in current Fig 4 legend | **CRITICAL** - Legend mismatch |

**Issue #12 (CRITICAL)**: All gene names for the midkine case study (MDK, PTN, MIF, etc.) appear in the v5 text under Figure 4 panels but have no corresponding content in the v4 figure legends (since v4 Figure 4 is about Module 4 program discovery, not the midkine case study).

### Moran's I Thresholds and Percentages

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| Module 1 Moran's I threshold | I > 0.1, p < 0.05 (L73) | I > 0.1, p < 0.05 (Fig 2A legend) | MATCH |
| Module 4 coherence threshold (text) | I > 0.2, p < 0.01 (L109) | I > 0.15, p < 0.01 (Fig 4C legend) | **WARNING** |
| % programs exceeding threshold | 68% (L109) | 57% (100/175) (Fig 4C legend) | **WARNING** |
| Fibroblast mean I | 0.28 (L109) | 0.28 (Fig 4B legend) | MATCH |

**Issue #13 (WARNING)**: Moran's I threshold discrepancy:
- V5 text (L109): "programs with significant positive spatial autocorrelation (I > 0.2, p < 0.01)"
- V4 legend (Fig 4C): "Threshold line at I = 0.15; programs above are considered spatially coherent"
- The threshold is 0.20 in text vs 0.15 in legend.

**Issue #14 (WARNING)**: Percentage of spatially coherent programs:
- V5 text: "68% of discovered programs exceeded the threshold for moderate spatial coherence"
- V4 legend: "57% of discovered programs (100/175) exceed this threshold"
- 68% vs 57%. This could be explained by different datasets (v5 describes the full 14-specimen cohort; v4 describes 5 Xenium regions with 175 programs), but the discrepancy should be clarified since they reference different underlying analyses.

### Conserved Relationships

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| Total conserved relationships | 191 (L115) | 191 (Fig 5 note) | MATCH |
| Co-localized | 26 / 13.6% (L115) | 26 / 13.6% (Fig 5 note) | MATCH |
| Mutually exclusive | 6 / 3.1% (L115) | 6 / 3.1% (Fig 5 note) | MATCH |

### Benchmarking Statistics (RMSE/MAE/JSD)

| Statistic | V5 Text | V4 Legend | Match? |
|-----------|---------|-----------|--------|
| CITEgeist RMSE | 0.167 +/- 0.006 (L85) | Not in legend | INFO |
| CITEgeist MAE | 0.115 +/- 0.005 (L85) | Not in legend | INFO |
| Cell2Location RMSE | 0.179 +/- 0.017 (L85) | Not in legend | INFO |
| RCTD RMSE | 0.177 +/- 0.004 (L85) | Not in legend | INFO |
| Simulated CITEgeist r | 0.95 (L83) | Not in legend | INFO |
| Simulated JSD | 0.16 (L83) | Not in legend | INFO |
| Wilcoxon p vs C2L | 0.31 (L85) | Not in legend | INFO |
| Wilcoxon p vs RCTD | 0.19 (L85) | Not in legend | INFO |

---

## D. V5 vs Old Unversioned Figures Comparison

### Content Migration Map

| Old Figure | Old Content | V5 Figure | Status |
|-----------|-------------|-----------|--------|
| **Fig 1** | Graphical abstract (trial, CITEgeist approach, simulation, wet lab validation) | **Fig 1** | REORGANIZED - Now pipeline overview schematic only (no trial/simulation/wetlab overview) |
| **Fig 2** | Simulated benchmarking (proportions) | **Fig 3B** (partial) | MERGED - Simulation results folded into benchmarking section; text mentions simulated results but no dedicated figure panel in v4 legends |
| **Fig 3** | GEX deconvolution metrics | **Fig 3C** | MERGED into unified benchmarking figure |
| **Fig 4** | Real cancer data validity (spatial plots, CellMarker, COMMOT) | **Fig 4A-F** (partial) | REORGANIZED - Spatial validation + COMMOT now part of midkine case study |
| **Fig 5** | Wet lab validation (ESR1, mutation sig, pathways, COMMOT, ELISA, IF) | **Fig 4A-H** | MERGED with old Fig 4 into unified midkine case study figure |
| **Fig S1** | ddPCR | Not referenced | **MISSING** - ddPCR supplementary figure dropped from v5 |
| **Fig S2** | Runtime comparison | Not referenced | **MISSING** - Runtime comparison dropped from v5 |
| **Fig S3** | T47D ELISA/IF | **Supp Fig S3** | RETAINED |

### New Content in V5 Not Present in Old Figures

| V5 Content | Description | Status |
|------------|-------------|--------|
| **Fig 2 (all)** | Module 1-2 automated profile discovery pipeline | NEW - Not in old manuscript |
| **Fig 5A** | Module 4 NMF program discovery with Moran's I validation | NEW - Modules 4+5 are entirely new |
| **Fig 5B** | Cross-sample UMAP integration (590 programs) | NEW |
| **Fig 5C** | Response-associated programs (responder vs progressor) | NEW |
| **Fig 5D** | PyDESeq2 differential expression volcano | NEW |
| V4 **Fig 4** (legends) | Module 4 spatial program discovery (NMF schematic, Moran's I box plot) | NEW - Did not exist in old manuscript |
| V4 **Fig 6** (legends) | Interoperability demonstration (scanpy, PyDESeq2, GSEApy, MDK workflow) | NEW but RETIRED in v5 |

### Content from Old Figures Missing in V5

| Old Content | Description | Assessment |
|-------------|-------------|------------|
| **Fig 1 graphical abstract** | Trial design overview, four-quadrant graphical abstract | DROPPED - Replaced by pipeline schematic. May be intentional for CRM format. |
| **Fig S1 ddPCR** | Droplet digital PCR validation | DROPPED - ddPCR mentioned in text (L91) as confirmation method but no figure. Consider whether reviewer may request this. |
| **Fig S2 runtime** | Computational runtime comparison across methods | DROPPED - No runtime benchmarking in v5. This was a competitive advantage; consider restoring as supplementary. |
| **CellMarker validation** | CellMarker database validation from old Fig 4 | DROPPED - Replaced by HPA secretion correlation (Spearman rho = 0.46) mentioned in text. |

### Supplementary Figure Comparison

| V5 Text Supp | Old Manuscript Supp | V4 Legends Supp | Status |
|--------------|--------------------|--------------------|--------|
| S3: T47D ELISA/IF | S3: T47D ELISA/IF | S3: Lambda sensitivity | THREE-WAY CONFLICT |
| *(none)* | S1: ddPCR | S1: Extended Mod 1-2 | CONFLICT |
| *(none)* | S2: Runtime | S2: Per-cell benchmarking | CONFLICT |
| *(none)* | *(none)* | S4: Per-sample Mod 4 | V4 only |
| *(none)* | *(none)* | S5: Extended DE/pathway | V4 only |

**Issue #15 (CRITICAL)**: All three versions (old manuscript, v5 text, v4 legends) have different supplementary figure numbering and content. The supplementary figures need to be reconciled across all three.

---

## E. Panel Letter Continuity

### Figure 1: A, B, C
- V4 legends: A, B, C -- **Continuous, no gaps**
- V5 text references: 1A, 1B-C -- **OK**

### Figure 2: A, B, C, D
- V4 legends: A, B, C, D -- **Continuous, no gaps**
- V5 text references: 2A-B, 2C, 2D -- **OK**

### Figure 3: A, B, C, D
- V4 legends: A, B, C, D -- **Continuous, no gaps**
- V5 text references: 3A, 3B only -- **WARNING: 3C and 3D not cited in text**

### Figure 4 (V4 legends): A, B, C, D
- V4 legends: A, B, C, D -- **Continuous, no gaps** (but wrong content for v5)
- V5 text references: 4A, 4B, 4C, 4D-E, 4F, 4G, 4H -- **Continuous A-H, no gaps**
- **CRITICAL**: V4 legends only have 4A-D; v5 text uses 4A-H. Legends need 4 new panels (E-H).

### Figure 5 (V4 legends): A, B, C, D
- V4 legends: A, B, C, D -- **Continuous, no gaps**
- V5 text references: 5A, 5B, 5C, 5D -- **Continuous, no gaps**
- **But**: V5 5A (Module 4 programs) does not correspond to v4 legend 5A (Integration schematic). Content shifted.

### Figure 6 (V4 legends): A, B, C, D, E
- V4 legends: A, B, C, D, E -- **Continuous, no gaps**
- V5 text: No references to Figure 6 -- **Entire figure retired**

---

## Summary of All Issues

### CRITICAL Issues (Must Fix)

| # | Issue | Location | Description |
|---|-------|----------|-------------|
| 3 | Figure 4 complete mismatch | Fig 4 | V5 text: midkine case study (4A-H); V4 legends: Module 4 program discovery (4A-D). Legends must be completely rewritten. |
| 4 | Figure 5A mismatch | Fig 5A | V5 text: Module 4 programs; V4 legends: Integration schematic. Content shifted by one panel. |
| 6 | Supplementary S3 conflict | Supp S3 | V5: T47D ELISA/IF; V4 legends: Lambda sensitivity analysis. Different content. |
| 9 | V4 Figure 4 orphaned | Fig 4 | Entire v4 Figure 4 (Module 4) content no longer matches any v5 figure. Needs reassignment. |
| 12 | Midkine genes missing from legends | Fig 4 | MDK, PTN, MIF and all case study content absent from v4 legends. |
| 15 | Supplementary figure numbering | All Supp | Three-way conflict between old manuscript, v5 text, and v4 legends for all supplementary figures. |

### WARNING Issues (Should Fix)

| # | Issue | Location | Description |
|---|-------|----------|-------------|
| 1 | Figure 1B-C content mismatch | Fig 1B-C | Text says interoperability; legend says spatial operations + resolution flexibility. |
| 2 | Simulated benchmarking missing panel | Fig 3B | Text describes simulation results under 3B but legend only covers Xenium. No panel for simulation benchmarks. |
| 5 | Figure 5B partial mismatch | Fig 5B | Schematic vs UMAP content split differs between text and legend. |
| 7 | Figure 3C unreferenced | Fig 3C | GEX benchmarking panel exists in legend but not cited by panel letter in text. |
| 11 | r = 0.39 aggregate scatter | Fig 3D | Value in legend only, never appears in text. |
| 13 | Moran's I threshold discrepancy | Fig 4C/5A | Text: I > 0.2; Legend: I > 0.15. Different thresholds. |
| 14 | % coherent programs discrepancy | Fig 4C/5A | Text: 68%; Legend: 57% (100/175). Different datasets but should clarify. |

### INFO Issues (Minor / Consider)

| # | Issue | Location | Description |
|---|-------|----------|-------------|
| 8 | Figure 3D orphaned | Fig 3D | Aggregate scatter not referenced in text. |
| 10 | Figure 6 entirely retired | Fig 6 | All 5 panels orphaned. Some content absorbed elsewhere. |
| -- | ddPCR figure dropped | Supp S1 (old) | ddPCR validation figure from old manuscript not in v5. |
| -- | Runtime figure dropped | Supp S2 (old) | Runtime comparison from old manuscript not in v5. |

---

## Recommended Actions

### Immediate (Pre-Submission)

1. **Rewrite v4 Figure 4 legends entirely** to match v5 midkine case study (panels A-H covering basal signatures, spatial distribution, EstroGene, pathways, COMMOT, ELISA, IF).

2. **Rewrite v4 Figure 5A legend** to cover Module 4 NMF program discovery content, and shift integration schematic to a combined 5A/5B or separate the Module 4 content into the figure design.

3. **Fix Figure 1B-C legend or text**: Either update the text to describe what 1B-C actually shows (spatial operations + resolution flexibility) or redesign the panels to include interoperability content.

4. **Reconcile supplementary figures**: Decide on final numbering. V5 text only references S3 (T47D). Determine which of the v4 supplementary panels (S1-S5) should be kept, renumbered, or dropped.

5. **Add panel references for 3C and 3D** in text, or remove panels from figure if they are no longer needed.

6. **Reconcile Moran's I thresholds**: Decide whether the threshold is 0.15 or 0.20 and make consistent across text and legends.

### Consider for Revision

7. **Add simulated benchmarking panel** (or relabel existing panels) so that scCube simulation results have a clear figure reference.

8. **Consider restoring runtime comparison** as supplementary material (competitive advantage over reference-based methods).

9. **Verify whether r = 0.39 aggregate scatter** (Fig 3D) is still desired or should be removed.

10. **Add explicit "(Figure 3C)" citation** in the gene expression deconvolution benchmarking paragraph.

---

*Generated: 2026-02-08*
