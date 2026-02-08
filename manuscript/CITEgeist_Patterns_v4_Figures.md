# CITEgeist Patterns Manuscript - Figures and Legends

**Manuscript**: CITEgeist: A Spatially-Native Framework for Multi-Modal Integration and Program Discovery in Spatial Transcriptomics

**Figure files**: `manuscript/figures/output/`

---

## Figure 1: CITEgeist Framework Overview

**File**: `figure1_pipeline_overview.pdf`

![Figure 1](figures/output/figure1_pipeline_overview.png)

**Legend:**
(A) Schematic of the five-module pipeline: Module 1 (Marker Interest Detection), Module 2 (Profile Assembly), Module 3 (Spatial Deconvolution), Module 4 (Program Discovery), Module 5 (Cross-Sample Integration). Arrows indicate data flow between modules. Input data consists of spatial transcriptomics with matched protein measurements (CITE-seq antibody capture). Each module produces interpretable outputs stored in standard AnnData format.
(B) Spatial-native operations highlighted at each stage: Moran's I spatial autocorrelation for marker detection (Module 1), colocalization networks for profile discovery (Module 2), Laplacian regularization for proportion estimation (Module 3), spatial coherence validation for programs (Module 4), and conserved relationship analysis across samples (Module 5).
(C) Resolution flexibility: the same framework operates on spot-level data (Visium, ~55 μm diameter, ~10-30 cells per spot) and single-cell resolution data (Xenium, CosMx). Module implementations automatically adapt to data resolution through spatial neighbor definitions that scale with technology.

---

## Figure 2: Modules 1-2 - Automated Profile Discovery

**File**: `figure2_profile_discovery.pdf`

![Figure 2](figures/output/figure2_profile_discovery.png)

**Legend:**
(A) Module 1: Marker interest detection using three statistical gates. Moran's I spatial autocorrelation identifies markers with organized spatial patterns (threshold I > 0.1, p < 0.05). Gaussian mixture modeling separates signal from background, computing signal-to-noise ratios (threshold SNR > 1.5). Kurtosis analysis detects peaked distributions indicative of localized expression (threshold κ > 2.0). Markers passing either the Moran's I gate or the kurtosis gate (combined with the GMM filter) are forwarded to Module 2.
(B) Module 2: Profile assembly from spatial colocalization. For each marker pair, we compute same-spot co-occurrence, expression correlation, adjacent-spot enrichment, and bivariate Moran's I. Significance-filtered relationships (permutation p < 0.05) form a colocalization network where connected components represent distinct cellular lineages. Hierarchical clustering with dynamic tree cutting extracts marker profiles, with optimal count determined by reconstruction accuracy.
(C) Xenium single-cell demonstration: Modules 1-2 applied to single-cell resolution data, showing that spatial colocalization analysis discovers biologically coherent profiles without reference data.
(D) Validation: discovered profiles correctly recover known cell type markers. CD3E and CD8A cluster together (T cells), CD68 and CD163 cluster together (macrophages), EPCAM and KRT markers cluster together (epithelial cells). This demonstrates that spatial colocalization provides sufficient signal for cell type definition.

---

## Figure 3: Module 3 - Spatial Deconvolution Benchmarking

**File**: `figure3_benchmarking.pdf`

![Figure 3](figures/output/figure3_benchmarking.png)

**Legend:**
(A) Two-pass deconvolution schematic. Pass 1 estimates cell type proportions from protein markers using quadratic programming with Laplacian regularization to enforce spatial coherence across neighborhoods. Neighborhood-aware finetuning adjusts estimates based on local context. Pass 2 deconvolves bulk gene expression into cell type-specific layers using the estimated proportions as constraints, producing per-spot, per-cell-type expression matrices.
(B) Benchmarking on Xenium pseudo-Visium spots. Single-cell Xenium data was aggregated into pseudo-Visium spots (55 μm diameter) to create ground truth cell type proportions. Methods compared: CITEgeist (same-slide proteomics, no reference), Cell2Location (reference-based), RCTD (reference-based), Tangram (label transfer), Seurat (label transfer). Metrics: Pearson correlation with ground truth (higher is better), Jensen-Shannon divergence (lower is better), RMSE (lower is better). CITEgeist achieves r = 0.60, comparable to reference-based methods (Cell2Location r = 0.61, RCTD r = 0.62), without requiring external reference data. n = 7,054 spots across 5 tissue regions.
(C) Gene expression deconvolution (Pass 2) benchmarking. Bar chart comparing CITEgeist performance across cell types, showing Pearson correlation between predicted and ground truth cell type-specific expression. Coverage percentage indicates proportion of spots with valid predictions for each cell type. CITEgeist achieves 100% spatial coverage across all cell types.
(D) Predicted versus ground truth scatter plot. Aggregated cell type proportions across all spots and cell types, demonstrating overall concordance (r = 0.39). Points colored by cell type.

---

## Figure 4: Module 4 - Spatial Program Discovery

**File**: `figure4_module4_programs.pdf`

![Figure 4](figures/output/figure4_module4_programs.png)

**Legend:**
(A) NMF-based program discovery schematic. For each cell type, the deconvolved expression matrix (spots × genes) weighted by cell type proportion is factorized using non-negative matrix factorization into program loadings (W matrix, genes × programs) and program activities (H matrix, spots × programs). Moran's I spatial autocorrelation on program activities validates that discovered programs exhibit spatial coherence.
(B) Program examples across cell types. Top genes for selected high-coherence programs with associated Moran's I values. Fibroblast programs show highest spatial coherence (mean I = 0.28), consistent with organized stromal architecture. Example: fibroblast program with HLA-DRA, FCGR3A, HLA-DRB1 (I = 0.37) suggests antigen-presenting fibroblast subpopulations.
(C) Moran's I validation. Box plot showing spatial coherence (Moran's I values) by cell type. Threshold line at I = 0.15; programs above are considered spatially coherent. 57% of discovered programs (100/175) exceed this threshold. Fibroblast programs show highest spatial coherence, consistent with organized stromal architecture. Xenium single-cell data, 5 regions, 7 cell types, 175 total programs.
(D) Summary statistics table. Per-cell-type breakdown showing number of programs (N), mean Moran's I, and maximum Moran's I for each cell type. Bivariate Moran's I is computed between all program pairs to identify co-localized programs (positive I) and mutually exclusive programs (negative I), providing a spatial interaction map complementing pathway analysis.

---

## Figure 5: Module 5 - Cross-Sample Integration

**File**: `figure5_integration.pdf`

![Figure 5](figures/output/figure5_integration.png)

**Legend:**
(A) Integration schematic. Programs from all samples (n = 14 breast cancer samples from 6 patients) are embedded in a shared latent space using PCA. Harmony batch correction aligns programs across samples, removing patient-specific technical variation. Hierarchical clustering identifies aligned programs—conserved transcriptional states appearing across patients. Response analysis compares program enrichment between responders (5 samples from patients P1, P5) and progressors (9 samples from patients P2-P4, P6).
(B) UMAP visualization of 590 individual programs from 14 samples, colored by cell type. Successful Harmony integration groups programs by biological identity rather than sample of origin. 73 aligned programs identified, with 5 appearing in >50% of samples.
(C) Response-associated programs. Bar plot showing 3 responder-enriched and 4 progressor-enriched aligned programs. Responder-enriched: macrophage programs with FABP4/HBA2/TNXB and VIM/TMSB4X/S100A6, cancer luminal program with MGP/MT-CO3/FOS. Progressor-enriched: dendritic cell program with FTL/FN1/TIMP1, cancer basal mitochondrial program, cancer luminal program with MGP/HBA2/LTF, endothelial program with KLHL17/MT-ND1/GPR101.
(D) PyDESeq2 differential expression volcano plot. Pseudo-bulk aggregation of deconvolved gene expression across 14 samples (5 responder, 9 progressor). 13,371 genes tested, 127 significant (adjusted p < 0.05). 122 genes upregulated in progressors (including MMP13, MMP3, ADAMTS4, ADAMTS15, MLKL, ALOX5AP, CLEC5A), 5 genes upregulated in responders (including TMEM38B, TRIM72). Top genes labeled.

*Module 5 also identifies conserved spatial relationships across samples: of 191 conserved relationships detected, 26 (13.6%) show consistent co-localization, 6 (3.1%) show consistent mutual exclusion, highlighting multicellular organization patterns persisting across patients.*

---

## Figure 6: Interpretable Outputs and Interoperability

**File**: `figure6_interoperability.pdf`

![Figure 6](figures/output/figure6_interoperability.png)

**Legend:**
(A) Workflow diagram showing CITEgeist outputs feeding into standard bioinformatics tools. Cell type proportions and deconvolved GEX layers (Module 3) enable PyDESeq2 differential expression within cell populations. Spatial programs (Module 4) enable GSEApy/Enrichr pathway enrichment. Program relationships (Module 5) enable network analysis with COMMOT and CellPhoneDB. All outputs in standard AnnData/CSV/Parquet formats compatible with scanpy workflows.
(B) Scanpy cluster visualization. UMAP embedding of CITEgeist deconvolved expression layers, colored by inferred cell type (Epithelial, Fibroblasts, Macrophages, Myofibroblast, Adipocyte, Endothelial). Demonstrates seamless integration with standard scanpy workflows for visualization and clustering.
(C) PyDESeq2 differential expression. Volcano plot showing differentially expressed genes from CITEgeist-derived cell type-stratified expression. CITEgeist outputs feed directly into standard DE frameworks, enabling discovery of condition-specific signatures within cell populations.
(D) GSEApy pathway enrichment. Horizontal bar chart showing enriched pathways from CITEgeist program gene loadings. Top enriched pathways include extracellular matrix organization, collagen metabolism, and cell-matrix adhesion, consistent with stromal remodeling in breast cancer.
(E) Discovery-to-validation workflow example: Midkine (MDK). Table showing how CITEgeist spatial transcriptomic discoveries can identify clinically actionable targets. Midkine was identified through spatial program analysis, validated via ChIP-seq and ELISA, and linked to age-related mammary tumorigenesis mechanisms, demonstrating the translational potential of CITEgeist outputs.

---

## Supplementary Figures

*(To be generated)*

- **Figure S1**: Extended Module 1-2 results for all Xenium regions
- **Figure S2**: Per-cell-type benchmarking breakdown
- **Figure S3**: Parameter sensitivity analysis (lambda regularization)
- **Figure S4**: Per-sample Module 4 programs
- **Figure S5**: Extended differential expression and pathway analysis

