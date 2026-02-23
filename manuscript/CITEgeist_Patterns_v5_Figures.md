# CITEgeist Patterns Manuscript - Figures and Legends (v5)

**Manuscript**: CITEgeist: Accurate deconvolution of spatial transcriptomics with same-slide proteomics reveals midkine as a secreted microenvironment modulator in ESR1 mutant breast cancer

**Figure files**: `manuscript/figures/output/`

---

## Figure 1: CITEgeist Framework Overview

**File**: `figure1_pipeline_overview.pdf`

![Figure 1](figures/output/figure1_pipeline_overview.png)

**Legend:**
(A) Schematic of the five-module pipeline: Module 1 (Marker Interest Detection), Module 2 (Profile Assembly), Module 3 (Spatial Deconvolution), Module 4 (Program Discovery), Module 5 (Cross-Sample Integration). Arrows indicate data flow between modules. Input data consists of spatial transcriptomics with matched protein measurements (CITE-seq antibody capture). Each module produces interpretable outputs stored in standard AnnData format.
(B) Spatial-native operations highlighted at each stage: Moran's I spatial autocorrelation for marker detection (Module 1), colocalization networks for profile discovery (Module 2), Laplacian regularization for proportion estimation (Module 3), spatial coherence validation for programs (Module 4), and conserved relationship analysis across samples (Module 5). All outputs in standard AnnData/CSV/Parquet formats compatible with scanpy workflows and interoperable with downstream tools such as PyDESeq2, GSEApy, and COMMOT.
(C) Resolution flexibility: the same framework operates on spot-level data (Visium, ~55 um diameter, ~10-30 cells per spot) and single-cell resolution data (Xenium, CosMx). Module implementations automatically adapt to data resolution through spatial neighbor definitions that scale with technology.

---

## Figure 2: Modules 1-2 - Automated Profile Discovery

**File**: `figure2_profile_discovery.pdf`

![Figure 2](figures/output/figure2_profile_discovery.png)

**Legend:**
(A) Module 1: Marker interest detection using three statistical gates. Moran's I spatial autocorrelation identifies markers with organized spatial patterns (threshold I > 0.1, p < 0.05). Gaussian mixture modeling separates signal from background, computing signal-to-noise ratios (threshold SNR > 1.5). Kurtosis analysis detects peaked distributions indicative of localized expression (threshold kappa > 2.0). Markers passing either the Moran's I gate or the kurtosis gate (combined with the GMM filter) are forwarded to Module 2.
(B) Module 2: Profile assembly from spatial colocalization. For each marker pair, we compute same-spot co-occurrence, expression correlation, adjacent-spot enrichment, and bivariate Moran's I. Significance-filtered relationships (permutation p < 0.05) form a colocalization network where connected components represent distinct cellular lineages. Hierarchical clustering with dynamic tree cutting extracts marker profiles, with optimal count determined by reconstruction accuracy.
(C) Xenium single-cell demonstration: Modules 1-2 applied to single-cell resolution data, showing that spatial colocalization analysis discovers biologically coherent profiles without reference data.
(D) Validation: discovered profiles correctly recover known cell type markers. CD3E and CD8A cluster together (T cells), CD68 and CD163 cluster together (macrophages), EPCAM and KRT markers cluster together (epithelial cells). This demonstrates that spatial colocalization provides sufficient signal for cell type definition.

---

## Figure 3: Benchmarking and Validation

**File**: `figure3_benchmarking.pdf`

![Figure 3](figures/output/figure3_benchmarking.png)

**Legend:**
(A) Profile discovery accuracy on Xenium RCC pseudo-Visium data. Modules 1-2 automatically discover cell type protein marker profiles without prior knowledge. Table comparing discovered profiles to known marker assignments for 7 cell types: B cells (CD20), CD4+ T cells (CD3E/CD4), CD8+ T cells (CD3E/CD8A), macrophages (CD68/CD163), endothelial cells (CD31), epithelial cells (PanCK), and fibroblasts (alphaSMA). 6 of 7 profiles exactly match known assignments; 1 is recovered as a superset.
(B) Xenium proportion benchmark. Cell type proportion accuracy on Xenium breast cancer pseudo-Visium spots (55 um diameter) across 5 tissue regions (n = 7,054 spots). Five methods compared: CITEgeist (same-slide proteomics, no reference), Cell2Location (reference-based), RCTD (reference-based), Tangram (label transfer), Seurat (label transfer). Three metrics shown (higher is better): Pearson r, 1-JSD, 1-RMSE. CITEgeist achieves r = 0.600 +/- 0.061, comparable to Cell2Location (r = 0.612) and RCTD (r = 0.619), without requiring external reference data. Tangram (r = 0.143) and Seurat (r = 0.173) perform substantially worse.
(C) Simulated proportion benchmark. scCube-generated synthetic Visium-like data from an ER+ breast cancer scRNA-seq atlas (Wu et al., 2021; 9 cell types, 5 replicates). Two tissue conditions: "high_seg" (homogeneous composition) and "mixed" (heterogeneous, cancer-like). CITEgeist maintains consistent performance across both conditions, while reference-based methods show degraded accuracy on mixed tissue.
(D) Gene expression deconvolution benchmark. Bar chart comparing GEX RMSE for the three methods capable of producing per-cell-type expression estimates (CITEgeist, Cell2Location, Tangram) across both simulated conditions. Fold-difference annotations indicate CITEgeist's advantage: 1.7x over Cell2Location and 4.3x over Tangram in high_seg; 5.1x over Cell2Location and 5.6x over Tangram in mixed tissue.
(E) Predicted versus ground truth scatter plot. Per-spot predicted versus actual cell type proportions for one Xenium tissue region, with points colored by cell type, diagonal reference line (dashed gray), and linear regression line (solid black) shown.
(F) Spatial comparison maps. 2x3 grid showing ground truth (top row) versus CITEgeist predictions (bottom row) for epithelial, macrophage, and T cell populations, demonstrating faithful recapitulation of spatial distributions.

---

## Figure 4: CITEgeist Identifies Upregulated Midkine Signaling in ESR1-Mutant Breast Cancer

**File**: `figure4_midkine_esr1.pdf`

![Figure 4](figures/output/figure4_midkine_esr1.png)

**Legend:**
(A) Basal cytokeratin spatial map on H&E tissue backdrop. Sample P4-S2_1i_rep showing cytokeratin expression across spots in the deconvolved cancer cell layer, overlaid on the histology image, consistent with prior work on ESR1-mutant breast cancers (Li et al., 2022).
(B) D538G mutation signal spatial map on H&E. The EstroGene signature is dispersed across the tissue section rather than clustered in a single location as might be expected from a single subclone.
(C) EstroGene validation. Left: split violin plot comparing EstroGene scores between D538G-high and wild-type spots (Mann-Whitney p < 0.0001). Right: 19-gene heatmap (15 upregulated, 4 downregulated EstroGene genes) with row-wise z-scores, confirming the ESR1 D538G molecular identity from the EstroGene database (Li et al., 2024).
(D) Combined pathway dot plot. Enrichr over-representation analysis (MSigDB_Hallmark_2020 + KEGG_2021_Human) of 165 significantly upregulated genes in D538G-positive cancer spots. Dot size represents gene overlap count, color represents -log10(p-value), and x-axis shows Enrichr combined score. Top pathways: Estrogen Response Early (538.37), Apoptosis (493.42), TNF-alpha Signaling via NF-kB (466.91), and Epithelial Mesenchymal Transition (401.97), consistent with known downstream effects of ESR1 D538G from prior work (Li et al., 2022).
(E) MDK-SDC4 ligand-receptor communication score on H&E spatial map. COMMOT-derived midkine signaling scores overlaid on tissue histology, showing that MDK-SDC4 communication spatially corresponds to D538G-high regions identified in panel B.
(F) COMMOT signaling analysis of D538G versus wild-type regions. Bar plot showing 5 significant ligand-receptor pairs with fold-change annotations: MDK-SDC4 (3.4x), MDK-NCL (2.7x), PTN-SDC4 (21.7x), PTN-NCL (20.3x), and MIF-CD74_CD44 (5.4x). All pairs FDR < 0.001 (Mann-Whitney U test with Benjamini-Hochberg correction). CITEgeist-derived secretion signals were validated against established cell type-specific secretion profiles from the Human Protein Atlas (Spearman rho = 0.46, p = 2.86e-37).
(G) ELISA validation of midkine secretion. Conditioned media from MCF7 wild-type versus D538G ESR1-mutant cell lines in estrogen-deprived and estradiol-supplemented (E2) conditions. The mutant cells secrete approximately double the midkine compared to wild-type in both conditions (p < 0.0001, unpaired t-test), consistent with constitutive receptor activity conferred by the D538G mutation.
(H) Immunofluorescence validation of midkine localization. MCF7 cells in estrogen-deprived environment mimicking a patient on endocrine therapy. D538G-mutant cell lines show approximately double the midkine at the cell membrane (p < 0.001, Mann-Whitney U test) and elevated intracellular midkine (p < 0.01), suggesting the mutation increases both midkine production and retention at the cell surface.

---

## Figure 5: Full Pipeline Application Reveals Endocrine Therapy Response Signatures

**File**: `figure5_full_pipeline.pdf`

![Figure 5](figures/output/figure5_full_pipeline.png)

**Legend:**
(A) Cell type proportion shift from biopsy to surgical specimens, grouped by response status. Grouped bar chart with blue bars for responders (n = 4 patients) and red bars for progressors (n = 2 patients), with individual patient values overlaid as dots. Progressors lose CD4 T cells (-30.7%) and gain cancer cells (+28.0% luminal in P1, +11.7% basal in P4), while responders gain fibroblasts (+7.8% average; P3: +31%) and lose cancer cells (-5.6%), reflecting divergent remodeling trajectories under endocrine therapy.
(B) Program conservation heatmap. Binary/intensity heatmap of 65 aligned programs (rows) across 12 samples (columns), with rows grouped by cell type (sidebar color) and columns annotated by patient, timepoint, and response status. Cancer_Basal programs are the most heterogeneous (13 aligned programs, mostly patient-specific), while aligned_032 (CD8 T cells) is perfectly conserved across all 12 samples. n = 490 individual programs aligned into 65 consensus programs via Harmony batch correction and hierarchical clustering.
(C) Response enrichment dot plot. Programs grouped by cell type on the y-axis, fraction of responder samples on the x-axis. Dot size indicates the number of samples containing each program; dot color indicates cell type. Four programs with significant differential enrichment are highlighted: one responder-enriched macrophage program (aligned_004; FABP4/HBA2/TNXB) and three progressor-enriched programs (aligned_008, CD4 T cells; aligned_002, Cancer Luminal; aligned_016, Monocytes). 8 responder specimens (4 patients) versus 4 progressor specimens (2 patients).
(D) Conserved relationship network. NetworkX spring layout of 65 program nodes colored by cell type, connected by 163 conserved spatial relationship edges. Green edges indicate co-localized programs (bivariate Moran's I > 0.15), red edges indicate mutually exclusive programs (I < -0.15), and gray edges indicate spatially independent programs. Key relationships labeled, including Fibroblast-CD4 T cell co-localization (I = 0.358) and Cancer-CD4 T cell exclusion (I = -0.194).
(E) Bivariate spatial co-localization in a responder (P3-S2). Two spatial maps showing fibroblast proportion and CD4 T cell proportion across tissue spots, with scatter inset showing positive spatial correlation. Bivariate Moran's I = 0.358, indicating that T cells and fibroblasts are spatially coordinated in responding tumors, suggesting that fibroblastic remodeling facilitates immune infiltration.
(F) Bivariate spatial exclusion in a progressor (P1-S2). Two spatial maps showing cancer luminal proportion and CD4 T cell proportion across tissue spots, with scatter inset showing negative spatial trend. Bivariate Moran's I = -0.194, indicating that T cells are excluded from tumor-dense regions in progressing tumors, consistent with immune evasion. This spatial organizational difference between responders and progressors cannot be detected by bulk or single-cell RNA-seq.
(G) Moran's I spatial coherence validation. Violin and strip plot showing Moran's I values by cell type across 490 NMF programs from 12 samples. Dashed line at threshold I = 0.15 (p < 0.01). Across the cohort, 47% of programs (231/490) exceed the spatial coherence threshold. CD4 T cell programs show the highest mean coherence (I = 0.33), followed by cancer and fibroblast programs (I = 0.25-0.26). B cell and endothelial programs show the lowest coherence (mean I = 0.08). This validates that NMF programs capture spatially organized biology rather than random noise.
(H) Spatial proportion comparison maps. 2x2 grid showing fibroblast and cancer luminal proportions for a representative responder (P3-S2) and progressor (P1-S2) on tissue spatial coordinates. Responder tissue shows fibrotic scarring with reduced cancer burden; progressor tissue shows tumor expansion with immune exclusion. Visual confirmation of the proportion shifts quantified in panel A and the spatial relationships demonstrated in panels E-F.

---

## Supplementary Figures

### Supplementary Figure S1: ddPCR Validation

*(ddPCR confirmation of ESR1 D538G mutation status)*

### Supplementary Figure S2: Differential Expression and Pathway Enrichment

**File**: `supp_figure2_de_pathway.pdf`

![Supplementary Figure S2](figures/output/supp_figure2_de_pathway.png)

**Legend:**
(A) PyDESeq2 differential expression volcano plot. Pseudo-bulk aggregation of deconvolved gene expression across 12 samples (8 responder from 4 patients, 4 progressor from 2 patients) with design formula ~condition + timepoint, where timepoint (biopsy vs surgical) accounts for systematic pre- vs post-treatment expression differences. 13,350 genes tested, 203 significant (adjusted p < 0.05). 120 genes upregulated in responders (including NEDD9, TMC5, GP2, ZNF655, DACH1, ACTG2), 83 genes upregulated in progressors (including VPREB3, COL4A4, CCNE1, FAM111B, PLK1). Top genes labeled by significance and fold change.
(B) Pathway enrichment dot plot. MSigDB Hallmark 2020 gene set over-representation analysis of differentially expressed genes, separated by direction. Responder-upregulated genes are enriched for Epithelial Mesenchymal Transition and TNF-alpha Signaling via NF-kB, suggesting stromal remodeling and immune activation in responding tumors. Progressor-upregulated genes are enriched for E2F Targets and G2-M Checkpoint, consistent with proliferative advantage in progressing tumors. Dot size represents gene overlap count; color represents -log10(adjusted p-value).

### Supplementary Figure S3: Runtime and Computational Performance

*(Runtime benchmarking across methods and dataset sizes)*

### Supplementary Figure S4: T47D Validation Results

**Legend:**
T47D ELISA and immunofluorescence results for midkine signaling in ESR1 D538G-mutant versus wild-type cells. T47D did not recapitulate the MCF7 results, consistent with previous work establishing that ESR1 mutation effects are context-specific (Arnesen et al., 2021; Li et al., 2022). MCF7 and T47D have different baseline ER expression levels and chromatin accessibility, leading to opposite downstream effects of the same mutation. These results demonstrate the cell line-specific nature of ESR1 mutation consequences and underscore the importance of validating computational findings across multiple experimental systems.
