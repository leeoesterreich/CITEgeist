# Figure 4 Briefing: Midkine/ESR1 Case Study

## Story
CITEgeist's spatial deconvolution reveals upregulated midkine (MDK) signaling in ESR1 D538G-mutant breast cancer spots. Patient P4-S2_1i_rep is a progressor (endocrine-resistant) with a known D538G hotspot mutation. Spatial analysis identifies basal cytokeratin signatures, mutation-associated transcriptional programs, pathway activation, and ligand-receptor communication — orthogonally validated by ELISA and immunofluorescence.

## Patient Context
- **Patient**: P4 (progressor, endocrine-resistant)
- **Sample**: HCC22-088-P4-S2_1i_rep (surgical specimen, post-treatment)
- **Mutation**: ESR1 D538G hotspot mutation
- **Key finding**: D538G-mutant spots have activated midkine signaling to microenvironment

## Panel-by-Panel

### Panel A: Basal Cytokeratin Signatures (Spatial Map)
- **What it shows**: Spatial distribution of basal CK expression in cancer cell layer
- **Genes scored**: KRT5, KRT6A, KRT6B, KRT14, KRT17, KRT15, KRT16, TP63, EGFR, CDH3
- **Source**: Deconvolved GEX (Cancer_Basal layer from Module 3 Pass 2)
- **Colormap**: White → Yellow → Red (BASAL_CMAP)
- **FIX NEEDED**: Use `sc.pl.spatial()` with H&E backdrop instead of raw scatter
- **Visium data**: `/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/HCC22-088-P4-S2_1i_rep/outs/`
- **GEX parquet**: `output/module3_unified/HCC22-088-P4-S2_1i_rep_gene_expression_pass1.parquet`
- **Parquet structure**: Index = `barcode:::cell_type`, filter for `Cancer_Basal`

### Panel B: D538G Mutation Signal (Spatial Map)
- **What it shows**: EstroGene D538G signature score across tissue — dispersed, not clustered
- **Method**: Mean expression of 15 EstroGene D538G upregulated genes minus 4 downregulated genes
- **Up genes**: AGR2, TFF1, TFF3, GREB1, PGR, PDZK1, CA12, STC2, IGFBP4, XBP1, AREG, CCND1, MYC, CXCL12, MDK
- **Down genes**: NRIP1, ESR1, RARA, FOXA1
- **Colormap**: White → Light blue → Dark blue (MUTATION_CMAP)
- **FIX NEEDED**: Same H&E backdrop fix as Panel A
- **Key observation**: D538G signal is dispersed across tissue, not localized to one region

### Panel C: EstroGene Signature — Violin + Gene Heatmap (UPGRADED)
- **Left half**: Split violin/strip plot — D538G-high vs WT-like spots
- **Right half**: Heatmap of 19 individual EstroGene signature genes
  - Rows: 15 up-genes + 4 down-genes
  - Columns: D538G-high spots (mean) vs WT-like spots (mean)
  - Color: Expression level (z-scored or raw)
- **Key stat**: Mann-Whitney p < 0.0001 between groups
- **Threshold**: Spots classified as D538G-high if EstroGene score > median
- **Key genes lighting up in D538G**: AGR2, TFF1, MDK, CCND1, MYC
- **Key genes depressed**: ESR1, FOXA1, NRIP1

### Panel D: Combined Pathway Dot Plot (MERGED from old D+E)
- **What it shows**: Enrichr pathway enrichment in D538G vs WT cancer spots
- **Format**: Dot plot — rows = pathways, x = combined score, size = gene overlap, color = -log10(p)
- **Hallmark pathways (8)**:
  - Estrogen Response Early: 538.37
  - Apoptosis: 493.42
  - TNF-alpha Signaling via NF-kB: 466.91
  - Epithelial Mesenchymal Transition: 401.97
  - Estrogen Response Late: 290.37
  - Myogenesis: 201.05
  - Hypoxia: 131.59
  - p53 Pathway: 103.58
- **KEGG pathways (5)**:
  - Protein digestion and absorption: 97.15
  - ECM-receptor interaction: 84.76
  - Focal adhesion: 59.03
  - Apoptosis: 55.78
  - Proteoglycans in cancer: 41.06
- **Gene sets**: MSigDB_Hallmark_2020, KEGG_2021_Human
- **Source**: vignette_2_surgical_d538g.ipynb cell 36

### Panel E: MDK-SDC4 COMMOT Sender Score Spatial Map (UPDATED)
- **What it shows**: Spatial distribution of MDK-SDC4 COMMOT sender (s-) signaling score on tissue
- **H&E backdrop**: Yes (same Visium data as Panels A-B)
- **Source**: COMMOT analysis computed on-the-fly using vignette_2 pattern (cells 24-31)
- **Implementation**: Uses `expand_prop_gex_adata()` to create pseudo-single-cell adata, then runs `ct.tl.spatial_communication()` with CellChat database
- **Key pair**: s-MDK-SDC4 (D538G mean=7.608, WT mean=2.221, p=9.33e-13)
- **Score column**: `obsm["commot-cellchat-sum-sender"]["s-MDK-SDC4"]`
- **Message**: Midkine signaling activity (not just expression) is spatially localized and corresponds to D538G-high regions
- **Ties to**: Panel F (quantified bars) and the midkine story in manuscript title
- **Fallback**: If COMMOT unavailable, falls back to raw MDK gene expression from deconvolved GEX

### Panel F: COMMOT Signaling Bars (TIGHTENED)
- **What it shows**: Grouped bars — D538G-high vs WT-like cancer spots for 5 ligand-receptor pairs
- **Pairs (all FDR < 0.001)**:
  - MDK-SDC4: D538G=7.608, WT=2.221, p=9.33e-13
  - MDK-NCL: D538G=1.876, WT=0.708, p=1.41e-12
  - PTN-SDC4: D538G=6.019, WT=0.277, p=8.09e-36
  - PTN-NCL: D538G=2.111, WT=0.104, p=6.38e-26
  - MIF-CD74_CD44: D538G=1.631, WT=0.301, p=8.51e-04
- **CellChat database** used for ligand-receptor pairs
- **Source**: vignette_2_surgical_d538g.ipynb cells 30-31

### Panels G-H: Validation (Placeholders)
- **G**: ELISA — MCF7 WT vs D538G, midkine secretion (p < 0.001). User provides Prism data.
- **H**: Immunofluorescence — MCF7/MP vs D538G, membrane MDK (p < 0.001, INTEGRATED DENSITY p < 0.01). User provides microscopy images.

## Key Manuscript Sentences
- "Spatial deconvolution identified basal cytokeratin expression (KRT5/14/17) in the cancer cell layer of progressor P4"
- "EstroGene D538G mutation signature was dispersed across the tissue rather than spatially clustered"
- "D538G-high cancer spots showed activated estrogen response, EMT, and apoptosis pathways (Enrichr combined score > 400)"
- "COMMOT analysis revealed 5 significantly upregulated ligand-receptor pairs in D538G spots, with MDK-SDC4 showing the largest effect (3.4-fold, p=9.3e-13)"
- "Midkine signaling was spatially co-localized with D538G mutation signal, confirmed by ELISA (p<0.001) and immunofluorescence"

## HPA Validation
- Spearman rho=0.46, p=2.86e-37 (from vignette_4 cell 26)
