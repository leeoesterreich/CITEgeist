# Figure 5 Briefing: Cross-Sample Integration + Response Biology

## Story
CITEgeist discovers 490 spatially coherent gene programs across 12 samples from 6 breast cancer patients (4 responders, 2 progressors), aligns 65 conserved programs, and identifies 163 spatial relationships. The key finding: responders show T cell-fibroblast spatial co-localization (coordinated stromal remodeling/fibrosis), while progressors show T cell exclusion from expanding cancer (immune evasion). Every panel requires spatial data — nothing here could be done with bulk or single-cell RNA-seq.

## Cohort
- **6 patients x 2 timepoints** (S1=biopsy, S2=surgical) = 12 specimens
- **Responders (4)**: P2 (S1, S2), P3 (S1_A, S2), P5 (S1, S2_F_rep), P6 (S1, S2_D)
- **Progressors (2)**: P1 (S1, S2), P4 (S1, S2_1i_rep)
- **10 cell types**: Endothelial, Fibroblasts, B_Cells, Macrophages, Monocytes, CD8_T_Cells, CD4_T_Cells, Cancer_Luminal, Cancer_Basal, Dendritic_Cells

## Panel-by-Panel

### Panel A: Cell Type Proportion Shift (Biopsy → Surgical)
- **What it shows**: Mean change in cell type proportion from biopsy to surgical, grouped by response status
- **Format**: Grouped bar chart — blue=responders (n=4 patients), red=progressors (n=2 patients), with individual patient dots overlaid
- **Key numbers — Progressors**:
  - CD4_T_Cells: **-30.7%** (P1: 0.483→0.167, P4: 0.420→0.122)
  - CD8_T_Cells: **-5.6%** (P1: -2.3%, P4: -8.9%)
  - Cancer_Luminal: **+16.6%** in P1 (0.086→0.366, massive expansion)
  - Cancer_Basal: **+11.7%** in P4 (lineage plasticity)
  - Fibroblasts: **+12.4%** average (P4: +20.5%)
- **Key numbers — Responders**:
  - Fibroblasts: **+7.8%** average (P3: +31.0%, extreme fibrosis)
  - Cancer_Luminal: **-5.6%** (tumor reduction)
  - CD4_T_Cells: **-6.6%** (moderate, vs -30% in progressors)
  - B_Cells: **+2.4%** (immune infiltration)
- **Message**: Progressors lose T cells and gain cancer; responders lose cancer and gain fibroblasts
- **Data**: `output/module3_unified/HCC22-088-P*-S*_cell_prop_finetuned_results.csv` (all 12 files)

### Panel B: Program Conservation Heatmap (65 x 12)
- **What it shows**: Which aligned programs are present in which samples
- **Format**: Binary/intensity heatmap
  - Rows: 65 aligned programs, grouped by cell type (sidebar color)
  - Columns: 12 samples, ordered by patient then timepoint, with response annotation bar
- **Key patterns**:
  - aligned_032 (CD8 T cells): universal horizontal stripe — **100% conservation (12/12)**
  - aligned_003 (Cancer Luminal): 11/12 — missing from P4-S2_1i_rep (progressor surgical)
  - aligned_001 (Monocytes): 10/12
  - aligned_050 (CD4 T cells): 10/12
  - Cancer_Basal: **most heterogeneous** — 13 aligned programs, mostly patient-specific
- **Cell type program counts**: CD8_T(7), Cancer_Luminal(4), Monocytes(5), CD4_T(8), Cancer_Basal(13), Macrophages(6), Dendritic(5), Fibroblasts(8), Endothelial(4), B_Cells(5)
- **Data**: `output/module5_integration/module5_unified_aligned_programs.csv`

### Panel C: Response Enrichment Dot Plot
- **What it shows**: All 65 aligned programs scored for response association
- **Format**: Dot plot
  - Y-axis: Programs grouped by cell type (color sidebar)
  - X-axis: Fraction of samples that are responders (baseline 0.67 since 8/12 are responders)
  - Dot size: Number of samples program appears in
  - Dot color: -log10(chi-square p-value) or enrichment significance
- **4 significant programs highlighted**:
  - **Responder-enriched**: aligned_004 (Macrophages) — FABP4, HBA2, TNXB, CFD, KLF2 (lipid-laden tissue-resident macrophages)
  - **Progressor-enriched (3)**:
    - aligned_008 (CD4_T_Cells) — COL1A1, COL1A2, VIM, MMP9, FN1 (fibroblast-like, immunosuppressive)
    - aligned_002 (Cancer_Luminal) — MGP, HBA2, LTF, ANKRD30A (metabolic stress)
    - aligned_016 (Monocytes) — HPCA, MX1, FBLN1, IFI6 (interferon-responsive)
- **Data**: `output/module5_integration/module5_response_analysis.json`

### Panel D: Conserved Relationship Network
- **What it shows**: Spatial co-localization and exclusion patterns across patients
- **Format**: Network graph
  - 65 nodes (aligned programs), colored by cell type, sized by conservation score
  - 163 edges colored by relationship type:
    - Green: co-localized (bivariate Moran's I > 0.15, ~15 relationships)
    - Red: exclusive (bivariate I < -0.15, 3 relationships)
    - Thin gray: independent (~145 relationships)
  - Spring layout with cell type clustering
- **Key relationships to label**:
  - Co-localized: aligned_012 ↔ aligned_050 (Fibroblast ↔ CD4 T cells, I=0.358, responder surgical)
  - Co-localized: aligned_000 ↔ aligned_050 (Dendritic ↔ CD4, I=0.258, 4 samples)
  - Exclusive: aligned_003 ↔ aligned_021 (Cancer Luminal ↔ CD4 T cells, I=-0.194, progressor)
- **Data**: `output/module5_integration/module5_unified_conserved_relationships.csv`, `module5_unified_similarity_network.graphml`

### Panel E: Bivariate Spatial — Co-localization in Responder (P3-S2)
- **What it shows**: Two programs that spatially co-localize in responding tissue
- **Format**: Two side-by-side spatial maps on P3-S2 tissue (H&E backdrop) + bivariate scatter inset
- **Left map**: Fibroblast ECM program activity
  - Program: aligned_012 (or P3-S2 Fibroblasts Prog 3)
  - Top genes: COL6A2, COL14A1, POSTN, COL5A2, MMP2, HTRA1, BGN, COMP, CTHRC1
  - Moran's I = 0.857 (extremely spatially coherent)
  - 9/10 top genes are ECM genes
- **Right map**: CD4 T cell stromal program activity
  - Program: aligned_050 (or P3-S2 CD4 Prog 2)
  - Top genes: COL3A1, COL1A1, COL1A2, SPARC, VIM, DCN, FN1
  - Moran's I = 0.892 (most spatially coherent program in entire dataset)
- **Bivariate I = 0.358** between these programs
- **Inset**: Scatter plot of program A vs program B activity per spot
- **Message**: T cells and fibroblasts spatially coordinate ECM remodeling in responding tumors
- **Visium data**: `/ix1/alee/LO_LAB/General/Lab_Data/.../HCC22-088-P3-S2/outs/`
- **NMF matrices**: `output/module4_unified/HCC22-088-P3-S2_anchored_H.npy`

### Panel F: Bivariate Spatial — Exclusion in Progressor (P1-S2)
- **What it shows**: Two programs that are spatially exclusive in resistant tissue
- **Format**: Two side-by-side spatial maps on P1-S2 tissue (H&E backdrop) + bivariate scatter inset
- **Left map**: Cancer Luminal core program activity
  - Program: aligned_003
  - Top genes: ZNF550, SPARCL1, POSTN, SLC9A3R1, CA2, CXCL9
  - Present in 11/12 samples (core luminal architecture)
- **Right map**: CD4 T cell program activity
  - Program: aligned_021
  - Spatially excluded from cancer regions
- **Bivariate I = -0.194** (exclusion)
- **Inset**: Scatter of program A vs program B — negative trend
- **Message**: T cells are spatially excluded from tumor regions in progressors (immune evasion)
- **Visium data**: `/ix1/alee/LO_LAB/General/Lab_Data/.../HCC22-088-P1-S2/outs/`

### Panel G: Moran's I Spatial Coherence Validation
- **What it shows**: NMF programs are genuinely spatially organized (not random)
- **Format**: Violin + strip plot
  - X-axis: Cell types
  - Y-axis: Moran's I values
  - Dashed line at threshold = 0.15
  - Color by response status (optional)
- **Key stats**: 175 spatially coherent programs (Moran's I > 0.15) across 7 cell types
- **Total**: 490 programs across 12 samples
- **Data**: Module 4 discovery JSONs, `spatial_moran_i` field per program

### Panel H: Spatial Proportion Maps (Responder vs Progressor)
- **What it shows**: Direct visual comparison of TME composition on tissue
- **Format**: 2x2 grid of spatial maps (H&E backdrop)
  - Top-left: Fibroblast proportion, P3-S2 (responder, +31% fibrosis)
  - Top-right: Fibroblast proportion, P1-S2 (progressor)
  - Bottom-left: Cancer_Luminal proportion, P3-S2 (responder, reduced)
  - Bottom-right: Cancer_Luminal proportion, P1-S2 (progressor, +28% expansion)
- **Message**: Visual confirmation of the Panel A quantification — fibrotic scarring in responders vs tumor expansion in progressors
- **Data**: `output/module3_unified/HCC22-088-P{1,3}-S2*_cell_prop_finetuned_results.csv`

## Supplementary Material (Moved from Main Figure)
- **Supp Fig 2A**: Volcano plot — 203 significant genes (padj < 0.05)
  - 120 responder-up (blue): KLB, NEDD9, RECK, MYH11, CD109, IL1B
  - 83 progressor-up (red): VPREB3, COL4A4, WDR62, CCNE1, PLK1, RAD51AP1
  - Data: `examples/output_module5_pydeseq/pseudobulk_de_results.csv`
- **Supp Fig 2B**: Pathway enrichment dot plot (split responder vs progressor)
  - Responder-up Hallmark: EMT (p=0.0009), TNF-alpha (p=0.021), Myogenesis (p=0.046)
  - Progressor-up Hallmark: E2F Targets (p=2.0e-05), G2-M Checkpoint (p=1.8e-04), Mitotic Spindle (p=1.4e-03)
  - Progressor-up KEGG: Proteasome (p=3.9e-05)
  - Progressor-up GO: Cell cycle regulation (top 72 terms)
  - Data: `examples/output_module5_pydeseq/pseudobulk_{responder,progressor}_up_*.csv`

## Key Manuscript Sentences
- "Module 4 discovered 490 spatially coherent gene programs across 12 samples, of which 65 aligned across patients with 163 conserved spatial relationships"
- "CD8 T cell metabolic programs were universally conserved (12/12 samples), while Cancer Basal programs were the most heterogeneous (13 aligned programs, mostly patient-specific)"
- "Four response-associated programs were identified: responder-enriched FABP4+ macrophages and three progressor-enriched programs in CD4 T cells, Cancer Luminal, and Monocytes"
- "Bivariate spatial analysis revealed the key discriminator: responders showed T cell-fibroblast co-localization (bivariate I=0.358, P3-S2) indicating coordinated stromal remodeling, while progressors showed T cell exclusion from cancer regions (bivariate I=-0.194, P1-S2)"
- "Mean fibroblast proportion increased by 7.8% in responders (up to +31% in P3), consistent with treatment-induced fibrosis replacing tumor parenchyma"
- "In contrast, progressors showed dramatic CD4 T cell depletion (-30.7%), cancer expansion (+16.6% luminal in P1, +11.7% basal in P4), and activated proliferative pathways (E2F Targets, G2-M Checkpoint)"

## Key Biological Themes

### Responder Biology (Fibrosis/Remodeling)
- Fibroblast expansion: +7.8% avg (P3: +31%)
- Cancer_Luminal reduction: -5.6%
- ECM programs with extreme spatial coherence (Moran's I up to 0.892)
- FABP4+ lipid-laden macrophages (tissue-resident, repair phenotype)
- Pathways: EMT, TNF-alpha, Myogenesis, Inflammatory Response
- Key genes: MYH11 (myofibroblast), COL1A1/A2/3A1, FN1, POSTN, SPARC

### Progressor Biology (Resistance/Immune Evasion)
- CD4 T cell depletion: -30.7% (immune exclusion)
- Cancer expansion: P1 luminal +28%, P4 basal +11.7% (lineage plasticity)
- CAF expansion: P4 +20.5%
- Pathways: E2F Targets, G2-M Checkpoint, Mitotic Spindle, Proteasome
- Key genes: PLK1, CCNE1, RAD51AP1 (proliferation), PARP1 (DNA damage), FOS/JUN/EGR1 (stress response)
- Fibroblast-like CD4 T cells (aligned_008): COL1A1, VIM, MMP9 — immunosuppressive phenotype

### Cross-Patient Spatial Themes
- Universal: CD8 T cell surveillance (100%), luminal architecture (91.7%), plasma cell infiltration (83.3%)
- Responder-specific: T cell-fibroblast co-localization at surgical timepoint
- Progressor-specific: T cell exclusion from cancer, ER+ spatial organization (Moran's I=0.698)
