# Figure 3 Briefing: Benchmarking

## Story
CITEgeist performs automatic cell type profile discovery (Modules 1-2) and competitive deconvolution of both proportions and gene expression (Module 3) across real and simulated benchmarks.

## Panel-by-Panel

### Panel A: Module 1+2 Profile Discovery Accuracy
- **What it shows**: CITEgeist auto-discovers cell type protein marker profiles without prior knowledge, matching known profiles
- **Dataset**: Xenium RCC (Renal Cell Carcinoma) pseudo-Visium
- **Cell types discovered**: B cells (CD20), CD4+ T cells (CD3E/CD4), CD8+ T cells (CD3E/CD8A), Macrophages (CD68/CD163), Endothelial (CD31), Epithelial (PanCK), Fibroblasts (alphaSMA)
- **Key stat**: 6/7 exact marker matches, 1 superset
- **Format**: Compact table or heatmap + deconvolution schematic inset (top-left)
- **Data**: Xenium ground truth from `Benchmarking/xenium_benchmarking/`

### Panel B: Proportion Benchmark — Xenium (Real Data)
- **What it shows**: 5-method comparison on real Xenium RCC data, 5 pseudo-Visium regions, 10 cell types
- **Methods**: CITEgeist, Cell2Location, RCTD, Tangram, Seurat
- **Metrics**: Pearson r, 1-JSD, 1-RMSE (higher = better for all)
- **Key results**:
  - CITEgeist: r=0.600+/-0.061, JSD=0.355+/-0.011, RMSE=0.167+/-0.006
  - Cell2Location: r=0.612+/-0.050, JSD=0.335+/-0.032
  - RCTD: r=0.619+/-0.030, JSD=0.347+/-0.014
  - Tangram: r=0.143+/-0.094 (poor)
  - Seurat: r=0.173+/-0.076 (poor)
- **Message**: CITEgeist is competitive with dedicated proportion methods (within 0.02 of best) while also doing GEX deconvolution
- **Data**: `Benchmarking/xenium_benchmarking/evaluation/results/method_summary.csv`

### Panel C: Proportion Benchmark — Simulated (high_seg + mixed)
- **What it shows**: Same 5 methods on scCube-generated synthetic data under two conditions
- **Conditions**: "high_seg" (homogeneous regions) and "mixed" (heterogeneous tissue)
- **Basis**: Wu et al. 2021 breast cancer reference, 9 cell types, 5 replicates each
- **Cell types**: B-cells, CAFs, Cancer Epithelial, Endothelial, Myeloid, Normal Epithelial, PVL, Plasmablasts, T-cells
- **Format**: Faceted bars or heatmap, two conditions side by side
- **Data**: `Benchmarking/simulation_benchmarking/Figures/prop_all_metrics_{highseg,mixed}_combined.csv`

### Panel D: GEX Benchmark — Simulated (CITEgeist Dominance)
- **What it shows**: Gene expression deconvolution accuracy — CITEgeist's biggest competitive advantage
- **Methods**: CITEgeist, Cell2Location, Tangram (3 methods with GEX output)
- **Key results (high_seg)**:
  - CITEgeist: RMSE=0.293, NRMSE=0.039
  - Cell2Location: RMSE=0.493, NRMSE=0.066
  - Tangram: RMSE=1.269, NRMSE=0.172
- **Key results (mixed)**:
  - CITEgeist: RMSE=0.257 (BEST)
  - Cell2Location: RMSE=1.303 (5x worse)
  - Tangram: RMSE=1.436 (5.6x worse)
- **Message**: CITEgeist is 2-5x better at GEX deconvolution than competitors, especially on heterogeneous tissue
- **Data**: `Benchmarking/simulation_benchmarking/Figures/all_GEX_metrics_{high_seg,mixed}.csv`

### Panel E: Predicted vs Ground Truth Scatter
- **What it shows**: Per-spot predicted vs actual proportions for one Xenium region
- **Format**: Scatter plot, color by cell type, regression line, r value annotation
- **Data**: `Benchmarking/xenium_benchmarking/xenium_pseudovisium/data/ground_truth/`
- **NOTE**: Current panel is BROKEN (Errno 2). Must fix file path.

### Panel F: Spatial Comparison Maps
- **What it shows**: Side-by-side spatial maps — ground truth vs CITEgeist predictions for 2-3 cell types on one Xenium region
- **Format**: Multi-panel spatial maps with consistent colorscale
- **Cell types to show**: Epithelial, Macrophages, T cells (high visual contrast)
- **Message**: Spatial patterns are faithfully recapitulated

## Key Manuscript Sentences
- "CITEgeist automatically discovered protein marker profiles matching 6 of 7 known cell type markers on Xenium RCC data"
- "Proportion deconvolution was competitive with state-of-the-art methods (Pearson r = 0.60 vs RCTD 0.62, Cell2Location 0.61)"
- "Gene expression deconvolution outperformed all methods by 2-5x on RMSE (0.257 vs Cell2Location 1.303 on mixed tissue)"
- "Performance was robust across both homogeneous (high_seg) and heterogeneous (mixed) tissue architectures"
