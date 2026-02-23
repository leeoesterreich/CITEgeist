# Figure Redesign Plan — 2026-02-09

## Motivation

Figures 3, 4, and 5 have excessive whitespace, low information density, and several broken or placeholder panels. The redesign packs real data into every panel, emphasizes spatial biology unique to CITEgeist, and tells a compelling biological story about treatment response and resistance.

## Design Principles

- Zero empty schematics as full panels (shrink to insets where needed)
- Every panel should show data that requires spatial multi-omics
- Spatial maps use `sc.pl.spatial()` with H&E backdrop (not raw scatter)
- Bulk-analyzable results (volcano, pathway bars) move to supplementary
- Bivariate spatial relationships are the unique selling point

---

## Figure 3: Benchmarking (6 panels, 12x14)

### Narrative
CITEgeist automatically discovers cell type profiles (Modules 1-2) and deconvolves both proportions and gene expression competitively against state-of-the-art methods on real (Xenium) and simulated data.

### Layout: 3 rows x 2 columns

| Panel | Content | Size |
|-------|---------|------|
| A | Module 1+2 profile discovery accuracy + deconvolution schematic inset | Top-left |
| B | Proportion benchmark — Xenium real data (5 methods, 3 metrics) | Top-right |
| C | Proportion benchmark — Simulated (high_seg + mixed, 5 methods) | Mid-left |
| D | GEX benchmark — Simulated (CITEgeist dominance, both conditions) | Mid-right |
| E | Scatter — predicted vs ground truth proportions (Xenium region) | Bot-left |
| F | Spatial comparison maps (GT vs CITEgeist predictions, 2-3 cell types) | Bot-right |

### Data Sources
- Xenium proportions: `Benchmarking/xenium_benchmarking/evaluation/results/method_summary.csv`
- Xenium per-region: `Benchmarking/xenium_benchmarking/evaluation/results/comparison_table.csv`
- Simulated proportions: `Benchmarking/simulation_benchmarking/Figures/prop_all_metrics_{highseg,mixed}_combined.csv`
- Simulated GEX: `Benchmarking/simulation_benchmarking/Figures/all_GEX_metrics_{high_seg,mixed}.csv`
- Xenium ground truth: `Benchmarking/xenium_benchmarking/xenium_pseudovisium/data/ground_truth/`

### Key Changes from Current
- Panel A: Shrink schematic to inset, add Module 1+2 discovery accuracy table
- Panel C: NEW — simulated proportion benchmark (was not shown)
- Panel D: NEW — simulated GEX benchmark (CITEgeist's biggest strength)
- Panel E: FIX broken scatter plot (Errno 2 file not found)
- Panel F: NEW — spatial comparison maps for visual validation

---

## Figure 4: Midkine/ESR1 Case Study (8 panels, 12x16)

### Narrative
CITEgeist reveals upregulated midkine signaling in ESR1 D538G-mutant breast cancer (P4-S2_1i_rep). Spatial deconvolution identifies basal CK signatures, mutation distribution, pathway activation, and ligand-receptor communication — validated by ELISA and immunofluorescence.

### Layout: 4 rows x 2 columns

| Panel | Content | Change |
|-------|---------|--------|
| A | Basal CK spatial map (H&E backdrop) | FIX: switch to sc.pl.spatial() |
| B | D538G mutation spatial map (H&E backdrop) | FIX: switch to sc.pl.spatial() |
| C | EstroGene violin + gene-level heatmap | UPGRADE: split panel, add 19-gene heatmap |
| D | Combined pathway dot plot (Hallmark + KEGG) | MERGE: was two separate bar panels |
| E | MDK signaling spatial map (H&E backdrop) | NEW: freed from D+E merge |
| F | COMMOT bars (tightened) | TIGHTEN layout |
| G | ELISA validation | Placeholder (user handles) |
| H | IF validation | Placeholder (user handles) |

### Key Fixes
1. **Panels A-B**: Replace `ax.scatter()` with `sc.pl.spatial()` + `sq.read.visium(..., load_images=True)` for H&E backdrop
2. **Panel C**: Add gene-level heatmap (15 up + 4 down EstroGene genes)
3. **Panels D+E merged → Panel D**: Dot plot (13 pathways, size=overlap, color=-log10 p)
4. **New Panel E**: Spatial map of MDK-SDC4 signaling intensity on tissue

### Data Sources
- Visium data: `/ix1/alee/LO_LAB/General/Lab_Data/20250210_CITEGeistPublicData_GEO_Alex/processed_files/HCC22-088-P4-S2_1i_rep/outs/`
- Proportions: `output/module3_unified/HCC22-088-P4-S2_1i_rep_cell_prop_finetuned_results.csv`
- Deconvolved GEX: `output/module3_unified/HCC22-088-P4-S2_1i_rep_gene_expression_pass1.parquet`
- Module 4: `output/module4_unified/HCC22-088-P4-S2_1i_rep_module4_discovery.json`
- Pathway/COMMOT values: hardcoded from vignette_2 analysis (cells 30-31, 36)

---

## Figure 5: Cross-Sample Integration + Response Biology (8 panels, 12x16)

### Narrative
CITEgeist discovers spatially coherent gene programs across 12 samples from 6 patients, reveals conserved spatial relationships, and identifies the key biological discriminator: responders show T cell-fibroblast co-localization (coordinated remodeling), progressors show T cell exclusion from cancer (immune evasion). Every panel requires spatial data — nothing here could be done with bulk or single-cell.

### Layout: 4 rows x 2 columns

| Panel | Content | Type |
|-------|---------|------|
| A | Cell type proportion shift (resp vs prog, biopsy→surgical) | Grouped bar + dots |
| B | Program conservation heatmap (65 aligned x 12 samples) | Heatmap |
| C | Response enrichment dot plot (65 programs) | Dot plot |
| D | Conserved relationship network (65 nodes, 163 edges) | Network graph |
| E | Bivariate spatial: co-localization (P3-S2 responder) | Tissue maps + scatter |
| F | Bivariate spatial: exclusion (P1-S2 progressor) | Tissue maps + scatter |
| G | Moran's I validation (violin + strip by cell type) | Violin plot |
| H | Spatial proportion maps (P3-S2 vs P1-S2 fibroblasts/cancer) | Tissue maps |

### Supplementary (moved from main)
- **Supp Fig 2A**: Volcano plot (203 sig genes, 120 resp-up, 83 prog-up)
- **Supp Fig 2B**: Pathway enrichment (Hallmark + GO, resp vs prog split)

### Data Sources
- Aligned programs: `output/module5_integration/module5_unified_aligned_programs.csv`
- Conserved relationships: `output/module5_integration/module5_unified_conserved_relationships.csv`
- Response analysis: `output/module5_integration/module5_response_analysis.json`
- Module 5 summary: `output/module5_integration/module5_summary.json`
- Network: `output/module5_integration/module5_unified_similarity_network.graphml`
- PyDESeq2: `examples/output_module5_pydeseq/pseudobulk_de_results.csv`
- Cell proportions: `output/module3_unified/HCC22-088-P*-S*_cell_prop_finetuned_results.csv`
- Module 4 programs: `output/module4_unified/HCC22-088-P*_module4_discovery.json`
- Visium spatial data: `/ix1/alee/LO_LAB/General/Lab_Data/.../processed_files/HCC22-088-P*-S*/outs/`

### Key Spatial Maps Needed
- **Panel E**: P3-S2 (responder surgical) — Fibroblast ECM program + CD4 stromal program
- **Panel F**: P1-S2 (progressor surgical) — Cancer Luminal program + CD4 T cell program
- **Panel H**: P3-S2 vs P1-S2 — fibroblast proportion + cancer luminal proportion

---

## Implementation Notes

### Spatial Plotting Pattern (for all spatial panels)
```python
import squidpy as sq
import scanpy as sc

adata = sq.read.visium(
    path=f"{DATA_ROOT}/{sample}/outs",
    counts_file='filtered_feature_bc_matrix.h5',
    load_images=True,
    gex_only=False
)
adata.obs['score'] = computed_values
sc.pl.spatial(adata, color='score', cmap=cmap, ax=ax, show=False)
```

### Conda Environment
- Path: `/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env`
- Requires: scanpy, squidpy, matplotlib, networkx, adjustText

### Execution
- All figure scripts should be sbatched (loading Visium images requires memory)
- Scripts: `manuscript/figures/generate_figure{3,4,5}.py`

---

## Timeline
1. Update generate_figure3.py — add simulated benchmarks, fix scatter, add spatial maps
2. Update generate_figure4.py — fix H&E backdrop, merge pathways, add MDK spatial
3. Rewrite generate_figure5.py — new 8-panel layout with bivariate spatial
4. Sbatch all three, export panels
5. Update manuscript text and figure legends
