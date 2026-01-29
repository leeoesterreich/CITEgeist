# Design: Attractor State Evidence Figures

**Date:** 2026-01-28
**Scripts:** `14_cistrome_divergence.py`, `15_proteostasis_setpoint.py`, `16_functional_consequence.py`
**Location:** `midkine/mdk_saturation_pipeline/scripts/`
**Purpose:** Grant/proposal supplement figures showing computational evidence that ESR1-D538G drives qualitatively different stable states in MCF7 vs T47D, supporting proposed wet lab experiments.

## Narrative Arc

1. **Figure 14 — Different Cistrome:** The same D538G mutation redistributes ER binding to different enhancer sets in MCF7 vs T47D, activating distinct transcriptional programs.
2. **Figure 15 — Different Proteostasis Setpoint:** The secretory module, UPR, and ERAD shift as coordinated blocks in opposite directions between cell lines, with MDK mRNA anti-correlating with secretion.
3. **Figure 16 — Functional Consequence:** Secretory capacity and stress response diverge to determine protein fate, bridging computational evidence to experimental predictions.

## Data Sources

All existing pipeline datasets and outputs (no new downloads):

| Dataset / Output | Type | Figures |
|-----------------|------|---------|
| GSE89888 | RNA-seq, MCF7/T47D WT vs D538G | 14, 15, 16 |
| GSE125117 | ER ChIP-seq, MCF7/T47D WT vs D538G | 14, 16 |
| GSE254216 | ATAC-seq, MCF7/T47D | 14 |
| Script 10 output | `unsupervised_hsp90b1_ranking.csv` (632 opposite-regulation genes) | 15 |
| Script 11 output | `venn_expression_binding.csv` (168 intersection genes) | 14 |
| Script 13 output | `mechanism_gene_signatures.csv`, `mechanism_discrimination_scorecard.csv` | 15, 16 |

## Figure 14: Different Cistrome

**Title:** ESR1-D538G engages distinct regulatory programs in MCF7 versus T47D
**File:** `fig14_cistrome_divergence.pdf/png`
**Layout:** 2 rows x 2 columns (4 panels)

### Panel A — Venn diagram

Regenerated from script 11 data with grant-quality styling. Two-circle Venn: opposite expression (632 genes) intersected with differential ER binding (168 genes). Clean labels, counts in each region.

- Data: `venn_expression_binding.csv`
- Key numbers: 632 total, 168 in intersection

### Panel B — ER binding change heatmap

Heatmap of ER ChIP-seq binding changes at the 168 intersection gene promoters.
- Rows: 168 genes, hierarchically clustered
- Columns: MCF7 (D538G - WT) binding diff, T47D (D538G - WT) binding diff
- Color: blue (binding loss) to white (no change) to red (binding gain)
- Side annotation: correlation_type (same_as_secretion = orange, inverse_repressor = purple)
- Data: `venn_expression_binding.csv` columns `MCF7_binding_diff`, `T47D_binding_diff`

### Panel C — ATAC-seq accessibility at intersection genes

Scatter plot of chromatin accessibility at the 168 gene promoters.
- X-axis: MCF7 ATAC peak count or signal at promoter (±10kb TSS)
- Y-axis: T47D ATAC peak count or signal at promoter (±10kb TSS)
- Points colored by correlation_type
- Diagonal y=x line for reference
- Data: GSE254216 ATAC-seq peaks, gene coordinates from mygene

### Panel D — GO functional enrichment

Two dot plots (or horizontal bar charts) side by side:
- Left: Top 5 GO Biological Process terms for "same_as_secretion" genes (n~100)
- Right: Top 5 GO terms for "inverse_repressor" genes (n~68)
- Dot size = gene count, color = -log10(p-value)
- Method: `scipy.stats.fisher_exact` or `statsmodels` with BH correction against full genome background
- Gene lists from `venn_expression_binding.csv` filtered by `correlation_type`

Note: GO enrichment requires gene-to-GO mapping. Use `mygene` API to fetch GO annotations, or `gprofiler` if available. If neither is available in the environment, fall back to manual pathway annotation using KEGG categories.

## Figure 15: Different Proteostasis Setpoint

**Title:** MCF7 and T47D occupy distinct secretory/stress states that diverge under ESR1-D538G
**File:** `fig15_proteostasis_setpoint.pdf/png`
**Layout:** 2 rows x 3 columns (6 panels)

### Panel A — Secretory module co-regulation

Dot plot showing fold-change of each secretory chaperone.
- X-axis: gene name (CALR, CANX, PDIA4, PDIA6, ERO1A, SEC61A1, SRP54, HSP90B1)
- Y-axis: log2(FC) D538G/WT
- Two series: MCF7 (blue dots) and T47D (red dots)
- HSP90B1 highlighted (gold, larger)
- Horizontal line at y=0
- Annotation: "93% concordant, p=0.0009" from Test 3
- Data: `mechanism_gene_signatures.csv`

### Panel B — UPR marker direction

Bar chart of UPR gene fold-changes in MCF7-D538G.
- Genes: XBP1, ATF6, ATF4, DDIT3, ERN1
- Y-axis: log2FC
- Color: red if log2FC > 0.15 (UP), blue if < -0.15 (DOWN), gray otherwise (NC)
- Annotation: "3/5 UP, 0 DOWN, p=0.04" from Test 2
- Data: `mechanism_gene_signatures.csv`

### Panel C — ERAD gene heatmap

Heatmap of ERAD gene fold-changes in both lines.
- Rows: EDEM1, EDEM2, EDEM3, OS9, SYVN1, SEL1L, HERPUD1, DERL1
- Columns: MCF7 log2FC, T47D log2FC
- Color: RdBu_r centered at 0
- Annotated with FC values
- Data: `mechanism_gene_signatures.csv`

### Panel D — MDK mRNA anti-correlation

Grouped bar chart showing MDK log2FC.
- Two bars: MCF7 (-0.45) and T47D (+0.48)
- Error bars from replicate SEM
- Significance brackets with p-values (0.0008 and 0.002)
- Arrows or text annotations: "secretion +83%" above MCF7 bar, "secretion -62%" above T47D bar
- Color: bars colored by direction (blue=down, red=up)
- Data: GSE89888 (individual replicate values for error bars)

### Panel E — Genome-wide opposite regulation

Scatter or MA-style plot of 632 opposite-regulation genes.
- X-axis: MCF7 log2FC
- Y-axis: T47D log2FC
- Points should cluster in quadrants II (MCF7 UP, T47D DOWN) and IV (MCF7 DOWN, T47D UP)
- Color by functional category: secretory module genes (orange), UPR genes (green), ERAD genes (purple), all other (gray)
- HSP90B1 and MDK labeled with text
- Data: `unsupervised_hsp90b1_ranking.csv` + `mechanism_gene_signatures.csv` for category labels

### Panel F — Chaperone-to-client ratio

Bar chart showing secretory chaperone capacity relative to client load.
- 4 bars: MCF7 WT, MCF7 D538G, T47D WT, T47D D538G
- Y-axis: mean TPM of 7 secretory chaperones (CALR, CANX, PDIA4, PDIA6, ERO1A, SEC61A1, SRP54)
- Error bars from gene-level variance across the 7 genes
- Alternative: ratio of mean chaperone TPM to mean UPR marker TPM (capacity / stress proxy)
- Data: GSE89888 via `mechanism_gene_signatures.csv`

## Figure 16: Functional Consequence

**Title:** Secretory capacity, UPR tone, and stress response diverge to determine protein fate
**File:** `fig16_functional_consequence.pdf/png`
**Layout:** 2 rows x 2 columns (4 panels)

### Panel A — Secretory module vs UPR trajectories

Paired dot/slope chart.
- Two metrics per cell line: mean secretory module log2FC and mean UPR log2FC
- X-axis: module (Secretory, UPR)
- Y-axis: mean log2FC
- Two lines connecting the dots: MCF7 (blue) and T47D (red)
- Shows MCF7 secretory UP with UPR UP (expanding under load), T47D secretory DOWN with UPR flat/mixed (contracting)
- Data: `mechanism_gene_signatures.csv`

### Panel B — Three-module heatmap

Combined heatmap of all 20 signature genes.
- Rows: 8 ERAD + 5 UPR + 7 secretory, grouped by module with horizontal dividers
- Columns: MCF7 log2FC, T47D log2FC
- Color: RdBu_r centered at 0
- Side annotation bar: module membership (ERAD=purple, UPR=green, Secretory=orange)
- Annotated with FC values
- No clustering — preserve module grouping
- Data: `mechanism_gene_signatures.csv`

### Panel C — MDK mRNA vs secretory module direction

Scatter plot with 4 points.
- X-axis: mean secretory chaperone log2FC for each condition
- Y-axis: MDK log2FC for each condition
- Points: MCF7 WT (reference, 0,0), MCF7 D538G, T47D WT (reference, 0,0), T47D D538G
- Since WT is the reference (log2FC=0), we have 2 informative points plus the origin
- Annotate each point with cell line + genotype
- Draw trend line if meaningful
- Alternative: use absolute TPM values across all 4 conditions (not fold-changes) to get 4 real data points
- Data: GSE89888

### Panel D — Integrated evidence table

Formatted table (rendered as matplotlib table or text).
- Rows: 3 "why's"
  1. Different cistrome (ER binds different loci)
  2. Different proteostasis setpoint (secretory vs stress modules)
  3. Module coordination determines protein fate
- Columns: Key finding | Dataset | Key statistic | Experimental prediction
- Row 1: "168/632 opposite genes have diff ER binding" | GSE125117+GSE89888 | "168 genes, p<0.05 (Fisher)" | "ChIP-seq in additional ER+ lines will show line-specific binding"
- Row 2: "93% secretory co-regulation, UPR activation" | GSE89888 | "p=0.0009 (binomial)" | "CHX chase: T47D-D538G shorter MDK half-life"
- Row 3: "MDK mRNA anti-correlates with secretion" | GSE89888 | "MCF7 p=0.0008, T47D p=0.002" | "Pulse-chase: MCF7-D538G higher fractional secretion"
- Color-coded by "why" (matching figure colors)

## Statistical Methods

- **GO enrichment:** Fisher's exact test per term, BH FDR correction. Background = all genes in GSE89888 dataset.
- **ATAC signal:** Peak count within ±10kb of TSS from GSE254216. Same approach as script 13 Test 7.
- **Error bars:** SEM across 4 biological replicates (GSE89888) where applicable.
- **All fold-changes:** log2(D538G/WT) from GSE89888.

## Output Files

| File | Description |
|------|-------------|
| `figures/fig14_cistrome_divergence.pdf` | Figure 14 (vector) |
| `figures/fig14_cistrome_divergence.png` | Figure 14 (raster, 300 dpi) |
| `figures/fig15_proteostasis_setpoint.pdf` | Figure 15 (vector) |
| `figures/fig15_proteostasis_setpoint.png` | Figure 15 (raster, 300 dpi) |
| `figures/fig16_functional_consequence.pdf` | Figure 16 (vector) |
| `figures/fig16_functional_consequence.png` | Figure 16 (raster, 300 dpi) |
| `tables/cistrome_intersection_go_enrichment.csv` | GO enrichment for Fig 14D |
| `tables/proteostasis_module_summary.csv` | Module-level summary stats for Fig 15/16 |

## Implementation Structure

Three standalone scripts following existing pipeline conventions:

```python
# Each script follows the pattern:
# - Load config.yaml
# - Load required data (GSE89888 TPM, pipeline outputs)
# - Compute panel-specific metrics
# - Generate multi-panel figure
# - Save outputs

# 14_cistrome_divergence.py
#   load_venn_data()          - Read script 11 output
#   load_atac_peaks()         - GSE254216 peaks
#   get_go_enrichment()       - mygene GO annotations + Fisher test
#   plot_figure_14()          - 4-panel figure

# 15_proteostasis_setpoint.py
#   load_signature_data()     - Read script 13 gene signatures
#   load_ranking_data()       - Read script 10 genome-wide ranking
#   load_replicate_data()     - GSE89888 for error bars
#   plot_figure_15()          - 6-panel figure

# 16_functional_consequence.py
#   load_signature_data()     - Read script 13 gene signatures
#   load_chipseq_data()       - GSE125117 for panel context
#   compute_module_means()    - Aggregate secretory/UPR/ERAD modules
#   plot_figure_16()          - 4-panel figure
```

## Style Requirements (Grant Quality)

- Font: Arial/Helvetica, 10pt minimum for axis labels, 8pt for annotations
- Colors: Consistent palette across all 3 figures
  - MCF7: steelblue
  - T47D: indianred
  - Secretory module: orange
  - UPR module: green
  - ERAD module: purple
  - HSP90B1: gold
- All panels labeled A, B, C, D (etc.) in bold 14pt
- Figure dimensions: each figure ~16x10 inches for high-resolution printing
- White background, minimal gridlines
- Significance annotations: brackets with p-values, not just stars

## Pipeline Integration

Add all three scripts to `run_pipeline.py` SCRIPTS list after script 13.
Add SLURM submission scripts to `slurm/`.
