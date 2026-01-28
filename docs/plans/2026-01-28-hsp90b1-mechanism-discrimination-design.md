# Design: HSP90B1 Mechanism Discrimination Script

**Date:** 2026-01-28
**Script:** `13_discriminate_hsp90b1_mechanism.py`
**Location:** `midkine/mdk_saturation_pipeline/scripts/`

## Problem

The existing MDK saturation pipeline (scripts 01-12) establishes that HSP90B1 is differentially regulated between MCF7 and T47D via ER chromatin saturation. What remains unresolved is the downstream mechanism: **how** does changing HSP90B1 abundance translate to changes in MDK secretion?

Two mechanistic models are plausible:

- **Model A (Direct folding/export):** HSP90B1 as ER chaperone improves MDK folding capacity, reducing ER retention and ERAD. Higher HSP90B1 directly enables more MDK to fold correctly and exit the ER.
- **Model B (LRP1 premature binding trap):** HSP90B1 matures LRP1 (an MDK receptor), which can trap MDK intracellularly during biosynthesis. Two sub-models exist: "rescue" (better chaperoning reduces aberrant MDK-LRP1 aggregation) and "exposure" (more LRP1 flux increases trapping opportunities).

## Approach

Hypothesis-driven gene signature analysis using existing pipeline datasets. Eight diagnostic tests produce an evidence scorecard scoring each model, with effect sizes and p-values. No forced winner — the reader interprets the tally.

## Data Sources

All existing pipeline datasets (no new downloads):

| Dataset | Type | Use |
|---------|------|-----|
| GSE89888 | RNA-seq, MCF7/T47D WT vs D538G | Primary expression signatures (all 8 tests) |
| GSE75329 | RNA-seq, MCF7 siER/siFOXA1/FOXA1-OE | Perturbation validation |
| GSE125117 | ER ChIP-seq, MCF7/T47D WT vs D538G | ER binding at LRP1 locus (Test 7) |
| GSE254216 | ATAC-seq, MCF7/T47D | Chromatin accessibility at LRP1 (Test 7) |

## Diagnostic Tests

### Model A (Direct folding/export)

**Test 1 — ERAD pathway suppression.**
When HSP90B1 goes UP (MCF7-D538G), ERAD genes should go DOWN or stay flat (less misfolding = less degradation needed). If ERAD goes UP with HSP90B1, that undermines the "better folding" narrative.

- Genes: EDEM1, EDEM2, EDEM3, OS9, SYVN1, SEL1L, HERPUD1, DERL1
- Score: mean fold-change of ERAD genes, one-sample t-test against zero
- Verdict: mean FC < 0 when HSP90B1 UP → Supports A; mean FC > 0 → Supports B; p > 0.05 → Ambiguous

**Test 2 — UPR stress relief.**
UPR sensors/effectors should NOT activate when HSP90B1 increases. If UPR markers decrease in MCF7-D538G, that supports Model A.

- Genes: XBP1, ATF6, ATF4, DDIT3, ERN1
- Score: directional concordance with "stress relief" prediction
- Verdict: majority of UPR genes decrease → Supports A; increase → Supports B; mixed → Ambiguous

**Test 3 — Secretory capacity co-regulation.**
Other ER folding chaperones and translocon components should co-regulate with HSP90B1 if this is a general secretory capacity shift.

- Genes: CALR, CANX, PDIA4, PDIA6, ERO1A, SEC61A1, SRP54
- Score: fraction of secretory genes changing in the same direction as HSP90B1, per cell line
- Verdict: fraction > 0.6 → Supports A; fraction < 0.4 → Ambiguous; n.s. → Ambiguous

**Test 4 — Known GRP94 client concordance.**
GRP94 has established clients. If these show expression patterns consistent with improved maturation when HSP90B1 is up, that supports the direct chaperoning model.

- Genes: IGF1, IGF2, ITGA4, ITGAL, TLR2, TLR4, TLR9
- Score: client gene fold-change correlation with HSP90B1 fold-change
- Verdict: positive correlation (r > 0.3) → Supports A; no correlation → Ambiguous; negative → Supports B

### Model B (LRP1 premature binding trap)

**Test 5 — LRP1 expression and regulation.**
Is LRP1 meaningfully expressed in MCF7/T47D? Does it change with D538G? If LRP1 is absent or very low, Model B is implausible.

- Genes: LRP1
- Score: LRP1 TPM and fold-change, minimum expression threshold (TPM > 1)
- Verdict: LRP1 not expressed → Supports A (by elimination); LRP1 expressed and regulated → Supports B; expressed but unchanged → Ambiguous

**Test 6 — RAP/LRPAP1 gating.**
LRPAP1 (RAP) prevents premature ligand-receptor binding in the secretory pathway. If LRPAP1 is highly expressed and unchanged, premature MDK-LRP1 binding is suppressed and Model B is unlikely.

- Genes: LRPAP1
- Score: LRPAP1 expression level and D538G fold-change
- Verdict: LRPAP1 high and stable → Supports A (trap suppressed); LRPAP1 low or decreasing → Supports B (trap open); intermediate → Ambiguous

**Test 7 — LRP1-HSP90B1 co-regulation.**
Model B predicts LRP1 and HSP90B1 should change together (rescue sub-model) or that their interaction determines MDK output.

- Data: GSE89888 expression + GSE125117 ER ChIP-seq + GSE254216 ATAC-seq
- Score: correlation of LRP1 and HSP90B1 across conditions; ER binding at LRP1 locus (±10kb TSS)
- Verdict: strong co-regulation + ER binding at LRP1 → Supports B; no co-regulation, no binding → Supports A; mixed → Ambiguous

### Consistency Check

**Test 8 — Opposite-model consistency.**
For each cell line separately, check whether the full pattern is internally consistent with one model. MCF7-D538G (HSP90B1 UP, MDK secretion UP) and T47D-D538G (HSP90B1 DOWN, MDK secretion DOWN) should both be explained by the same model without contradiction.

- Input: results from Tests 1-7
- Score: does the winning model explain both cell lines?
- Verdict: consistent across both lines → reinforces leading model; inconsistent → flags that the mechanism may be context-dependent

## Figure Design

**File:** `fig13_mechanism_discrimination.pdf/png`
**Layout:** 3 rows × 3 columns (9 panels)

| Panel | Content |
|-------|---------|
| A | ERAD gene heatmap: fold-changes (D538G/WT) for MCF7 and T47D, annotated with direction arrows |
| B | UPR marker bar chart: fold-changes per cell line, colored by stress-relief vs stress-activation |
| C | Secretory capacity dot plot: HSP90B1 vs other ER chaperones, showing co-regulation fraction |
| D | GRP94 client scatter: client gene fold-changes vs HSP90B1 fold-change, per cell line |
| E | LRP1/LRPAP1 expression bars: TPM values in WT and D538G for both lines, with significance brackets |
| F | LRP1 locus epigenomic strip: ER ChIP-seq signal + ATAC-seq accessibility at LRP1 promoter region |
| G | Correlation matrix: HSP90B1, MDK, LRP1, LRPAP1 pairwise correlations across all conditions |
| H | Cell-line consistency diagram: flow chart showing which model explains both MCF7 and T47D without contradiction |
| I | Scorecard summary: 8 tests as rows, columns for verdict/effect size/p-value, colored blue (Model A) / orange (Model B) / gray (Ambiguous) |

## Statistical Methods

- **Expression fold-changes:** log2(D538G/WT) from GSE89888, per-gene Welch's t-test across 4 replicates
- **Gene set tests:** One-sample t-test on mean fold-change of gene set against zero (Tests 1-3)
- **Correlation:** Pearson correlation of fold-changes across cell lines/conditions (Tests 4, 7)
- **Expression thresholds:** Gene considered "expressed" if mean TPM > 1 across replicates (Tests 5, 6)
- **ChIP-seq binding:** Peak presence within ±10kb of TSS from GSE125117 (Test 7)
- **ATAC-seq accessibility:** Normalized signal at gene promoter from GSE254216 (Test 7)
- **Multiple testing:** Benjamini-Hochberg FDR within scorecard; individual tests report raw p-values
- **Effect sizes:** Cohen's d for expression comparisons, Pearson r for correlations

## Output Files

| File | Description |
|------|-------------|
| `tables/mechanism_discrimination_scorecard.csv` | Per-test verdicts with effect sizes and p-values |
| `tables/mechanism_gene_signatures.csv` | Expression data for all diagnostic genes across conditions |
| `figures/fig13_mechanism_discrimination.pdf` | Multi-panel figure (vector) |
| `figures/fig13_mechanism_discrimination.png` | Multi-panel figure (raster) |

## Implementation Structure

```python
@dataclass
class TestResult:
    test_id: int
    test_name: str
    model_tested: str          # "A", "B", or "both"
    verdict: str               # "Supports A", "Supports B", "Ambiguous"
    effect_size: float         # Cohen's d or Pearson r
    p_value: float
    details: str               # One-sentence interpretation
    gene_level_data: pd.DataFrame

# Functions:
# load_config()              - Read config.yaml
# load_expression_data()     - GSE89888 RNA-seq
# load_perturbation_data()   - GSE75329 perturbations
# load_chipseq_data()        - GSE125117 ER ChIP-seq
# load_atac_data()           - GSE254216 ATAC-seq
# define_gene_signatures()   - Returns dict of gene lists per model
# run_test_1_erad()          - ERAD suppression
# run_test_2_upr()           - UPR stress relief
# run_test_3_secretory()     - Secretory capacity co-regulation
# run_test_4_clients()       - GRP94 client concordance
# run_test_5_lrp1()          - LRP1 expression/regulation
# run_test_6_rap()           - RAP/LRPAP1 gating
# run_test_7_coregulation()  - LRP1-HSP90B1 co-regulation + epigenomics
# run_test_8_consistency()   - Cross-cell-line consistency
# compile_scorecard()        - Aggregate TestResults into CSV
# plot_figure()              - 3x3 multi-panel figure
# main()                     - Orchestrator
```

## Pipeline Integration

- Add script 13 to `run_pipeline.py` supplementary script list (after script 10)
- Add to `config.yaml`:
  - Gene coordinates for LRP1 and LRPAP1 (hg38 and hg19)
  - Gene signature lists under a `mechanism_discrimination` key
- No new data downloads required
- No new Python package dependencies

## Gene Signature Definitions

```yaml
mechanism_discrimination:
  erad_genes: [EDEM1, EDEM2, EDEM3, OS9, SYVN1, SEL1L, HERPUD1, DERL1]
  upr_genes: [XBP1, ATF6, ATF4, DDIT3, ERN1]
  secretory_chaperones: [CALR, CANX, PDIA4, PDIA6, ERO1A, SEC61A1, SRP54]
  grp94_clients: [IGF1, IGF2, ITGA4, ITGAL, TLR2, TLR4, TLR9]
  lrp1_axis: [LRP1, LRPAP1]

gene_coordinates_hg38:
  LRP1: ["chr12", 57170000, 57255000]
  LRPAP1: ["chr4", 3410000, 3440000]

gene_coordinates_hg19:
  LRP1: ["chr12", 57500000, 57585000]
  LRPAP1: ["chr4", 3440000, 3470000]
```
