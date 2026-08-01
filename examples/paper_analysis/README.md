# Paper analysis code

This directory contains the analysis scripts that produced the numbers and figure panels
reported in the CITEgeist manuscript. It is a **curated** set, not a dump of the authors'
working tree: every file below either produced a value that appears in the paper, generated
a main-text figure sub-panel, or is imported by something that does.

If you are looking for how to *use* CITEgeist, start at [`../README.md`](../README.md) and
[`../../docs/quickstart_real_visium.md`](../../docs/quickstart_real_visium.md) instead. This
directory is for checking the paper.

---

## Data

All input datasets are public.

| Dataset | Used for | Source |
|---|---|---|
| Visium + CITE-seq clinical-trial samples (this study, 12 specimens from 6 patients) | Figures 5, 6; the competitor-GEX validation | NCBI GEO **GSE289326** |
| Wu et al. 2021 breast cancer scRNA-seq atlas | Reference for every reference-based competitor method; source cells for the scCube simulation | NCBI GEO **GSE176078** |
| Xenium Human Kidney Dataset (renal cell carcinoma) | Figure 4 real-tissue benchmark | 10x Genomics, publicly available |

Two inputs are **not** covered by the three sources above, and the manuscript does not record
an accession for them:

- The **breast cancer Xenium** dataset used for Figure 2. The manuscript describes it as
  breast cancer Xenium data aggregated to pseudo-Visium spots; it is not listed in the
  paper's Deposited Data table.
- **GSE89888** (MCF7/T47D ESR1-D538G RNA-seq), read directly by the Figure 6B panel function
  from a locally downloaded TPM table.

### Paths are placeholders

Every script carries hardcoded input/output roots as module-level constants, following the
same convention as [`../README.md`](../README.md). They ship as `/path/to/...` placeholders
and **must be edited before the script will run**:

- `/path/to/CITEgeist_analysis` — root of your local copy of the analysis outputs. The
  scripts expect the layout mirrored by this directory: `benchmarks/…`, `mdk_analysis/…`,
  `figures/…`, plus an `output/` tree holding per-module CITEgeist results.
- `/path/to/CITEgeist_public_data/processed_files` — the GSE289326 SpaceRanger output
  directories (`<specimen>/outs/`).
- `/path/to/Wu_2021_BRCA_scRNASeq` — the GSE176078 download.
- `/path/to/Xenium_BreastCancer` — the breast Xenium download.
- `/path/to/imaging_data` — raw immunofluorescence TIFFs for the Figure 6C montage.

Specimen identifiers (`HCC22-088-P1-S1`, `HCC22-088-P4-S2_1i_rep`, …) are the anonymized
labels used in the manuscript and are left intact so the scripts read the way the paper does.

### Environment

`pip install citegeist` gives you the CITEgeist package, but these analysis scripts also need
the competitor methods and a few analysis libraries that are **not** CITEgeist dependencies:

| Needed by | Extra packages |
|---|---|
| `benchmarks/simulation_benchmarking/` wrappers | `scikit-learn`; RCTD/CARD/Seurat runs are driven from R (`spacexr`, `CARD`, `Seurat`) — the Python wrappers score their outputs |
| `benchmarks/competitor_gex_validation/run_c2l_*.py` | `cell2location` (GPU strongly recommended) |
| `benchmarks/competitor_gex_validation/run_tangram_mapping.py` | `tangram-sc` |
| `benchmarks/competitor_gex_validation/run_commot.py`, `mdk_analysis/scripts/*commot*.py` | `commot`, `statsmodels` |
| `benchmarks/competitor_gex_validation/run_bivariate_morans.py` | `esda`, `libpysal` |
| `mdk_analysis/scripts/analyze_p4s2_commot.py` | `gseapy`, `statsmodels` |
| `figures/` | `matplotlib`, `adjustText`, `lxml`, `squidpy`, `tifffile`; `rsvg-convert` (or `cairosvg`) for PNG/PDF rasterization |

**GPU:** nothing in this directory runs the CITEgeist QP solver, so no GPU and no cuOPT is
required to re-run the *evaluations*. Regenerating the CITEgeist predictions that these
evaluations consume does require a GPU with cuOPT — see [`../README.md`](../README.md).

The manuscript figures were typeset in Arial. Set `CITEGEIST_FIGURE_FONT=/path/to/Arial.ttf`
to reproduce them exactly; without it matplotlib's default sans-serif is used and text metrics
will differ slightly.

---

## Script → figure map

### Figure 2 — breast cancer Xenium benchmark (competitors only; CITEgeist is absent because that dataset has no antibody capture)

| Script | Produces | Feeds |
|---|---|---|
| `benchmarks/xenium_benchmarking/evaluation/src/evaluate_benchmark.py` | `results_breast/<method>_results.json` (run once per method with `--tissue breast --method <M>`) | Fig 2A, 2B |
| `benchmarks/xenium_benchmarking/evaluation/src/compare_breast_gex_methods.py` | `results_breast/gex_comparison.json` | Fig 2C, 2D |
| `figures/breast_benchmarking/generate.py` | panels A–D | Figure 2 |

Supporting: `benchmarks/xenium_benchmarking/breast_constants.py` (tissue configs and the
Xenium-cluster → atlas cell-type mapping), imported by `evaluate_benchmark.py`.

### Figure 3 — simulated (scCube) benchmark

| Script | Produces | Feeds |
|---|---|---|
| `benchmarks/simulation_benchmarking/src/compare_prop_methods.py` | proportion metrics for all six methods | Fig 3A |
| `benchmarks/simulation_benchmarking/src/compare_gex_methods.py` | per-cell-type GEX metrics | Fig 3B |
| `benchmarks/simulation_benchmarking/src/{citegeist,cell2location,RCTD,CARD,seurat,tangram}_bench_wrapper.py` | per-method scoring of each method's prediction directory | inputs to the two comparison scripts above |
| `figures/simulated_benchmarking/generate.py` | panels A–E | Figure 3 |

Supporting: `benchmarking_spot_deconv.py` (proportion metrics: Pearson *r*, RMSE, MAE,
Jensen–Shannon) and `benchmarking_gex.py` (GEX RMSE/NRMSE/MAE) — every wrapper imports one
or both.

### Building the pseudo-Visium benchmark (run before Figures 2 and 4)

Both Xenium benchmarks score methods on *pseudo-Visium* spots aggregated from Xenium
single-cell data, against ground truth derived from the same cells. `benchmarks/xenium_pseudovisium/src/`
builds both. Run this before the evaluators below — they read its output and will fail
without it.

| Step | Script | Produces |
|---|---|---|
| 1 | `load_xenium.py` | Xenium cells → AnnData |
| 2 | `define_cell_types.py`, `rna_cell_types.py` | per-cell type assignments (protein-gated and RNA-based) |
| 3 | `create_pseudo_spots.py` | hexagonal grid matching Visium geometry (**55 µm spot diameter, 100 µm center-to-center**) and cell→spot count aggregation |
| 4 | `split_regions.py` | 5 non-overlapping regions used as pseudo-replicates |
| 5 | `create_protein_gt.py` | `data_protein_gt/` — `cell_to_spot_mapping.csv`, `cell_type_assignments.csv`, per-region proportion CSVs |
| 6 | `generate_protein_gt_gex.py` | `ground_truth_gex/` — per-cell-type expression matrices |

`generate_ground_truth.py`, `generate_expression_weighted_gt.py` and
`generate_gex_ground_truth.py` are the shared proportion/expression GT builders these steps call.

For the **breast** arm (Figure 2), `benchmarks/xenium_benchmarking/prepare_breast_data.py` runs
the equivalent pipeline end-to-end (extract → annotate → spots → 5 regions → GT), and
`convert_breast_h5ad.py` exports per-region counts/coords CSVs for the R-based methods.

The 6-type and 7-type SingleR ground truth used for the *reported* Figure 4 numbers is built by
`evaluation/src/generate_singler_gex_gt.py`, listed below.

### Figure 4 — Xenium renal cell carcinoma real-tissue benchmark

| Script | Produces | Feeds |
|---|---|---|
| `benchmarks/xenium_benchmarking/evaluation/src/generate_singler_gex_gt.py` | SingleR-based per-cell-type GEX ground truth | ground truth for the GEX comparison |
| `benchmarks/xenium_benchmarking/evaluation/src/compare_all_methods.py` | `method_comparison_singler_6type/` and `…_7type/` result JSONs | Fig 4A, 4B |
| `benchmarks/xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py` | `xenium_gex_comparison_7type_summary.csv` | Fig 4C, 4D, 4E |
| `figures/xenium_benchmarking/generate.py` | panels A–E | Figure 4 |

Supporting: `evaluate_benchmark.py` (the per-region metric engine `compare_all_methods.py`
calls) and `benchmark_constants.py` (the 6-type and 7-type ground-truth collapse mappings).

**Reported values to compare against** (5 Xenium regions, CellTypist→SingleR annotated
ground truth; Pearson *r* on spot-level proportions):

| Method | 6-type *r* | 7-type *r* |
|---|---|---|
| **CITEgeist** | **0.763** | **0.722** |
| Cell2Location | 0.639 | 0.624 |
| RCTD | 0.517 | 0.589 |
| Seurat | 0.392 | 0.442 |
| CARD | 0.385 | −0.190 |
| Tangram | 0.289 | 0.262 |

GEX, expressed genes: flattened *r* — CITEgeist 0.676, Cell2Location 0.625. RMSE
(per-block ground-truth-scale normalized) — CITEgeist 0.420, Cell2Location 0.554,
Tangram 0.674. On Cell2Location's own 180 matched genes: CITEgeist 0.727 vs 0.626.

Use the `*_6type_results.json` files for the 6-type ground truth and `*_7type_results.json`
for the 7-type one. Mixing them (7-type predictions against 6-type ground truth) splits the
T-cell compartment against itself and understates every method.

### Figure 5 — MDK spatial co-localization

| Script | Produces | Feeds |
|---|---|---|
| `figures/mdk_analysis/generate_panel_A_rasters.py` | spatial raster grid (cancer proportion + early/late E2 response, 3×4) | Fig 5A |
| `mdk_analysis/scripts/commot_hurdle_test.py` | `hurdle_results_summary.csv`, `hurdle_results.json` | Fig 5B |
| `mdk_analysis/scripts/analyze_p4s2_commot.py` | per-cell COMMOT ligand–receptor analysis + enrichment for specimen P4-S2 | Fig 5B, Supplementary MDK/COMMOT panels |
| `benchmarks/competitor_gex_validation/` (whole pipeline, below) | `morans_comparison.csv`, `failure_mode_summary.csv` | Fig 5C, 5D |
| `figures/mdk_analysis/generate.py --figure spatial` | panels A–E | Figure 5 |

Directional oestrogen-response scoring for panel A lives in `figures/mdk_analysis/_e2_scoring.py`
(`score_genes(UP) − score_genes(DN)` over the four MSigDB Li 2023 EstroGene sets).

The COMMOT data are heavily zero-inflated, so Figure 5B uses a hurdle model — Fisher exact on
sender participation plus Mann–Whitney on the non-zero senders. MDK–NCL is the strongest pair
(2.1×, *p* = 9.2 × 10⁻¹⁰).

#### The competitor-GEX validation pipeline (Figure 5C, 5D)

Run in order; each step is a separate script because each was a separate array job:

1. `prepare_erpos_reference.py` — filter the Wu et al. atlas to ER+ patients, remap cell types
2. `run_c2l_reference_training.py` — Cell2Location regression model on that reference
3. `run_c2l_mapping.py` — Cell2Location spot-level mapping, per specimen
4. `run_tangram_mapping.py` — Tangram cell-to-space mapping, per specimen
5. `expand_pseudocells.py` — expand each method's spot-level output into pseudo-single-cells
6. `run_bivariate_morans.py` — bivariate Moran's I of MDK against the secretory composite, using the *same* metric CITEgeist is scored on
7. `run_commot.py` — COMMOT ligand–receptor scores on the pseudo-single-cell data
8. `run_module4_programs.py` — standalone NMF program discovery per cell-type layer (deliberately does not import the CITEgeist package, so it needs no cuOPT)
9. `aggregate_comparison.py` — collects 6–8, classifies failure modes, writes the comparison tables

`constants.py` holds the shared specimen list, cell types, candidate genes and secretory gene
set for all nine.

### Figure 6 — ESR1-D538G regulation of MDK secretion

| Script | Produces | Feeds |
|---|---|---|
| `figures/mdk_analysis/generate.py --figure mechanism` | ER ChIP-seq heatmap (A), secretory RNA log2 fold-change (B), ELISA (D) | Figure 6 |
| `figures/mdk_analysis/generate_if_panel.py --label C` | immunofluorescence micrograph montage | Fig 6C |

Panel B is computed inside `generate.py` directly from the GSE89888 TPM table; panels A and D
read tabulated ChIP-seq and ELISA measurements.

### Bivariate Moran's I producers (Figures 5C, 5E, and supplementary)

These compute the spatial-correlation tables the figure scripts above consume. Run them before
the figures that read their output.

| Script | Produces | Feeds |
|---|---|---|
| `mdk_analysis/scripts/persist_bivariate_morans.py` | `if_analysis/bivariate_morans_results.csv` — cancer-only bivariate Moran's I of MDK vs the secretory score, per specimen (mean 0.254) | Fig 5C (the CITEgeist baseline competitors are compared against) |
| `mdk_analysis/scripts/spatial_morans.py` | `spatial/trajectory_deltas.csv` — biopsy→surgery change in Moran's I per patient | Fig 5E |
| `mdk_analysis/scripts/compute_pergene_bivariate_morans.py` | `if_analysis/pergene_bivariate_morans.csv` — MDK vs each of 11 secretory genes across 12 specimens, FDR-corrected, with 11 random non-secretory controls | the supplementary heatmap below |
| `mdk_analysis/scripts/generate_pergene_morans_heatmap.py` | per-gene bivariate Moran's I heatmap (MDK vs each secretory gene, per specimen) | Supplementary |

`spatial_morans.py` reads the sample inclusion and pairing manifests shipped at
`mdk_analysis/data/{sample_manifest.csv,sample_pairs.csv}`; override with `--sample-manifest` /
`--sample-pairs`.

> **One naming caveat, worth stating plainly.** `spatial_morans.py` also writes
> `cancer_cell_morans.csv`, and despite the name that composite is computed over **all cell
> types** (mean Moran's I 0.290) — not cancer cells only. The cancer-only value reported in the
> manuscript (mean 0.254) is the one from `persist_bivariate_morans.py`. The same note is in the
> script's own docstring. Use `bivariate_morans_results.csv` for the cancer-only figure.

---

## What is deliberately not here

Being explicit, because an omission that looks accidental is worse than one that is stated:

- **The figure assembly engine.** These scripts generate individual *sub-panels* as SVG/PNG.
  Stitching panels into the final composite figures, journal font enforcement and PDF export
  are publication plumbing, not analysis, and are not included. Panel sizes ship as
  `panel_sizes.json` next to each generator so panels render at their published dimensions.
- **Superseded and exploratory scripts.** The working tree contains many more analysis
  attempts than the paper reports. Only the versions that produced reported numbers are here.
  Where two scripts could plausibly have produced a figure, the one that is *not* read by the
  shipped figure generator was dropped rather than shipped alongside it.
- **SLURM launchers.** Every script is a plain command-line program; the `sbatch` wrappers
  were specific to the authors' cluster and would not run anywhere else.
