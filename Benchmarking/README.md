# Benchmarking

This directory holds the code behind the benchmark figures and tables of the CITEgeist paper
(CITEgeist 3.0.0). It compares CITEgeist with RCTD, Cell2Location, Seurat, CARD and Tangram on a
simulated benchmark and on real Xenium tissue, and holds the protein-sensitivity and midkine
analyses.

It is not part of the `pip install citegeist` package. Use a clone of the repository.

The scripts are SLURM jobs written for one cluster. Adjust the partition, account and module lines
in the `.sh` and `.sbatch` files for your site.

## Figure and table map

| Item | Run | Score |
|---|---|---|
| CITEgeist proportions and expression (Xenium) | `xenium_benchmarking/spots_baseline/write_marker_share_preds_p138.sbatch`, `sace_props_cpu.sbatch` | `xenium_benchmarking/gt_rebuild/score_methods.py`, `xenium_benchmarking/evaluation/src/evaluate_benchmark.py`, `xenium_benchmarking/spots_baseline/score_p138_strict.py`, `xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py` |
| RCTD, Cell2Location, Seurat, CARD, Tangram (Xenium) | `xenium_benchmarking/<Method>/slurm/xenium_benchmark_singler_{6,7}type.sh` | same scorers as above |
| Xenium scoring and ground truth | `xenium_benchmarking/gt_rebuild/rnaprot.py`, `record_inputs.py` | `xenium_benchmarking/gt_rebuild/score_methods.py`, `xenium_benchmarking/evaluation/src/{evaluate_benchmark,compare_xenium_gex_methods,build_native_and_split_layers,matched_gene_r,build_canonical_results}.py` |
| Simulated benchmark | `simulation_benchmarking/regen/slurm/`, `regen/src/`, `regen/diag/d1/` | `simulation_benchmarking/score_calsim_per_type.py` |
| Protein-sensitivity analysis | `protein_sensitivity/run_sensitivity.py`, `slurm/sensitivity_v41p.sbatch` | `protein_sensitivity/aggregate_v41p.py` |
| Midkine analysis | `midkine_remediation/src/` | `midkine_remediation/src/a12_9_real.py`, `a12_10_real.py` |
| GeoMx (legacy, see below) | `geomx_benchmarking/CITEgeist/src/run_pseudo_geomx.py` | `geomx_benchmarking/evaluation/slurm/evaluate_all.sh` |

## Environment variables

Paths come from `_paths.sh`, `_paths.py` and `_paths.R`. Each root is read from an environment
variable and falls back to a default.

| Variable | Meaning | Default |
|---|---|---|
| `CITEGEIST_REPO_ROOT` | Root of this repository | the clone containing `Benchmarking/` |
| `CITEGEIST_DATA_ROOT` | Input data (layout below) | `<repo>/data` |
| `CITEGEIST_REF_ROOT` | Reference data | `<repo>/ref` |
| `CITEGEIST_OUT_ROOT` | Outputs written by the runners | `<repo>/out` |
| `CITEGEIST_ENVS_ROOT` | Directory holding the conda environments | `$(conda info --base)/envs` |

`CITEGEIST_ENVS_ROOT` needs conda on `PATH` (or the variable set). If neither is available, scripts
that need an environment stop with a message telling you to set it. Make sure `conda` is on `PATH`
inside your sbatch jobs.

`_paths.sh` also reads `~/.citegeist_local.env` if it exists, so you can keep your overrides there.

## Data layout

Put the downloaded data under `CITEGEIST_DATA_ROOT` like this:

```
$CITEGEIST_DATA_ROOT/
├── GSE289326/processed_files/<specimen>/outs/   # Visium + CITE-seq (SpaceRanger output)
├── 10x_PublicData/
│   ├── Xenium_RNA_Proteomic_RenalCellCarcinoma/ # Xenium RCC download
│   └── Xenium_RNA_BreastCancer/                 # breast Xenium download
├── BC-Datasets/Wu_2021/Wu_etal_2021_BRCA_scRNASeq/   # GSE176078
├── Wu_Visium/Simulations/scCube_12k/            # simulation replicates and per-method outputs
├── simulated_benchmark_cache/                   # cached expected means for simulation checks
└── midkine_outputs/                             # midkine analysis tables (e.g. if_analysis/)
```

The GSE156632 RCC reference is processed by the scripts in
`xenium_benchmarking/reference_data/GSE156632/`.

## Data sources

| Dataset | Used for | Source |
|---|---|---|
| Visium + CITE-seq clinical-trial samples (this study) | Midkine analyses; competitor expression validation | NCBI GEO [GSE289326](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289326) |
| RCC scRNA-seq reference | Reference for the reference-based methods on Xenium | NCBI GEO [GSE156632](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156632) |
| Wu et al. 2021 breast cancer scRNA-seq atlas | Reference for the simulation and breast benchmark | NCBI GEO [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) |
| Xenium RNA + Proteomic Renal Cell Carcinoma | Real-tissue benchmark | public 10x Genomics dataset, [10x Genomics](https://www.10xgenomics.com/datasets) |

## Known limits for public reproducers

- The public repository holds code only. The Xenium ground truth (GT) comes from the paper's data
  bundle, together with the signed GT sign-off files. Use the bundle's GT. Do not rebuild it.
- The GT-rebuild stages (`gt_rebuild/rnaprot.py`, `record_inputs.py` and the other builders) need the
  lab's sign-off. You cannot re-sign it outside the lab, so these stages will refuse to run.
- The scorers (`gt_rebuild/score_methods.py`, the Xenium evaluation scripts,
  `spots_baseline/score_p138_strict.py` and the `protein_sensitivity` scorers) call
  `gt_rebuild/manifest.py::require_rnagt_signoff` first. They run as-is when the bundle's files are in
  place and unchanged.
- Where the guard looks: paths are relative to the repository checkout, not to a `CITEGEIST_*_ROOT`
  variable. Unpack the bundle's GT into `Benchmarking/xenium_benchmarking/`:
  - `ground_truth_xenium_clusters/manifest.json`, `reports/integrity_report.json` and
    `reports/blinded_review.md` (under that same folder), plus every other report file the
    manifest lists under `signoff.reports_seen_sha256`
  - `ground_truth_xenium_rnaprot/manifest.json` and `_COMPLETE.json` (for the `rnaprot` label sets)
  - the GT tables the scorers read, `ground_truth_xenium_{clusters,rnaprot}/ground_truth_{6,7}type/ground_truth/Xenium_region_<r>_prop.csv`
    (and `ground_truth_gex/` for expression scoring)
  - Every file must be byte-identical to the signed one. The guard checks their SHA-256 hashes.
- The frozen midkine protocol JSON files keep some internal document paths in text fields so their
  hashes stay stable.

## Environments

Per-method conda environment files: `xenium_benchmarking/reference_data/envs/` (RCTD and Seurat share
`r_deconv_env.yml`; Cell2Location; Tangram) and `environments/card_env.yml` (CARD). CITEgeist, the scoring
scripts, the protein-sensitivity analysis and the midkine analysis use the `citegeist` package environment
(see the repository `pyproject.toml`).

## Running CITEgeist

Xenium: `xenium_benchmarking/spots_baseline/write_marker_share_preds_p138.sbatch` writes the marker-share
proportions, then `sace_props_cpu.sbatch` runs the SACE GEX allocation. GeoMx: `geomx_benchmarking/CITEgeist/src/run_pseudo_geomx.py`
and `run_geomx_hnscc.py`.

## Running RCTD

`xenium_benchmarking/RCTD/slurm/xenium_benchmark_singler_{6,7}type.sh` (after the reference preparation in
`xenium_benchmarking/reference_data/GSE156632/`).

## Running Cell2Location

`xenium_benchmarking/Cell2Location/slurm/train_reference*.sh`, then `xenium_benchmark_singler_{6,7}type.sh`.

## Running Seurat

`xenium_benchmarking/Seurat/slurm/xenium_benchmark_singler_{6,7}type.sh`.

## Running CARD

`xenium_benchmarking/CARD/slurm/xenium_benchmark_singler_{6,7}type.sh`.

## Running Tangram

`xenium_benchmarking/Tangram/slurm/xenium_benchmark_singler_{6,7}type.sh`.

## Scoring the Xenium benchmark

`xenium_benchmarking/gt_rebuild/score_methods.py` (proportions, joint RNA+protein ground truth) and
`xenium_benchmarking/evaluation/src/compare_xenium_gex_methods.py` (expression). See "Known limits"
above for the ground-truth stages you cannot re-run.

## Scoring the GeoMx benchmark

> **Legacy.** Not reported in the v41 manuscript; uses the earlier QP (cuOPT) proportion solver (method pinned to qp); provided for completeness.

`geomx_benchmarking/evaluation/slurm/evaluate_all.sh`; the real GeoMx HNSCC evaluation is
`geomx_benchmarking/CITEgeist/src/evaluate_geomx_hnscc.py`.

## Scoring the simulated benchmark

`simulation_benchmarking/score_calsim_per_type.py`. The competitor predictions it reads were produced by
the runners in `simulation_benchmarking/regen/slurm/` and `regen/src/` (Tangram, Cell2Location) and
`regen/diag/d1/` (RCTD, Seurat, CARD with the full-gene reference), then gathered by
`regen/src/promote_canonical.py`.

## Protein-sensitivity analysis

`protein_sensitivity/run_sensitivity.py`, scored by `protein_sensitivity/aggregate_v41p.py`.

## Midkine analysis

Stage scripts under `midkine_remediation/src/` (marker-share proportions, single-cell rebuild, SACE layers,
compartments, COMMOT, NMF), with the Moran's I and hurdle statistics in `a12_10_real.py` and `a12_9_real.py`.
