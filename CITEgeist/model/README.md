# CITEgeist Model Package

Spatial multi-omics deconvolution framework: CITE-seq antibody profiles + Visium GEX → cell-type proportions, per-cell expression, and spatial programs.

## Package Structure

```
model/
├── citegeist_model.py          # Main CitegeistModel orchestrator
├── utils.py                    # Shared utilities (logging, I/O, memory)
├── checkpoints.py              # Checkpoint manager
├── unified_config.py           # Pipeline configuration constants
├── module2_proposal_builder.py # M2 profile proposal construction
├── proposal_review_loader.py   # Proposal review result loader
│
├── discovery/                  # M1 + M2: Marker scoring & profile building
│   ├── marker_interest.py      #   Kurtosis, GMM SNR, Moran's I
│   └── spatial_colocalization.py #  Colocalization, profile discovery/selection
│
├── deconvolution/              # M3: QP solver & cell-type detection
│   ├── qp_solver.py            #   QP-based cell proportion solver (production)
│   ├── detection.py            #   GMM-based detection gating
│   ├── detection_refinement.py #   GEX-informed detection fusion
│   ├── emission_init.py        #   Marker config, beta init
│   └── reference_model.py      #   Optional NB reference refinement
│
├── assignment/                 # M3-post: Spot → cell resolution
│   ├── cell_assignment.py      #   Bayesian/Hungarian assignment (production)
│   ├── hungarian_assignment.py #   Cost-matrix solver
│   ├── module3b_nucleus_assignment.py
│   ├── cellularity_utils.py    #   Cellularity scaling helpers
│   └── single_cell_output.py   #   AnnData assembly
│
├── gex/                        # M3-gex: Per-cell expression allocation
│   ├── sace_gex.py             #   SACE Poisson-multinomial EM (production)
│   ├── cell_level_gex.py       #   Simple proportional distribution
│   └── gex_modules.py          #   Anchor genes, module-aware enrichment
│
├── annotation/                 # M3.5: Protein gating & subtype splitting
│   ├── functional_annotation.py #  3-stage GMM gating
│   ├── subtype_splitting.py    #   split_by_protein_gates()
│   ├── module3_5_benchmark.py  #   Aggregation helpers
│   ├── module3_5_projection.py #   Enrichment projection checks
│   └── coverage_check.py       #   Module coverage validation
│
├── programs/                   # M4 + M5: Programs & integration
│   ├── anchored_program_discovery.py  # NMF programs, bivariate analysis
│   └── cross_sample_integration.py    # Harmony alignment, conservation
│
├── qc/                         # Quality control checks
│   ├── canonical_markers.py    #   Canonical marker validation
│   ├── gex_qc.py              #   GEX quality metrics
│   ├── marker_enrichment.py   #   Marker enrichment analysis
│   ├── proportion_qc.py       #   Proportion QC checks
│   ├── report.py              #   QC report generation
│   └── single_cell_qc.py     #   Single-cell QC
│
├── morphology/                 # Vision, segmentation, SSL
│   ├── segmentation.py         #   StarDist patchwise (production)
│   ├── morphology_backbone.py  #   DAPI/HE backbone wrappers
│   ├── morphology_features.py  #   Nucleus/cell feature extraction
│   ├── morphology_prior.py     #   LLP morphology prior
│   ├── patch_extraction.py     #   Image patch extraction
│   ├── prototype_contrastive.py #  Prototype-contrastive LLP
│   ├── simclr.py               #   SimCLR training backbone
│   ├── soft_label_classifier.py
│   ├── vit_encoder.py          #   ViT-Small encoder
│   └── vit_extractor.py        #   ViT feature extraction (UNI)
│
└── _archive/                   # Dead-end experiments (git-excluded)
```

## Pipeline Flow

```
Discovery (M1→M2)
    │  marker scoring → colocalization → profile building
    ▼
Deconvolution (M3)
    │  QP proportions (cuOPT GPU) + detection gating
    ▼
Assignment (M3-post)          GEX Allocation (M3-gex)
    │  spot→cell typing           │  SACE per-cell expression
    ▼                             ▼
Annotation (M3.5)
    │  protein gating → subtype splitting
    ▼
Programs (M4→M5)
    │  NMF gene programs → cross-sample integration
```

## Entry Point

All public API is accessible via `from CITEgeist.model import ...` (lazy-loaded).

For direct imports: `from CITEgeist.model.<subpackage>.<module> import ...`

Example:
```python
from CITEgeist.model import CitegeistModel  # via __init__.py
from CITEgeist.model.deconvolution.qp_solver import optimize_cell_proportions_per_marker  # direct
```
