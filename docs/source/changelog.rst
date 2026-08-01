Changelog
=========

All notable changes to CITEgeist will be documented in this file.

Version 2.0.0 (2026-07-31)
---------------------------

CITEgeist is published on PyPI (see :doc:`installation`). This section describes the state
of the 2.0 codebase.

* Production stack: **cuOPT QP** for cell-type proportion estimation (Module 3) +
  **SACE GEX** for per-cell gene-expression deconvolution (Module 3-GEX)
* Requires Python 3.10+ and an NVIDIA GPU with CUDA for the QP solver (no CPU fallback)
* Runtime dependencies are now fully declared in ``pyproject.toml``
* The installed distribution is named ``citegeist`` (lowercase); the importable package is ``CITEgeist``

Features
~~~~~~~~

* **CitegeistModel**: Main class orchestrating all modules (Modules 1–5)
* **Cell Proportion Optimization**: GPU-accelerated QP deconvolution via cuOPT
* **SACE GEX Deconvolution**: Single-pass Poisson-multinomial per-cell GEX allocation
* **Per-cell Type Assignment**: StarDist nuclei segmentation + Hungarian matching from spot proportions (optional Bayesian assignment when morphology scores are precomputed)
* **Spatial Gene Program Discovery**: NMF-based programs with bivariate Moran's I
* **Cross-sample Integration**: Harmony-style program alignment across specimens

API Reference
~~~~~~~~~~~~~

* ``CITEgeist.CitegeistModel``: Top-level model class (lazy re-export)
* ``CITEgeist.model.citegeist_model.CitegeistModel``: Full import path
* ``CITEgeist.model.deconvolution.qp_solver``: cuOPT QP solver (GPU required)
* ``CITEgeist.model.gex.sace_gex``: SACE GEX deconvolution
* ``CITEgeist.model.assignment.cell_assignment``: Per-cell type assignment

Dependencies
~~~~~~~~~~~~

* Python 3.10+
* numpy, pandas, scipy, scanpy, squidpy, anndata
* scikit-learn, scikit-image, matplotlib, networkx, tqdm, joblib
* esda, libpysal
* cuOPT (NVIDIA RAPIDS/NGC — not on PyPI, GPU required)
