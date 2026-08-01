Changelog
=========

All notable changes to CITEgeist will be documented in this file.

Version 2.0.1 (2026-08-01)
---------------------------

A documentation and packaging-correctness release. No changes to the deconvolution
algorithms or to any numerical result.

Fixed
~~~~~

* ``CITEgeist.model.qc`` raised a ``NotImplementedError`` claiming the report module
  was "not yet available" and referencing an internal task number. The module ships
  and works; the message was stale and is baked into the 2.0.0 artifact.
* Module 3.5 (functional annotation) imported ``torch``, which is not a CITEgeist
  dependency, producing a bare ``ImportError``. It is now gated behind a new optional
  extra with an actionable message: ``pip install citegeist[functional]``. Module 3.5
  is experimental and is not part of the production QP + SACE pipeline.
* Removed ``register_gurobi`` and the ``gurobipy`` pin. Gurobi is an archived dead end;
  the production solver is cuOPT.
* Corrected five documentation sites that stated ``CITEgeist_env.yml`` includes cuOPT.
  It does not, and neither does it include StarDist — both must be installed separately.
  There is no CPU fallback for the proportion solver.
* Rewrote the Sphinx quickstart, tutorial and examples pages, which documented an API
  that does not exist (a ``run_cell_expression_pass1(radius=...)`` signature, the wrong
  cell-profile format, and a tissue type CITEgeist was never validated on). Every
  signature is now verified against the installed code.
* Corrected the documented default for single-cell assignment: ``use_morphology=True``
  routes through the morphology-informed Bayesian posterior; Hungarian matching is
  available via ``use_morphology=False``.

Added
~~~~~

* Data availability: the Visium clinical-trial data are public at NCBI GEO accession
  ``GSE289326``; the reference atlas is ``GSE176078``.
* ``examples/paper_analysis/`` — the curated analysis code behind the paper's reported
  numbers, with a map from each script to the figure it produces.
* A Read the Docs build that fits within the hosted build limits.

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
* **Per-cell Type Assignment**: StarDist nuclei segmentation, then morphology-informed Bayesian assignment from spot proportions by default (``use_morphology=True``); count-constrained Hungarian matching via ``use_morphology=False``
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
