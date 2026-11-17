Changelog
=========

All notable changes to CITEgeist will be documented in this file.

Version 3.0.0
-------------

Breaking changes
~~~~~~~~~~~~~~~~

* ``run_cell_proportion_model()`` now defaults to ``method="auto"``, which resolves to
  ``"rna_refined"``: SPOTS marker share (Ben-Chetrit et al. 2023) followed by
  reference-free RNA refinement (multinomial, c = 0.1, 3 passes, no anchor). CPU only, no
  GPU needed. The output is an RNA compartment (an RNA-mass share), not a cell count.
* ``auto`` needs a GEX AnnData with a ``raw_counts`` layer and per-spot nuclei counts
  (``nuclei_counts=`` as a barcode-indexed Series, or ``obs["n_nuclei"]``). Missing inputs
  raise ``ValueError``; there is no silent fallback. It also raises with ``use_gating=True``.
* To get 2.0.0 behaviour, pass ``method="qp"`` (cuOPT, GPU only; same API and results).
  ``method="marker_share"`` (CPU, no nuclei needed) is also available explicitly.
* ``run_sace_allocation(use_morphology=...)`` now defaults to ``False``
  (count-constrained Hungarian). ``use_morphology=True`` is unchanged.
* ``filter_gex(keep_unfiltered=False)``: keeping the unfiltered GEX is opt-in.

Added / fixed
~~~~~~~~~~~~~

* ``brightness_method`` argument (``None``/"plain" default, ``"gmm"``, RNA-gated methods);
  ``brightness_correction=True`` is a deprecated alias for ``"gmm"``.
* Marker share reads raw ADT counts (``layers["raw_counts"]``) when present.
* New dependency ``scikit-misc``. New experimental extra ``citegeist[functional]``
  (``torch``). ``citegeist[imaging]`` is StarDist.
* Experimental dictionary discovery and guided state discovery APIs.

See ``CHANGELOG.md`` in the repository for details.

Version 2.0.2 (2026-11-01)
---------------------------

A bug-fix release. Results are unchanged for any dataset whose spots were already in
sorted barcode order, which includes standard SpaceRanger output.

Fixed
~~~~~

* ``run_cell_proportion_model`` stores proportions sorted by barcode, but the SACE
  steps read them by position. When the spots in the GEX AnnData were not in sorted
  barcode order and the row counts matched, each spot silently received another spot's
  proportions. ``run_sace_allocation`` (layers and single-cell modes), ``assign_cells``
  and ``run_sace_protein`` now align proportions to the GEX spots by barcode.
* ``assign_cells`` reused the stored ``cell_counts_integer`` table in antibody spot order
  while reading ``cell_to_spot`` as a position in the (sorted) proportions, so cells could
  be given another spot's types. It is now aligned by barcode.
* ``run_module3_5_functional_annotation`` paired antibody counts with sorted proportions and
  GEX coordinates by position. It now aligns them by barcode.
* The GEX-informed detection step of ``run_cell_proportion_model`` paired GEX rows with
  antibody rows by position. It now pairs them by barcode, so proportions are correct even
  when the antibody and GEX objects list their spots in different orders.
* ``assign_cells`` accepts an optional ``cell_prop`` argument (proportions in the row
  order that ``cell_to_spot`` indexes). For the Bayesian method, the stored detection
  mask and proportion prior are aligned to those rows by barcode.

Changed
~~~~~~~

* A mismatch between proportion barcodes and GEX barcodes (missing, extra or duplicated
  barcodes) now raises ``ValueError`` instead of being paired by position. In
  ``run_sace_allocation`` and ``run_sace_protein`` the proportion columns must also be
  the same set of cell types as the cell profile dictionary.
* ``method="per_type_beta"`` can name its cell types differently from your profile dictionary
  (it uses its own fixed type names when the panel has the extra markers). Calling
  ``run_sace_allocation`` or ``run_sace_protein`` after such a run now raises ``ValueError``
  ("type columns differ; missing [...], extra [...]") instead of continuing silently.
  Rename the columns of ``model.results["cell_prop"]`` to your profile names first if the
  types correspond.

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
