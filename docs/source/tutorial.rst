Tutorial
========

This tutorial walks through the CITEgeist pipeline on real Visium CITE-seq data.
It mirrors
`docs/quickstart_real_visium.md <https://github.com/leeoesterreich/CITEgeist/blob/main/docs/quickstart_real_visium.md>`_,
which carries the same workflow in full with troubleshooting notes.

.. warning::

   Modules 3 and 3-post require an **NVIDIA GPU**. The proportion solver uses
   cuOPT, which is not on PyPI and is not in ``CITEgeist_env.yml``; nucleus
   segmentation uses StarDist (``pip install citegeist[imaging]``). There is no CPU
   fallback for the proportion solver. See :doc:`installation`.

Data Preparation
----------------

CITEgeist reads a SpaceRanger output directory containing both gene-expression and
antibody-capture features:

* ``filtered_feature_bc_matrix.h5`` loaded with ``gex_only=False``
* a ``spatial/`` folder; for single-cell output, an H&E image at
  ≤ 1.0 µm/px (``tissue_fullres_image.tif``, with CytAssist/hires fallbacks)
* spatial coordinates in ``adata.obsm["spatial"]`` — note that
  ``sq.read.visium(load_images=False)`` does **not** populate this, so either load
  images or inject the coordinates from ``spatial/tissue_positions.csv``
* feature types in ``adata.var["feature_types"]``

Antibody feature names in real SpaceRanger output carry a ``-1`` suffix (for example
``CD68-1``). Simulated data does not.

No external scRNA-seq reference dataset is needed.

Workflow Overview
-----------------

1. **Profile discovery** (Modules 1–2, optional) — find spatially informative
   antibody markers and group them into candidate profiles.
2. **Preprocessing** — split the AnnData and normalize both modalities.
3. **Spot-level proportions** (Module 3) — cuOPT QP against the antibody profiles.
4. **Single-cell assignment and GEX** (Modules 3-post / 3-gex) — StarDist nuclei,
   per-nucleus type assignment, and SACE allocation of spot counts to cells.

Step 1: Profile discovery (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modules 1–2 are a functional API, not methods on the model class. Run once per
cohort or tissue type; the resulting profiles are reusable across samples.

.. code-block:: python

   import numpy as np
   import scipy.sparse
   import squidpy as sq

   from CITEgeist.model import (
       identify_interesting_markers,
       analyze_marker_colocalization,
       select_profiles,
   )
   from CITEgeist.model.discovery.spatial_colocalization import discover_profiles_continuous

   adata = sq.read.visium(
       "/path/to/spaceranger/output/sample_P1_S1",
       counts_file="filtered_feature_bc_matrix.h5",
       load_images=True,
       gex_only=False,
   )

   ab_adata = adata[:, adata.var_names.str.endswith("-1")]
   X = ab_adata.X.toarray() if scipy.sparse.issparse(ab_adata.X) else np.asarray(ab_adata.X)
   coords = np.asarray(adata.obsm["spatial"], dtype=np.float64)
   marker_names = list(ab_adata.var_names)

   m1 = identify_interesting_markers(
       X=X, coords=coords, marker_names=marker_names, morans_k=8, morans_n_perm=99
   )
   coloc = analyze_marker_colocalization(
       X=X, coords=coords, marker_names=marker_names,
       markers_to_analyze=m1.interesting_markers,
       neighbor_k=6, n_permutations=999,
   )
   discovered = discover_profiles_continuous(colocalization_result=coloc, top_k=3)
   selection = select_profiles(
       X=X, coords=coords, marker_names=marker_names,
       profiles=discovered.profiles,
       _interesting_markers=m1.interesting_markers,
       _colocalization_result=coloc,
       min_spatial_explained=0.90,
       min_protein_explained=0.90,
   )

The output is a set of *unlabeled* marker groups. Assigning a biological cell-type
name to each group is a manual step and produces the profile dict used in Step 3.

Step 2: Initialize and preprocess
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from CITEgeist.model.citegeist_model import CitegeistModel

   model = CitegeistModel(
       sample_name="sample_P1_S1",
       adata=adata,
       output_folder="output/module3",
       simulation=False,
   )

   model.split_adata()                       # required for real patient data
   model.filter_gex(
       nonzero_percentage=0.01,
       mean_expression_threshold=1.1,
       min_counts=100,                       # 25 for low-depth surgical samples
   )
   model.copy_gex_to_protein_adata()
   model.preprocess_gex(target_sum=10000)
   model.preprocess_antibody()

Surgical samples may carry NaN spatial coordinates; filter those explicitly after
``split_adata()``.

Step 3: Load cell-type profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   cell_profiles = {
       "Endothelial":     {"Major": ["PECAM1-1"]},
       "Fibroblasts":     {"Major": ["ACTA2-1"]},
       "B_Cells":         {"Major": ["CD19-1"]},
       "Macrophages":     {"Major": ["CD68-1", "CD163-1"]},
       "CD8_T_Cells":     {"Major": ["CD3E-1", "CD8A-1"]},
       "CD4_T_Cells":     {"Major": ["CD4-1", "CD3E-1"]},
       "Cancer_Luminal":  {"Major": ["EPCAM-1"]},
       "Cancer_Basal":    {"Major": ["KRT5-1", "SDC1-1", "EPCAM-1"]},
   }
   model.load_cell_profile_dict(cell_profiles)

See :doc:`quickstart` for the accepted keys (``Major``, optional ``Soft``) and the
common mistake of passing marker weights instead of marker groups.

Step 4: Spot-level proportions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   global_props, finetuned_props = model.run_cell_proportion_model(
       method="qp",
       use_detection_gating=True,
       use_gex_detection=True,
       gex_detection_k=10,
       validation_warn_only=True,
   )

   model.append_proportions_to_adata(key="finetuned")
   model.get_adata().write_h5ad("output/module3/sample_P1_S1_module3_results.h5ad")

All arguments are keyword-only. ``radius`` defaults to ``None``, which auto-detects
the neighbourhood radius from the spatial coordinates.

Step 5: Single-cell assignment and gene expression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gene-expression deconvolution runs *after* assignment: SACE allocates each spot's
counts to the nuclei assigned within that spot.

.. code-block:: python

   spot_type_gex, cell_adata, diagnostics = model.run_sace_allocation(
       output_mode="single_cell",
   )
   cell_adata.write_h5ad("output/module3/sample_P1_S1_single_cell.h5ad")

``output_mode="single_cell"`` orchestrates StarDist segmentation, per-nucleus
assignment, and SACE projection behind an H&E resolution gate (≤ 1.0 µm/px). It
defaults to ``use_morphology=True``, which scores nucleus morphology
(``CITEgeist.model.morphology.nucleus_scores``) and assigns each nucleus by the
Bayesian posterior. Pass ``use_morphology=False`` for count-constrained Hungarian
assignment, which reproduces the QP spot proportions exactly.

``output_mode="layers"`` (the default) returns spot-level per-type expression
compartments and needs no nuclei or H&E image:

.. code-block:: python

   spot_type_gex, diagnostics = model.run_sace_allocation(output_mode="layers")

Results
-------

* ``{sample_name}_cell_prop_finetuned_results.csv`` — spots × cell-type proportions
* an AnnData with proportion columns appended to ``adata.obs``
* in ``single_cell`` mode, a per-cell AnnData with ``.X`` = allocated counts,
  ``.obs["cell_type"]``, and ``.obsm["spatial"]``

For the runnable scripts, see :doc:`examples`.
