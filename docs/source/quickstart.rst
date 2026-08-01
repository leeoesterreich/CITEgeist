Quick Start
===========

This page is the shortest path from a SpaceRanger output directory to spot-level
cell-type proportions. For the full end-to-end walkthrough — profile discovery,
proportions, and single-cell output — see
`Running on Real Visium Data <https://github.com/leeoesterreich/CITEgeist/blob/main/docs/quickstart_real_visium.md>`_.

.. warning::

   **A GPU is required.**
   :meth:`~CITEgeist.model.citegeist_model.CitegeistModel.run_cell_proportion_model`
   solves the proportion QP with NVIDIA cuOPT. cuOPT is **not on PyPI and is not
   included in** ``CITEgeist_env.yml`` — install it separately
   (``pip install cuopt-cu12 --extra-index-url https://pypi.nvidia.com``) or use a
   RAPIDS/NGC container. **There is no CPU fallback**; on a CPU-only machine the
   proportion solver cannot run. See :doc:`installation`.

No external scRNA-seq reference is required — cell-type profiles are built from the
antibody panel measured on the same slide.

Minimal example
---------------

.. code-block:: python

   import squidpy as sq
   from CITEgeist import CitegeistModel

   # Load SpaceRanger output. gex_only=False keeps the antibody-capture features.
   adata = sq.read.visium(
       "/path/to/spaceranger/output/sample_P1_S1",
       counts_file="filtered_feature_bc_matrix.h5",
       load_images=True,
       gex_only=False,
   )

   model = CitegeistModel(
       sample_name="sample_P1_S1",
       adata=adata,
       output_folder="output/module3",
       simulation=False,
   )

   # Real (non-simulated) data must be split before any preprocessing.
   model.split_adata()
   model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1, min_counts=100)
   model.copy_gex_to_protein_adata()
   model.preprocess_gex(target_sum=10000)
   model.preprocess_antibody()

   # Cell-type profiles map a cell-type name to its marker groups.
   # Antibody names in real SpaceRanger output carry a "-1" suffix.
   cell_profiles = {
       "Endothelial":    {"Major": ["PECAM1-1"]},
       "Macrophages":    {"Major": ["CD68-1", "CD163-1"]},
       "CD8_T_Cells":    {"Major": ["CD3E-1", "CD8A-1"]},
       "Cancer_Luminal": {"Major": ["EPCAM-1"]},
   }
   model.load_cell_profile_dict(cell_profiles)

   # Spot-level proportions (cuOPT QP — GPU only).
   global_props, finetuned_props = model.run_cell_proportion_model()

   model.append_proportions_to_adata(key="finetuned")
   result_adata = model.get_adata()

The finetuned proportions are also written to
``{sample_name}_cell_prop_finetuned_results.csv`` under ``output_folder`` and kept on
the model as ``model.results["cell_prop"]`` for downstream steps.

Cell-profile format
-------------------

``load_cell_profile_dict`` expects a mapping of cell-type name to a dict of marker
*groups* — not a marker-to-weight mapping:

.. code-block:: python

   {
       "CD8_T_Cells":  {"Major": ["CD3E-1", "CD8A-1"]},
       "Cancer_Basal": {"Major": ["KRT5-1", "SDC1-1", "EPCAM-1"]},
   }

The solver reads two keys. ``"Major"`` is a list of discriminative marker names.
``"Soft"`` is an optional list of ``(marker_name, weight)`` tuples that receive
graded rather than binary assignment weights. A legacy ``"Minor"`` key is accepted
by the loader but ignored by the solver.

``validate_cell_profile_dict`` only checks the outer ``dict[str, dict]`` shape, so a
marker-to-weight mapping passes validation and then fails deep inside the solver.
Use the group format.

Key parameters
--------------

All arguments to
:meth:`~CITEgeist.model.citegeist_model.CitegeistModel.run_cell_proportion_model`
are keyword-only.

* ``radius`` — spatial neighbourhood radius. ``None`` (the default) auto-detects it
  from the spatial coordinates (three spot rings).
* ``lambda_reg`` — regularization strength on the proportion vector (default ``1``).
* ``alpha`` — L1/L2 trade-off, ``0`` = L2, ``1`` = L1 (default ``0.5``).
* ``use_gex_detection`` — GEX-informed detection refinement (default ``True``).
* ``method`` — ``"qp"``, the only production backend.

Next steps
----------

* :doc:`tutorial` — the same workflow explained step by step, through per-cell
  gene-expression deconvolution.
* :doc:`examples` — the runnable scripts that ship with the package.
