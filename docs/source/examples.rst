Examples
========

.. toctree::
   :maxdepth: 2
   :caption: Interactive Tutorials:

   notebooks

Runnable scripts
----------------

Seven end-to-end runner scripts ship with the repository under ``examples/scripts/``
(and inside the installed package as ``CITEgeist/examples/``). Each takes
command-line arguments; run any of them with ``--help`` to see its options.

.. list-table::
   :header-rows: 1
   :widths: 34 46 20

   * - Script
     - Module(s)
     - GPU
   * - ``run_module12_discovery.py``
     - Modules 1–2: marker interest + protein-profile discovery
     - no
   * - ``run_module3_unified.py``
     - Module 3: cuOPT QP proportions + SACE deconvolved GEX layers
     - **yes**
   * - ``run_cuopt_qp_patient.py``
     - Module 3 proportions only, single patient sample (alternative to the above)
     - **yes**
   * - ``run_single_cell_assignment.py``
     - Module 3-post/3-gex: StarDist nuclei, cell assignment, SACE per-cell GEX
     - Stage 1
   * - ``run_module4_discovery.py``
     - Module 4: anchored spatial gene programs
     - no
   * - ``run_module4b_bivariate.py``
     - Module 4b: bivariate Moran's I between programs
     - no
   * - ``run_module5_integration.py``
     - Module 5: cross-sample program integration
     - no

Input paths in these scripts are placeholders that you must edit for your own data
layout; ``examples/README.md`` documents the expected directory structure and the
order the scripts are meant to be run in.

Worked walkthrough
------------------

For a narrated version of the same pipeline with real API calls, parameter guidance
and a troubleshooting table, see
`docs/quickstart_real_visium.md <https://github.com/leeoesterreich/CITEgeist/blob/main/docs/quickstart_real_visium.md>`_,
or the condensed :doc:`tutorial` page.

Processing several samples
--------------------------

There is no batch entry point; each sample is an independent model instance. On a
cluster, submit one job per sample rather than looping in a single process — Module 3
holds the full antibody matrix and the cuOPT problem in memory.

.. code-block:: python

   import squidpy as sq
   from CITEgeist.model.citegeist_model import CitegeistModel

   def run_one(sample_name, spaceranger_dir, cell_profiles):
       adata = sq.read.visium(
           spaceranger_dir,
           counts_file="filtered_feature_bc_matrix.h5",
           load_images=True,
           gex_only=False,
       )
       model = CitegeistModel(
           sample_name=sample_name,
           adata=adata,
           output_folder=f"output/{sample_name}",
           simulation=False,
       )
       model.split_adata()
       model.filter_gex(min_counts=100)
       model.copy_gex_to_protein_adata()
       model.preprocess_gex()
       model.preprocess_antibody()
       model.load_cell_profile_dict(cell_profiles)
       model.run_cell_proportion_model()
       model.append_proportions_to_adata(key="finetuned")
       model.cleanup()
       return model.get_adata()

Visualizing proportions
-----------------------

``append_proportions_to_adata`` writes one ``adata.obs`` column per cell type, so the
results plot with the standard squidpy spatial plotting call. Always pass an explicit
``crop_coord`` — squidpy does not auto-crop to the spot bounding box.

.. code-block:: python

   import numpy as np
   import squidpy as sq

   result_adata = model.get_adata()

   xy = np.asarray(result_adata.obsm["spatial"])
   crop = (xy[:, 0].min(), xy[:, 1].min(), xy[:, 0].max(), xy[:, 1].max())

   sq.pl.spatial_scatter(
       result_adata,
       color="Macrophages",       # a cell-type column added by append_proportions_to_adata
       crop_coord=crop,
       cmap="YlOrRd",
       alpha=0.7,
   )
