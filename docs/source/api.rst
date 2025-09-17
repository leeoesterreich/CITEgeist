API Reference
=============

This page provides a comprehensive reference for CITEgeist's public API, organized by workflow and functionality.

Quick Start
-----------

For most users, the typical workflow involves:

1. **Initialize**: :class:`~citegeist_model.CitegeistModel`
2. **Preprocess**: :meth:`~citegeist_model.CitegeistModel.split_adata`, :meth:`~citegeist_model.CitegeistModel.filter_gex`
3. **Analyze**: :meth:`~citegeist_model.CitegeistModel.run_cell_proportion_model`

Core Classes
------------

CitegeistModel
~~~~~~~~~~~~~~

The main analysis class for spatial transcriptomics deconvolution.

.. autoclass:: citegeist_model.CitegeistModel
   :members: __init__
   :show-inheritance:

**Initialization & Setup**

.. autoclass:: citegeist_model.CitegeistModel
   :members: load_cell_profile_dict, register_gurobi, split_adata
   :noindex:

**Data Preprocessing** 

.. autoclass:: citegeist_model.CitegeistModel  
   :members: copy_gex_to_protein_adata, filter_gex, preprocess_antibody, preprocess_gex
   :noindex:

**Core Analysis Methods**

.. autoclass:: citegeist_model.CitegeistModel
   :members: compute_expression_prior, run_cell_expression_pass1, run_cell_proportion_model
   :noindex:

**Results & Export**

.. autoclass:: citegeist_model.CitegeistModel
   :members: append_gex_to_adata, append_proportions_to_adata, get_adata
   :noindex:

**Validation & Utilities**

.. autoclass:: citegeist_model.CitegeistModel
   :members: cleanup, validate_neighborhood_size
   :noindex:

**Static Methods**

.. autoclass:: citegeist_model.CitegeistModel
   :members: global_clr, row_normalize, winsorize
   :noindex:

Optimization Engine
-------------------

Low-level optimization functions (advanced users).

**Cell Proportion Optimization**

.. automodule:: gurobi_impl
   :members: deconvolute_local_cell_proportions, finetune_cell_proportions, optimize_cell_proportions
   :show-inheritance:

**Gene Expression Deconvolution**

.. automodule:: gurobi_impl  
   :members: deconvolute_spot_with_neighbors_with_prior, optimize_gene_expression
   :noindex:

**Prior Computation**

.. automodule:: gurobi_impl
   :members: compute_global_prior, validate_prior_effect
   :noindex:

**Data Processing**

.. automodule:: gurobi_impl
   :members: map_antibodies_to_profiles, normalize_counts, scale_genes, unscale_genes
   :noindex:

**Analysis & Logging**

.. automodule:: gurobi_impl
   :members: log_marker_gene_patterns
   :noindex:

Utility Functions
-----------------

Helper functions for analysis and validation.

**Neighborhood Analysis**

.. automodule:: utils
   :members: assert_neighborhood_size, find_fixed_radius_neighbors, get_neighbors_with_fixed_radius, plot_neighbors_with_fixed_radius
   :show-inheritance:

**Performance Evaluation**

.. automodule:: utils
   :members: benchmark_cell_proportions, calculate_expression_metrics
   :noindex:

**Data Export & Management**

.. automodule:: utils
   :members: export_anndata_layers, save_results_to_output
   :noindex:

**System Utilities**

.. automodule:: utils
   :members: cleanup_memory, setup_logging, validate_cell_profile_dict
   :noindex:


