API Reference
=============

This page provides a comprehensive reference for CITEgeist's public API, organized by workflow and module.

Quick Start
-----------

For most users, the typical workflow involves:

1. **Initialize**: :class:`~CITEgeist.model.citegeist_model.CitegeistModel`
2. **Preprocess**: :meth:`~CITEgeist.model.citegeist_model.CitegeistModel.split_adata`, :meth:`~CITEgeist.model.citegeist_model.CitegeistModel.filter_gex`
3. **Analyze**: :meth:`~CITEgeist.model.citegeist_model.CitegeistModel.run_cell_proportion_model`

Core Model
----------

CitegeistModel
~~~~~~~~~~~~~~

The main analysis class for spatial transcriptomics deconvolution.

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: __init__
   :show-inheritance:

**Initialization & Setup**

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: load_cell_profile_dict, split_adata
   :noindex:

**Data Preprocessing**

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: copy_gex_to_protein_adata, filter_gex, preprocess_antibody, preprocess_gex
   :noindex:

**Core Analysis Methods**

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: assign_cells, compute_expression_prior, run_cell_expression_pass1, run_cell_proportion_model, run_sace_allocation
   :noindex:

**Results & Export**

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: append_gex_to_adata, append_proportions_to_adata, get_adata
   :noindex:

**Validation & Utilities**

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: cleanup, validate_neighborhood_size
   :noindex:

**Static Methods**

.. autoclass:: CITEgeist.model.citegeist_model.CitegeistModel
   :members: global_clr, row_normalize, winsorize
   :noindex:

Configuration
~~~~~~~~~~~~~

.. automodule:: CITEgeist.model.unified_config
   :members:
   :show-inheritance:

Checkpoints
~~~~~~~~~~~

.. automodule:: CITEgeist.model.checkpoints
   :members:
   :show-inheritance:

Deconvolution (Module 3 — QP Proportions)
------------------------------------------

.. automodule:: CITEgeist.model.deconvolution.qp_solver
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.deconvolution.detection
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.deconvolution.detection_refinement
   :members:
   :show-inheritance:

GEX Deconvolution (SACE, Module 3-gex)
---------------------------------------

.. automodule:: CITEgeist.model.gex.sace_gex
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.gex.gex_modules
   :members:
   :show-inheritance:

Cell Assignment (Module 3-post)
--------------------------------

.. automodule:: CITEgeist.model.assignment.cell_assignment
   :members:
   :show-inheritance:

Discovery (Modules 1, 2a/b/c)
------------------------------

.. automodule:: CITEgeist.model.discovery.marker_interest
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.discovery.spatial_colocalization
   :members:
   :show-inheritance:

Programs (Modules 4, 4b, 5)
-----------------------------

.. automodule:: CITEgeist.model.programs.anchored_program_discovery
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.programs.cross_sample_integration
   :members:
   :show-inheritance:

Annotation (Module 3.5)
------------------------

.. automodule:: CITEgeist.model.annotation.functional_annotation
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.annotation.subtype_splitting
   :members:
   :show-inheritance:

Morphology
----------

.. automodule:: CITEgeist.model.morphology.segmentation
   :members:
   :show-inheritance:

Quality Control
---------------

.. automodule:: CITEgeist.model.qc.proportion_qc
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.qc.gex_qc
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.qc.report
   :members:
   :show-inheritance:

Utility Functions
-----------------

.. automodule:: CITEgeist.model.utils
   :members:
   :show-inheritance:

.. automodule:: CITEgeist.model.module2_proposal_builder
   :members:
   :show-inheritance:
