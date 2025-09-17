Changelog
=========

All notable changes to CITEgeist will be documented in this file.

Version 1.0.0 (Unreleased)
--------------------------

* Initial release of CITEgeist
* Spatial transcriptomics deconvolution framework
* Cell type proportion estimation
* Gene expression profile deconvolution
* Gurobi-based optimization
* Comprehensive documentation

Features
~~~~~~~~

* **CitegeistModel**: Main class for spatial transcriptomics analysis
* **Cell Proportion Optimization**: Estimate cell type proportions using antibody capture data
* **Gene Expression Deconvolution**: Estimate cell-type-specific gene expression profiles
* **Spatial Regularization**: Incorporate spatial neighborhood information
* **Checkpointing**: Save and resume long-running optimization processes
* **Parallel Processing**: Multi-core optimization support

API Reference
~~~~~~~~~~~~~

* `citegeist_model.CitegeistModel`: Main model class
* `gurobi_impl`: Core optimization functions
* `utils`: Helper functions and utilities
* `checkpoints.CheckpointManager`: Checkpoint management

Dependencies
~~~~~~~~~~~~

* Python 3.7+
* numpy
* pandas
* scanpy
* scipy
* pyarrow
* gurobipy (requires Gurobi license)
