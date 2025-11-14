Quick Start
===========

This guide will help you get started with CITEgeist for spatial transcriptomics analysis.

Basic Usage
-----------

Here's a minimal example of how to use CITEgeist:

.. code-block:: python

   import scanpy as sc
   from citegeist_model import CitegeistModel
   
   # Load your spatial transcriptomics data
   adata = sc.read_h5ad("your_data.h5ad")
   
   # Initialize CITEgeist model
   model = CitegeistModel(
       sample_name="sample_1",
       adata=adata,
       output_folder="./results"
   )
   
   # Split data into gene expression and antibody capture
   model.split_adata()
   
   # Load cell type profiles
   cell_profiles = {
       "T_cell": {...},  # Your cell type profiles
       "B_cell": {...},
       # ... more cell types
   }
   model.load_cell_profile_dict(cell_profiles)
   
   # Preprocess data
   model.filter_gex()
   model.preprocess_gex()
   model.preprocess_antibody()
   
   # Run cell proportion optimization
   global_props, finetuned_props = model.run_cell_proportion_model(radius=50)
   
   # Run gene expression deconvolution
   model.run_cell_expression_pass1(radius=50)

Key Parameters
--------------

* **radius**: Spatial neighborhood radius for optimization
* **lambda_reg**: Regularization strength for cell proportions
* **alpha**: L1-L2 tradeoff factor (0 = L2, 1 = L1)
* **max_workers**: Number of parallel workers for optimization

For more detailed examples, see the :doc:`tutorial` section.
