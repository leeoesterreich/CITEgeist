Tutorial
========

This tutorial provides a comprehensive guide to using CITEgeist for spatial transcriptomics analysis.

Data Preparation
----------------

CITEgeist requires spatial transcriptomics data with both gene expression and antibody capture information. Your data should be in AnnData format with:

* Gene expression counts in the main matrix
* Antibody capture data in a separate layer or as part of the main matrix
* Spatial coordinates in `adata.obsm['spatial']`
* Feature types specified in `adata.var['feature_types']`

Workflow Overview
-----------------

The typical CITEgeist workflow consists of:

1. **Data Loading and Splitting**: Separate gene expression and antibody capture data
2. **Preprocessing**: Filter and normalize both data types
3. **Cell Proportion Optimization**: Estimate cell type proportions using antibody data
4. **Gene Expression Deconvolution**: Estimate cell-type-specific gene expression profiles

Detailed Steps
--------------

Step 1: Initialize the Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   model = CitegeistModel(
       sample_name="your_sample",
       adata=your_adata,
       output_folder="./results"
   )

Step 2: Split and Preprocess Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Split data
   model.split_adata()
   
   # Filter gene expression
   model.filter_gex(nonzero_percentage=0.01, mean_expression_threshold=1.1)
   
   # Normalize data
   model.preprocess_gex(target_sum=10000)
   model.preprocess_antibody()

Step 3: Load Cell Type Profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   cell_profiles = {
       "T_cell": {"CD3": 0.8, "CD4": 0.6, ...},
       "B_cell": {"CD19": 0.9, "CD20": 0.7, ...},
       # ... define profiles for each cell type
   }
   model.load_cell_profile_dict(cell_profiles)

Step 4: Run Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Cell proportion optimization
   global_props, finetuned_props = model.run_cell_proportion_model(
       radius=50,
       lambda_reg=1.0,
       alpha=0.5
   )
   
   # Gene expression deconvolution
   model.run_cell_expression_pass1(radius=50)

Results
-------

CITEgeist outputs:

* Cell type proportion estimates for each spot
* Cell-type-specific gene expression profiles
* Spatial maps of cell type distributions
* Quality metrics and diagnostic information

For more examples, see the :doc:`examples` section.
