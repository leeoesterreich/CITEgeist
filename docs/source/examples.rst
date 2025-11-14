Examples
========

This section provides practical examples of using CITEgeist for different types of spatial transcriptomics analysis.

.. toctree::
   :maxdepth: 2
   :caption: Interactive Tutorials:

   notebooks

Basic Analysis
--------------

A complete example of CITEgeist analysis:

.. code-block:: python

   import scanpy as sc
   import pandas as pd
   from citegeist_model import CitegeistModel
   
   # Load data
   adata = sc.read_h5ad("spatial_data.h5ad")
   
   # Initialize model
   model = CitegeistModel(
       sample_name="brain_section_1",
       adata=adata,
       output_folder="./citegeist_results"
   )
   
   # Configure Gurobi (if needed)
   model.register_gurobi("/path/to/gurobi.lic")
   
   # Process data
   model.split_adata()
   model.filter_gex()
   model.preprocess_gex()
   model.preprocess_antibody()
   
   # Define cell type profiles
   cell_profiles = {
       "Neuron": {"NeuN": 0.9, "MAP2": 0.8},
       "Astrocyte": {"GFAP": 0.9, "S100B": 0.7},
       "Microglia": {"Iba1": 0.9, "CD68": 0.6},
       "Oligodendrocyte": {"Olig2": 0.8, "MBP": 0.7}
   }
   model.load_cell_profile_dict(cell_profiles)
   
   # Run analysis
   global_props, finetuned_props = model.run_cell_proportion_model(radius=75)
   model.run_cell_expression_pass1(radius=75)
   
   # Get results
   results_adata = model.get_adata()
   print("Analysis complete!")

Parameter Tuning
----------------

Example of parameter optimization:

.. code-block:: python

   # Test different radius values
   for radius in [25, 50, 75, 100]:
       print(f"Testing radius: {radius}")
       
       # Run with different parameters
       props, _ = model.run_cell_proportion_model(
           radius=radius,
           lambda_reg=0.5,
           alpha=0.3
       )
       
       # Evaluate results
       print(f"Mean cell proportions: {props.mean().mean():.3f}")

Batch Processing
----------------

Process multiple samples:

.. code-block:: python

   samples = ["sample1", "sample2", "sample3"]
   
   for sample in samples:
       print(f"Processing {sample}")
       
       # Load sample data
       adata = sc.read_h5ad(f"{sample}.h5ad")
       
       # Initialize model
       model = CitegeistModel(
           sample_name=sample,
           adata=adata,
           output_folder=f"./results/{sample}"
       )
       
       # Run analysis
       # ... (same workflow as above)

Visualization
-------------

Basic visualization of results:

.. code-block:: python

   import matplotlib.pyplot as plt
   
   # Plot cell type proportions
   fig, axes = plt.subplots(2, 2, figsize=(12, 10))
   
   cell_types = ["Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"]
   
   for i, cell_type in enumerate(cell_types):
       ax = axes[i//2, i%2]
       sc.pl.spatial(
           results_adata, 
           color=cell_type, 
           ax=ax,
           title=f"{cell_type} Proportions"
       )
   
   plt.tight_layout()
   plt.show()
