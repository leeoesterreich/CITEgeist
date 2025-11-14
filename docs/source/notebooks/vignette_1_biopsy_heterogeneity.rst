Vignette 1: Biopsy Heterogeneity Analysis
=========================================

This notebook demonstrates how to use CITEgeist to analyze spatial transcriptomics data from a cancer biopsy sample. We'll explore cellular heterogeneity and interactions within the tumor microenvironment using Visium spatial data.

Overview
--------

This vignette covers:

* Loading and preprocessing Visium spatial transcriptomics data
* Setting up CITEgeist for biopsy analysis
* Exploring cellular heterogeneity in cancer samples
* Analyzing spatial patterns of cell type distributions
* Visualizing results and interpreting findings

Key Learning Objectives
-----------------------

After completing this notebook, you will understand:

* How to prepare Visium data for CITEgeist analysis
* The importance of radius parameter selection for spatial neighborhoods
* How to interpret cell type proportion results
* Methods for visualizing spatial transcriptomics results

Run Interactive Notebook
-------------------------

.. raw:: html

    <div style="margin: 20px 0;">
        <a href="https://mybinder.org/v2/gh/leeoesterreich/dal413/HEAD?labpath=CITEgeist%2FCITEgeist%2FJupyter%2Fvignette_1_biopsy_heterogeneity.ipynb" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <a href="https://colab.research.google.com/github/leeoesterreich/dal413/blob/main/CITEgeist/CITEgeist/Jupyter/vignette_1_biopsy_heterogeneity.ipynb" target="_blank">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
        </a>
    </div>

    <div class="admonition tip">
    <p class="admonition-title">🚀 Run Interactively</p>
    <p><strong>Binder</strong>: Click above to launch this notebook in a live, executable environment with all dependencies pre-installed. No setup required!</p>
    <p><strong>Colab</strong>: Run in Google Colab with free GPU access. You may need to install some packages.</p>
    <p><strong>Download</strong>: <a href="vignette_1_biopsy_heterogeneity.ipynb">Download the notebook</a> to run locally.</p>
    </div>

**Prerequisites**:
* CITEgeist installed (see :doc:`../installation`)
* Gurobi license configured
* Visium spatial transcriptomics data
* Required Python packages: scanpy, squidpy, matplotlib, seaborn

**Data Requirements**:
* Visium output folder with spatial coordinates
* Gene expression and antibody capture data
* Properly formatted feature types

Code Example
------------

Here's a preview of the key steps from the notebook:

.. code-block:: python

   # Load dataset
   path_to_visium_folder = os.path.join(DATA_FOLDER, "HCC22-088-P1-S2/outs")
   
   # Initialize CITEgeist model
   model = CitegeistModel(
       sample_name="biopsy_sample",
       adata=adata,
       output_folder="./results"
   )
   
   # Split and preprocess data
   model.split_adata()
   model.filter_gex()
   model.preprocess_gex()
   model.preprocess_antibody()
   
   # Run analysis
   global_props, finetuned_props = model.run_cell_proportion_model(radius=50)

Expected Outputs
----------------

This notebook will generate:

* Cell type proportion estimates for each spot
* Spatial maps showing cell type distributions
* Quality metrics and diagnostic plots
* Processed data files for downstream analysis

Next Steps
----------

After completing this vignette:

1. Try the notebook with your own Visium data
2. Experiment with different radius parameters
3. Explore the other vignettes for advanced analyses
4. Check out the :doc:`../tutorial` for detailed explanations

Related Documentation
---------------------

* :doc:`../tutorial` - Comprehensive workflow guide
* :doc:`../examples` - Additional code examples
* :doc:`../api` - Complete API reference

