Vignette 2: ESR1 D538G Mutations in Surgical Samples
====================================================

This vignette demonstrates the analysis of spatial transcriptomics data from surgical samples containing ESR1 D538G mutations, a common mutation in breast cancer. The analysis combines CITEgeist deconvolution with pathway analysis and mutation signature validation.

Overview
--------

This vignette covers:

* Analysis of surgical specimens with known mutations
* Integration of mutation data with spatial transcriptomics
* Pathway analysis and signature validation
* Comparison of mutation-positive and mutation-negative regions
* Advanced visualization techniques

Key Learning Objectives
-----------------------

After completing this notebook, you will understand:

* How to analyze surgical specimens with CITEgeist
* Methods for integrating mutation data with spatial analysis
* Approaches to pathway analysis in spatial context
* Techniques for validating mutation signatures spatially

Run Interactive Notebook
-------------------------

.. raw:: html

    <div style="margin: 20px 0;">
        <a href="https://mybinder.org/v2/gh/leeoesterreich/CITEgeist/dev?labpath=docs%2Fcore_scripts%2FJupyter%2Fvignette_2_surgical_d538g.ipynb" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <a href="https://colab.research.google.com/github/leeoesterreich/CITEgeist/blob/dev/docs/core_scripts/Jupyter/vignette_2_surgical_d538g.ipynb" target="_blank">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
        </a>
    </div>

    <div class="admonition tip">
    <p class="admonition-title">🚀 Run Interactively</p>
    <p><strong>Binder</strong>: Click above to launch this notebook in a live, executable environment with all dependencies pre-installed. No setup required!</p>
    <p><strong>Colab</strong>: Run in Google Colab with free GPU access. You may need to install some packages.</p>
    <p><strong>Download</strong>: <a href="vignette_2_surgical_d538g.ipynb">Download the notebook</a> to run locally.</p>
    </div>

**Prerequisites**:
* CITEgeist installed (see :doc:`../installation`)
* Gurobi license configured
* Surgical sample spatial transcriptomics data
* Mutation annotation data
* Required Python packages: scanpy, pandas, matplotlib, seaborn

**Data Requirements**:
* Spatial transcriptomics data from surgical specimens
* Mutation annotation (ESR1 D538G status)
* Gene expression and antibody capture data
* Spatial coordinates and tissue morphology

Code Example
------------

Here's a preview of the key analysis steps:

.. code-block:: python

   # Load surgical sample data with mutation annotations
   adata = sc.read_visium(path_to_surgical_sample)
   
   # Initialize CITEgeist model
   model = CitegeistModel(
       sample_name="surgical_d538g",
       adata=adata,
       output_folder="./surgical_results"
   )
   
   # Preprocess and analyze
   model.split_adata()
   model.filter_gex()
   model.preprocess_gex()
   model.preprocess_antibody()
   
   # Run deconvolution with mutation context
   global_props, finetuned_props = model.run_cell_proportion_model(radius=75)
   
   # Integrate with mutation data for pathway analysis
   mutation_enriched_spots = identify_mutation_regions(adata, mutation_data)

Expected Outputs
----------------

This notebook will generate:

* Cell type proportions stratified by mutation status
* Spatial maps highlighting mutation-positive regions
* Pathway enrichment analysis results
* Differential expression between mutation contexts
* Validation plots for mutation signatures

Advanced Analysis
-----------------

This vignette includes:

* **Mutation Integration**: Methods for incorporating genetic variant data
* **Pathway Analysis**: Gene set enrichment in spatial context
* **Comparative Analysis**: Mutation-positive vs negative region comparison
* **Signature Validation**: Confirming mutation-associated expression patterns

Next Steps
----------

After completing this vignette:

1. Apply the workflow to your own surgical specimens
2. Explore different mutation types and their spatial patterns
3. Integrate with clinical outcome data
4. Consider multi-sample comparative analyses

Related Documentation
---------------------

* :doc:`vignette_1_biopsy_heterogeneity` - Basic CITEgeist workflow
* :doc:`vignette_3_responder_macrophages` - Treatment response analysis
* :doc:`../tutorial` - Comprehensive methodology guide
* :doc:`../api` - Complete API reference
