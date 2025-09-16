Vignette 3: Macrophage Behavior in Responding Patients
======================================================

This vignette examines macrophage populations and their behavior in responding patients, comparing core biopsies with surgical specimens and focusing on spatial distribution and gene expression patterns of tumor-associated macrophages (TAMs).

Overview
--------

This vignette covers:

* Analysis of macrophage populations in treatment-responsive patients
* Comparison between core biopsy and surgical specimen samples
* Spatial distribution analysis of tumor-associated macrophages (TAMs)
* Gene expression profiling of macrophage subtypes
* Integration with clinical response data

Key Learning Objectives
-----------------------

After completing this notebook, you will understand:

* How to analyze macrophage populations with CITEgeist
* Methods for comparing different sample types (biopsy vs surgical)
* Approaches to spatial analysis of immune cell populations
* Techniques for integrating clinical response data

Run Interactive Notebook
-------------------------

.. raw:: html

    <div style="margin: 20px 0;">
        <a href="https://mybinder.org/v2/gh/leeoesterreich/CITEgeist/dev?labpath=docs%2Fcore_scripts%2FJupyter%2Fvignette_3_responder_macrophages.ipynb" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <a href="https://colab.research.google.com/github/leeoesterreich/CITEgeist/blob/dev/docs/core_scripts/Jupyter/vignette_3_responder_macrophages.ipynb" target="_blank">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
        </a>
    </div>

    <div class="admonition tip">
    <p class="admonition-title">🚀 Run Interactively</p>
    <p><strong>Binder</strong>: Click above to launch this notebook in a live, executable environment with all dependencies pre-installed. No setup required!</p>
    <p><strong>Colab</strong>: Run in Google Colab with free GPU access. You may need to install some packages.</p>
    <p><strong>Download</strong>: <a href="vignette_3_responder_macrophages.ipynb">Download the notebook</a> to run locally.</p>
    </div>

**Prerequisites**:
* CITEgeist installed (see :doc:`../installation`)
* Gurobi license configured
* Multiple sample types (biopsy and surgical)
* Clinical response data
* Required Python packages: scanpy, squidpy, matplotlib, seaborn, pandas

**Data Requirements**:
* Paired biopsy and surgical specimens from responding patients
* Clinical response annotations
* Gene expression and antibody capture data
* Spatial coordinates for both sample types

Code Example
------------

Here's a preview of the comparative analysis workflow:

.. code-block:: python

   # Load paired samples
   biopsy_adata = sc.read_visium(path_to_biopsy)
   surgical_adata = sc.read_visium(path_to_surgical)
   
   # Initialize CITEgeist models for comparison
   biopsy_model = CitegeistModel(
       sample_name="responder_biopsy",
       adata=biopsy_adata,
       output_folder="./biopsy_results"
   )
   
   surgical_model = CitegeistModel(
       sample_name="responder_surgical", 
       adata=surgical_adata,
       output_folder="./surgical_results"
   )
   
   # Process both samples
   for model in [biopsy_model, surgical_model]:
       model.split_adata()
       model.filter_gex()
       model.preprocess_gex()
       model.preprocess_antibody()
   
   # Run deconvolution with macrophage focus
   biopsy_props = biopsy_model.run_cell_proportion_model(radius=60)
   surgical_props = surgical_model.run_cell_proportion_model(radius=60)
   
   # Compare macrophage populations
   compare_macrophage_populations(biopsy_props, surgical_props, response_data)

Expected Outputs
----------------

This notebook will generate:

* Comparative analysis of macrophage populations between sample types
* Spatial distribution maps of TAM subtypes
* Gene expression profiles of responding vs non-responding regions
* Statistical comparisons of immune cell infiltration
* Predictive models for treatment response based on macrophage patterns

Advanced Analysis
-----------------

This vignette includes:

* **Multi-sample Comparison**: Systematic comparison between biopsy and surgical specimens
* **Immune Cell Profiling**: Detailed analysis of macrophage subtypes (M1/M2 polarization)
* **Spatial Statistics**: Quantification of immune cell spatial relationships
* **Response Correlation**: Integration with clinical response data
* **Predictive Modeling**: Building models to predict treatment response

Clinical Relevance
------------------

Key findings from this analysis:

* **Treatment Impact**: How therapy affects macrophage populations
* **Spatial Reorganization**: Changes in immune cell distribution post-treatment
* **Biomarker Discovery**: Identification of predictive immune signatures
* **Therapeutic Targets**: Potential targets for immunomodulatory therapy

Next Steps
----------

After completing this vignette:

1. Apply to your own paired sample datasets
2. Explore other immune cell populations beyond macrophages
3. Integrate with additional clinical variables
4. Consider longitudinal analysis with multiple timepoints

Related Documentation
---------------------

* :doc:`vignette_1_biopsy_heterogeneity` - Basic spatial analysis workflow
* :doc:`vignette_2_surgical_d538g` - Mutation-focused analysis
* :doc:`vignette_4_external_tools` - Integration with other tools
* :doc:`../tutorial` - Comprehensive methodology guide
* :doc:`../api` - Complete API reference
