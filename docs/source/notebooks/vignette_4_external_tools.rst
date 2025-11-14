Vignette 4: Integration with External Tools
============================================

This vignette demonstrates how to integrate CITEgeist results with external analysis tools and workflows for comprehensive spatial transcriptomics analysis. It shows how to export CITEgeist outputs and use them with other popular bioinformatics tools.

Overview
--------

This vignette covers:

* Exporting CITEgeist results for external analysis
* Integration with communication analysis tools (COMMOT)
* Data format conversion and compatibility
* Workflow integration strategies
* Quality control and validation

Key Learning Objectives
-----------------------

After completing this notebook, you will understand:

* How to export CITEgeist results in various formats
* Methods for integrating with external analysis tools
* Data format requirements for different tools
* Strategies for building comprehensive analysis workflows

Run Interactive Notebook
-------------------------

.. raw:: html

    <div style="margin: 20px 0;">
        <a href="https://mybinder.org/v2/gh/leeoesterreich/dal413/HEAD?labpath=CITEgeist%2FCITEgeist%2FJupyter%2Fvignette_4_external_tools.ipynb" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <a href="https://colab.research.google.com/github/leeoesterreich/dal413/blob/main/CITEgeist/CITEgeist/Jupyter/vignette_4_external_tools.ipynb" target="_blank">
            <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab">
        </a>
    </div>

    <div class="admonition tip">
    <p class="admonition-title">🚀 Run Interactively</p>
    <p><strong>Binder</strong>: Click above to launch this notebook in a live, executable environment with all dependencies pre-installed. No setup required!</p>
    <p><strong>Colab</strong>: Run in Google Colab with free GPU access. You may need to install some packages.</p>
    <p><strong>Download</strong>: <a href="vignette_4_external_tools.ipynb">Download the notebook</a> to run locally.</p>
    </div>

**Prerequisites**:
* CITEgeist installed (see :doc:`../installation`)
* Gurobi license configured
* External analysis tools (COMMOT, etc.)
* Required Python packages: scanpy, squidpy, pandas, matplotlib

**Data Requirements**:
* Processed CITEgeist results
* Spatial transcriptomics data
* Cell type annotations
* Communication ligand-receptor databases

Code Example
------------

Here's a preview of the key steps from the notebook:

.. code-block:: python

   # Load CITEgeist results
   model = CitegeistModel(
       sample_name="integration_example",
       adata=processed_adata,
       output_folder="./integration_results"
   )
   
   # Export results for external tools
   results_adata = model.get_adata()
   
   # Prepare data for COMMOT analysis
   # Merge cell types and filter
   merged_adata = merge_celltypes(results_adata, threshold=0.2)
   
   # Run communication analysis
   sq.gr.ligrec(merged_adata, n_perms=1000)
   sq.gr.ligrec(merged_adata, cluster_key="cell_type")

Expected Outputs
----------------

This notebook will generate:

* Exported CITEgeist results in multiple formats
* Communication analysis results
* Integrated analysis visualizations
* Quality control reports
* Workflow documentation

Next Steps
----------

After completing this vignette:

1. Explore other external tool integrations
2. Build custom analysis pipelines
3. Develop automated workflows
4. Share your integration strategies with the community

Related Documentation
---------------------

* :doc:`../tutorial` - Comprehensive workflow guide
* :doc:`vignette_3_responder_macrophages` - Multi-sample analysis
* :doc:`../api` - Complete API reference
* :doc:`../contributing` - Guidelines for sharing workflows

