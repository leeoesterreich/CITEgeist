Interactive Tutorials and Vignettes
====================================

This section contains interactive Jupyter notebooks that demonstrate various aspects of CITEgeist analysis. These notebooks provide hands-on examples and can be run interactively in your browser or downloaded locally.

.. raw:: html

    <div class="admonition tip" style="margin: 20px 0;">
    <p class="admonition-title">🚀 Try CITEgeist Interactively!</p>
    <p>Click the <strong>Binder</strong> or <strong>Colab</strong> buttons on any vignette page to run the notebooks instantly in your browser - no installation required!</p>
    <div style="margin: 15px 0;">
        <a href="https://mybinder.org/v2/gh/leeoesterreich/CITEgeist/dev?labpath=docs%2Fcore_scripts%2FJupyter" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <span style="color: #666;">← Launch all notebooks in Binder</span>
    </div>
    </div>

.. toctree::
   :maxdepth: 1
   :caption: Tutorial Notebooks:

   notebooks/vignette_1_biopsy_heterogeneity
   notebooks/vignette_2_surgical_d538g
   notebooks/vignette_3_responder_macrophages
   notebooks/vignette_4_external_tools

Getting Started with Notebooks
------------------------------

To run these notebooks locally:

1. **Install CITEgeist**: Follow the :doc:`installation` guide
2. **Download the notebooks**: Click the download button on any notebook page
3. **Set up your environment**: Ensure you have the required data and dependencies
4. **Run the notebooks**: Execute cells step by step or run all cells

Notebook Descriptions
---------------------

**Vignette 1: Biopsy Heterogeneity Analysis**
   Demonstrates how to use CITEgeist to analyze spatial transcriptomics data from a cancer biopsy sample, exploring cellular heterogeneity and interactions within the tumor microenvironment using Visium spatial data.

**Vignette 2: ESR1 D538G Mutations in Surgical Samples**
   Shows the analysis of spatial transcriptomics data from surgical samples containing ESR1 D538G mutations, combining CITEgeist deconvolution with pathway analysis and mutation signature validation.

**Vignette 3: Macrophage Behavior in Responding Patients**
   Examines macrophage populations and their behavior in responding patients, comparing core biopsies with surgical specimens and focusing on spatial distribution and gene expression patterns of tumor-associated macrophages (TAMs).

**Vignette 4: Integration with External Tools**
   Demonstrates how to integrate CITEgeist results with external analysis tools and workflows for comprehensive spatial transcriptomics analysis.

Data Requirements
-----------------

These notebooks use example datasets that should be available in your CITEgeist installation. Make sure to:

* Update the `DATA_FOLDER` path in each notebook to point to your data
* Configure the `LICENSE_FILE` path for Gurobi
* Ensure all required Python packages are installed

For more information about data preparation, see the :doc:`tutorial` section.
