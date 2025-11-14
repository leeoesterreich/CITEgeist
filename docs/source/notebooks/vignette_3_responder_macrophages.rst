Vignette 3: Macrophage Behavior in Responding Patients
======================================================

This vignette examines macrophage populations and their behavior in responding patients, comparing core biopsies with surgical specimens. The analysis focuses on spatial distribution and gene expression patterns of tumor-associated macrophages (TAMs).

Overview
--------

This vignette covers:

* Analysis of macrophage populations in responding patients
* Comparison between core biopsies and surgical specimens
* Spatial distribution analysis of tumor-associated macrophages
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
        <a href="https://mybinder.org/v2/gh/leeoesterreich/dal413/HEAD?labpath=CITEgeist%2FCITEgeist%2FJupyter%2Fvignette_3_responder_macrophages.ipynb" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <a href="https://colab.research.google.com/github/leeoesterreich/dal413/blob/main/CITEgeist/CITEgeist/Jupyter/vignette_3_responder_macrophages.ipynb" target="_blank">
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
* Clinical response annotation
* Required Python packages: scanpy, pandas, matplotlib, seaborn

**Data Requirements**:
* Core biopsy samples from responding patients
* Surgical specimens from the same patients
* Clinical response annotations
* Gene expression and antibody capture data
* Spatial coordinates for all samples

Code Example
------------

Here's a preview of the key steps from the notebook:

.. code-block:: python

   # Define sample paths
   biopsy_responding_paths = [
       os.path.join(DATA_FOLDER, "HCC22-088-P2-S1/outs"),
       os.path.join(DATA_FOLDER, "HCC22-088-P3-S1_A/outs"),
       # ... more biopsy samples
   ]
   
   surg_responding_paths = [
       os.path.join(DATA_FOLDER, "HCC22-088-P2-S2/outs"),
       os.path.join(DATA_FOLDER, "HCC22-088-P3-S2/outs"),
       # ... more surgical samples
   ]
   
   # Process multiple samples
   for sample_path in biopsy_responding_paths:
       model = CitegeistModel(
           sample_name=f"biopsy_{sample_path.split('/')[-3]}",
           adata=load_sample(sample_path),
           output_folder="./macrophage_analysis"
       )
       # ... run analysis

Expected Outputs
----------------

This notebook will generate:

* Macrophage proportion comparisons between sample types
* Spatial distribution maps of macrophage populations
* Gene expression profiles of macrophage subtypes
* Clinical response correlation analysis
* Comparative visualizations between biopsies and surgical samples

Next Steps
----------

After completing this vignette:

1. Apply the analysis to your own multi-sample datasets
2. Explore other immune cell populations
3. Integrate with additional clinical variables
4. Try the external tools integration in :doc:`vignette_4_external_tools`

Related Documentation
---------------------

* :doc:`../tutorial` - Comprehensive workflow guide
* :doc:`vignette_2_surgical_d538g` - Surgical sample analysis
* :doc:`../api` - Complete API reference

