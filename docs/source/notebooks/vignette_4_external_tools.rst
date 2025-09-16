Vignette 4: Integration with External Tools
===========================================

This vignette demonstrates how to integrate CITEgeist results with external analysis tools and workflows for comprehensive spatial transcriptomics analysis. Learn to export data, interface with popular tools, and build integrated analysis pipelines.

Overview
--------

This vignette covers:

* Exporting CITEgeist results in various formats
* Integration with external spatial analysis tools
* Data format requirements for different platforms
* Building comprehensive analysis workflows
* Best practices for tool interoperability

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
        <a href="https://mybinder.org/v2/gh/leeoesterreich/CITEgeist/dev?labpath=docs%2Fcore_scripts%2FJupyter%2Fvignette_4_external_tools.ipynb" target="_blank">
            <img src="https://mybinder.org/badge_logo.svg" alt="Launch Binder" style="margin-right: 10px;">
        </a>
        <a href="https://colab.research.google.com/github/leeoesterreich/CITEgeist/blob/dev/docs/core_scripts/Jupyter/vignette_4_external_tools.ipynb" target="_blank">
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
* Required Python packages: scanpy, squidpy, pandas, numpy

**Data Requirements**:
* CITEgeist analysis results from previous vignettes
* Spatial transcriptomics data with coordinates
* Cell type proportion estimates
* Gene expression profiles

Code Example
------------

Here's a preview of the integration workflow:

.. code-block:: python

   # Export CITEgeist results for external tools
   from model.utils import export_anndata_layers
   
   # Load CITEgeist results
   model = CitegeistModel(sample_name="integration_example")
   model.load_results("./citegeist_results/")
   
   # Export in various formats
   export_anndata_layers(model.get_adata(), "./exports/", pass_number=2)
   
   # Export for specific tools
   export_for_commot(model.get_adata(), "./commot_input/")
   export_for_cellphonedb(model.get_adata(), "./cellphonedb_input/")
   export_for_squidpy(model.get_adata(), "./squidpy_input/")
   
   # Integration with COMMOT for cell communication
   import commot as ct
   adata_commot = ct.read_spatial_data("./commot_input/")
   ct.tl.spatial_communication(adata_commot)
   
   # Integration with CellPhoneDB for ligand-receptor analysis
   run_cellphonedb_analysis("./cellphonedb_input/")
   
   # Combine results back into comprehensive analysis
   integrated_results = combine_analysis_results(
       citegeist_results=model.get_adata(),
       commot_results=adata_commot,
       cellphonedb_results="./cellphonedb_output/"
   )

External Tool Integration
-------------------------

**Supported Integrations**:

* **COMMOT**: Cell communication and signaling pathway analysis
* **CellPhoneDB**: Ligand-receptor interaction analysis  
* **Squidpy**: Additional spatial statistics and neighborhood analysis
* **Seurat**: R-based spatial transcriptomics workflows
* **STUtility**: Spatial transcriptomics utilities and visualization
* **SpaceRanger**: 10x Genomics spatial pipeline integration

**Export Formats**:

* **AnnData**: Standard scanpy format with layers
* **CSV**: Cell type proportions and gene expression matrices
* **H5AD**: Compressed AnnData format
* **MTX**: Matrix Market format for sparse matrices
* **Parquet**: High-performance columnar format

Advanced Workflows
------------------

This vignette demonstrates:

* **Multi-tool Pipelines**: Chaining different analysis tools
* **Format Conversion**: Converting between different data formats
* **Result Integration**: Combining outputs from multiple tools
* **Quality Control**: Validating integration results
* **Visualization**: Creating comprehensive plots from integrated data

Integration Examples
--------------------

**Example 1: CITEgeist + COMMOT**
```python
# Use CITEgeist for deconvolution, COMMOT for cell communication
cell_proportions = citegeist_deconvolution(adata)
communication_scores = commot_analysis(adata, cell_proportions)
```

**Example 2: CITEgeist + CellPhoneDB**
```python
# Combine spatial deconvolution with ligand-receptor analysis
spatial_context = citegeist_spatial_analysis(adata)
lr_interactions = cellphonedb_analysis(spatial_context)
```

**Example 3: Multi-tool Validation**
```python
# Cross-validate results using multiple approaches
citegeist_results = run_citegeist_analysis(adata)
squidpy_results = run_squidpy_analysis(adata)
validate_consistency(citegeist_results, squidpy_results)
```

Expected Outputs
----------------

This notebook will generate:

* Exported data files in multiple formats
* Integration results from external tools
* Comparative analysis plots
* Workflow documentation and best practices
* Validation metrics for integration quality

Best Practices
--------------

**Data Export Guidelines**:
* Always validate exported data integrity
* Document format specifications and requirements
* Include metadata and analysis parameters
* Version control for reproducibility

**Integration Strategies**:
* Start with simple pairwise integrations
* Validate results at each integration step
* Document tool versions and parameters
* Create reproducible workflows

Next Steps
----------

After completing this vignette:

1. Explore additional external tools in your domain
2. Develop custom integration workflows
3. Contribute integration methods back to the community
4. Build automated analysis pipelines

Related Documentation
---------------------

* :doc:`vignette_1_biopsy_heterogeneity` - Basic CITEgeist workflow
* :doc:`vignette_2_surgical_d538g` - Advanced analysis techniques
* :doc:`vignette_3_responder_macrophages` - Multi-sample analysis
* :doc:`../tutorial` - Comprehensive methodology guide
* :doc:`../api` - Complete API reference

External Resources
------------------

* `COMMOT Documentation <https://commot.readthedocs.io/>`_
* `CellPhoneDB Documentation <https://cellphonedb.readthedocs.io/>`_
* `Squidpy Documentation <https://squidpy.readthedocs.io/>`_
* `Seurat Spatial Vignettes <https://satijalab.org/seurat/articles/spatial_vignette.html>`_
