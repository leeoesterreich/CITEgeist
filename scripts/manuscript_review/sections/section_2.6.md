# Section 2.6: Interpretable Outputs Enable Downstream Analysis

**Referenced Figures**: Figure 6 (A, B, C)

## Section Text

A key design principle of CITEgeist is generating outputs compatible with established bioinformatics tools. All intermediate and final results are stored in standard formats (AnnData, CSV, JSON), enabling seamless integration with the broader computational biology ecosystem (Figure 6A).

**Differential expression.** Cell type-specific expression layers from Module 3 can be directly analyzed using PyDESeq2 or other DE frameworks. By filtering spots based on cell type proportions, researchers can perform DE analysis within specific populations—for example, comparing macrophage-enriched spots between conditions. The Module 5 response analysis demonstrated this workflow, identifying 127 DE genes between responder and progressor samples using standard PyDESeq2.

**Pathway analysis.** Gene lists from DE analysis or program loadings from Module 4 integrate directly with GSEApy/Enrichr for pathway enrichment. In our cohort, progressor-enriched genes showed enrichment for extracellular matrix organization (GO:0030198), collagen catabolic process (GO:0030574), and the MSigDB Hallmark epithelial-mesenchymal transition signature.

**Cell-cell communication.** CITEgeist's spatial neighbor graphs and cell type proportions are compatible with tools like COMMOT and CellPhoneDB for inferring ligand-receptor interactions. Spatial coordinates and cell type stratification enable analysis of distance-dependent signaling, such as identifying ligand-receptor pairs enriched at specific cellular interfaces (Figure 6B).

**Clustering and visualization.** Standard scanpy workflows (Leiden clustering, UMAP, rank_genes_groups) operate directly on CITEgeist outputs. The integration with Harmony enables joint visualization of programs across samples, as demonstrated in Module 5.

This interoperability is essential for practical adoption. Rather than requiring users to learn new analysis paradigms, CITEgeist extends existing workflows with spatially-aware preprocessing and cell type stratification. Researchers can apply familiar methods to the enhanced data structures CITEgeist provides.

## Figure 6 Legend (for reference)

(A) Workflow diagram showing CITEgeist outputs feeding into standard bioinformatics tools. Cell type proportions and deconvolved GEX layers (Module 3) enable PyDESeq2 differential expression within cell populations. Spatial programs (Module 4) enable GSEApy/Enrichr pathway enrichment. Program relationships (Module 5) enable network analysis with COMMOT and CellPhoneDB. All outputs in standard AnnData/CSV/Parquet formats compatible with scanpy workflows.

(B) Example application demonstrating the discovery-to-validation workflow. CITEgeist analysis (Modules 1-5) produces spatial programs that feed into downstream bioinformatics tools (PyDESeq2, GSEA, COMMOT), generating biological hypotheses for experimental validation. Demonstrated capabilities include: 3 responder-enriched and 4 progressor-enriched aligned programs (Module 5), 127 significant differentially expressed genes between responders and progressors, ECM organization and EMT pathway enrichment in progressors, and spatially-resolved ligand-receptor interactions via COMMOT.

(C) Output format compatibility table. CITEgeist exports all results in standard formats: cell type proportions (CSV, AnnData), deconvolved GEX layers (h5ad, parquet), spatial programs (CSV, AnnData with gene loadings and H matrices), bivariate relationships (CSV, JSON), and spatial coordinates (AnnData.obsm). Compatible with scanpy, seaborn, PyDESeq2, GSEA tools, networkx, Cytoscape, squidpy, and COMMOT.
