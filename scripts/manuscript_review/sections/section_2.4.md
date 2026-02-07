# Section 2.4: Module 4 - Discovery of Spatially Coherent Gene Expression Programs

**Referenced Figures**: Figure 4 (A, B, C, D, E)

## Section Text

A central goal of spatial transcriptomics is discovering gene expression programs that exhibit organized spatial patterns. Module 4 applies non-negative matrix factorization (NMF) to the cell type-specific expression layers from Module 3, discovering programs within each cell type population (Figure 4A). Crucially, we validate that discovered programs are spatially coherent using Moran's I, ensuring they represent genuine spatial biology rather than technical variation.

For each cell type, we extract the deconvolved expression matrix (spots × genes) weighted by that cell type's proportion. NMF factorizes this matrix into program loadings (W matrix, genes × programs) and program activities (H matrix, spots × programs). We discover K programs per cell type, where K is set based on data complexity (typically 5-10).

To distinguish biologically meaningful programs from noise, we compute Moran's I spatial autocorrelation on each program's activity pattern across tissue space. Programs with significant positive Moran's I (I > 0.2, p < 0.01) exhibit spatial clustering—cells running these programs tend to co-localize. Programs with near-zero Moran's I have random spatial distributions and are less likely to represent coordinated tissue-level processes.

Applied to 5 Xenium tissue regions with 7 cell types, Module 4 discovered 175 total programs (25 programs per region, 5 per cell type) (Figure 4B-C). Moran's I values ranged from 0.04 to 0.37, with 68% of programs exceeding the I > 0.15 threshold for moderate spatial coherence. The highest spatial coherence was observed in fibroblast programs (mean I = 0.28), consistent with the known spatial organization of stromal cells in tissue architecture. T cell programs showed more variable coherence (I range 0.08-0.26), reflecting heterogeneity in immune infiltration patterns.

Example high-coherence programs included a fibroblast program enriched for HLA-DRA, FCGR3A, and HLA-DRB1 (Moran's I = 0.37), suggesting antigen-presenting fibroblast subpopulations with organized spatial distribution. A CD4+ T cell program with VIM, CD3E, and PTEN (I = 0.26) identified spatially clustered activated T cells (Figure 4D).

Module 4 also computes bivariate Moran's I between all program pairs, identifying co-localized programs (positive I) that may represent coordinated multicellular activities and mutually exclusive programs (negative I) that may represent competing cellular states (Figure 4E). This relationship network provides a spatial interaction map complementing conventional pathway analysis.

## Figure 4 Legend (for reference)

(A) NMF-based program discovery schematic. For each cell type, the deconvolved expression matrix (spots × genes) weighted by cell type proportion is factorized using non-negative matrix factorization into program loadings (W matrix, genes × programs) and program activities (H matrix, spots × programs). Moran's I spatial autocorrelation on program activities validates that discovered programs exhibit spatial coherence.

(B) Program examples across cell types. Top genes for selected high-coherence programs with associated Moran's I values. Fibroblast programs show highest spatial coherence (mean I = 0.28), consistent with organized stromal architecture. Example: fibroblast program with HLA-DRA, FCGR3A, HLA-DRB1 (I = 0.37) suggests antigen-presenting fibroblast subpopulations.

(C) Moran's I validation. Scatter plot of program activity variance versus spatial coherence (Moran's I values). Programs above I > 0.15 threshold (dashed line) are considered spatially coherent. 68% of discovered programs exceed this threshold. Xenium single-cell data, 5 regions, 7 cell types, 175 total programs.

(D) Spatial visualization of program activities for selected high-coherence programs, demonstrating organized spatial patterns rather than random noise.

(E) Bivariate relationships heatmap. Bivariate Moran's I computed between all program pairs identifies co-localized programs (positive I, programs active in same spatial regions) and mutually exclusive programs (negative I, programs with anti-correlated spatial patterns). This relationship network provides a spatial interaction map complementing pathway analysis.
