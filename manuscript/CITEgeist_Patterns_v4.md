# CITEgeist: A Spatially-Native Framework for Multi-Modal Integration and Program Discovery in Spatial Transcriptomics

Alexander Chih-Chieh Chang^1,2^, Brent T. Schlegel^1^, Neil Carleton^1^, Priscilla F. McAulife^1^, Steffi Oesterreich^1,2^, Russell Schwartz^2,3^, Adrian V. Lee^1,2,*^

^1^ Women's Cancer Research Center, UPMC Hillman Cancer Center and Magee-Women's Research Institute, Pittsburgh, PA
^2^ Department of Computational and Systems Biology, School of Medicine, University of Pittsburgh, Pittsburgh, PA
^3^ Computational Biology Department, School of Computer Science, Carnegie Mellon University, Pittsburgh, PA

\* Corresponding author: leeav@upmc.edu

---

## Abstract

Spatial transcriptomics technologies have transformed our ability to study tissue architecture, yet current computational methods often treat spatial information as an afterthought—applying single-cell analysis frameworks with coordinates appended as metadata. Here we present CITEgeist (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends), a framework designed from the ground up to be spatially-native, integration-first, and resolution-agnostic. By leveraging same-slide proteomics to anchor all downstream discovery, CITEgeist eliminates dependency on external single-cell atlases that can introduce artifacts from non-spatial reference data. We demonstrate CITEgeist across five integrated modules: automated marker detection and profile discovery using spatial colocalization; protein-anchored cell type deconvolution with neighborhood-aware refinement; spatially coherent gene expression program discovery; and cross-sample integration that identifies conserved relationships across patient cohorts. Applied to a 14-patient breast cancer cohort with matched CITE-seq spatial transcriptomics, CITEgeist identified response-associated transcriptional programs distinguishing immunotherapy responders from progressors. Benchmarking against established methods (Cell2Location, RCTD, Tangram, Seurat) on Xenium data with single-cell ground truth demonstrates competitive accuracy. CITEgeist's interpretable outputs integrate seamlessly with standard analysis tools including PyDESeq2 for differential expression, GSEApy for pathway analysis, and COMMOT for cell-cell communication inference. As spatial multi-omics technologies mature toward subcellular resolution with expanded protein panels, CITEgeist provides a principled framework for extracting biological insight from integrated spatial data.

---

## 1. Introduction

Spatial transcriptomics has emerged as a transformative technology for understanding tissue organization and cellular heterogeneity in situ. Technologies ranging from spot-based platforms like 10x Genomics Visium to single-cell resolution methods such as Xenium and CosMx now enable measurement of thousands of transcripts while preserving spatial context. Yet despite rapid technological advancement, computational methods have largely adapted existing single-cell analysis frameworks rather than developing approaches that treat spatial information as a first-class citizen.

The dominant paradigm for analyzing spatial transcriptomics data relies on transferring cell type annotations from external single-cell RNA sequencing (scRNA-seq) atlases. Methods such as Cell2Location, RCTD, and Tangram learn cell type signatures from reference datasets and project these onto spatial measurements. While powerful, this approach has fundamental limitations. Reference atlases are derived from dissociated tissues, losing the spatial context that may influence gene expression programs. Atlas composition may not match the tissue being studied—a reference built from healthy tissue may poorly represent tumor heterogeneity, and batch effects between technologies introduce systematic biases. Most critically, these methods cannot discover spatial patterns that don't exist in the reference, effectively "hallucinating" spatial organization from non-spatial data.

An alternative approach leverages the protein layer available from same-slide measurements. CITE-seq technology, originally developed for single-cell analysis, has been adapted for spatial transcriptomics platforms, enabling simultaneous measurement of transcripts and surface proteins from identical tissue sections. Unlike external references, same-slide proteomics captures the actual cellular composition of the tissue being studied. Protein markers with well-characterized cell type associations (CD3 for T cells, CD68 for macrophages, EPCAM for epithelial cells) provide direct anchors for downstream analysis without requiring reference transfer.

We present CITEgeist (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends), a computational framework built on three core principles:

**Spatially-native analysis.** Every module in CITEgeist incorporates spatial information as an integral component rather than an optional feature. Marker detection uses Moran's I spatial autocorrelation to identify proteins with organized expression patterns. Profile discovery leverages spatial colocalization—the tendency of certain markers to appear together in the same spots and adjacent neighborhoods. Proportion estimation applies Laplacian smoothing and neighborhood-aware refinement. Program discovery validates spatial coherence through bivariate spatial statistics. This contrasts with methods that merely use coordinates as additional features in otherwise non-spatial algorithms.

**Integration-first design.** CITEgeist was architected from inception for multi-modal spatial data combining proteomics and transcriptomics. Rather than analyzing modalities separately and attempting post-hoc integration, the framework uses protein measurements to inform transcriptomic analysis at every stage. Cell type profiles are derived from protein colocalization patterns. Proportions are estimated using protein markers as anchors. Gene expression is deconvolved using protein-derived cell type assignments. This tight integration ensures consistent interpretation across modalities.

**Resolution-agnostic architecture.** The same analytical framework operates on spot-level data (Visium, ~55 μm diameter, 1-10 cells per spot) and single-cell resolution data (Xenium, CosMx). Module implementations automatically adapt to the data structure—spatial neighbor definitions scale with technology resolution, and program discovery handles both aggregate and individual cell measurements. This flexibility positions CITEgeist for emerging technologies like Visium HD that bridge spot and single-cell scales.

CITEgeist comprises five integrated modules. Module 1 identifies spatially-variable protein markers worth analyzing. Module 2 discovers cell type profiles through spatial colocalization network analysis. Module 3 estimates cell type proportions and deconvolves gene expression using protein-derived profiles. Module 4 discovers gene expression programs within cell type-specific expression layers, validating spatial coherence through Moran's I. Module 5 integrates programs across multiple samples, identifying conserved relationships and response-associated signatures.

We demonstrate CITEgeist on a cohort of 14 breast cancer samples with matched CITE-seq spatial transcriptomics, identifying transcriptional programs that distinguish immunotherapy responders from progressors. Benchmarking against established deconvolution methods on Xenium data with single-cell ground truth shows competitive accuracy. We further demonstrate interoperability with standard bioinformatics tools, enabling PyDESeq2 differential expression analysis on CITEgeist-derived cell type stratifications and COMMOT-based cell-cell communication inference. All outputs are provided in standard scanpy/AnnData formats, facilitating integration into existing analysis workflows.

---

## 2. Results

### 2.1 Framework Overview

CITEgeist is organized as a modular pipeline where each stage leverages spatial information and passes interpretable outputs to subsequent modules (Figure 1A). The framework accepts spatial transcriptomics data with matched protein measurements in standard AnnData format, supporting both spot-level (Visium) and single-cell (Xenium, CosMx) resolutions.

**Module 1: Marker Interest Detection** evaluates each protein marker for spatial organization using three complementary statistics. Moran's I spatial autocorrelation identifies markers with non-random spatial patterns. Gaussian mixture modeling separates signal from background, computing signal-to-noise ratios. Kurtosis analysis detects peaked distributions indicative of localized expression. Markers passing these gates are forwarded to profile discovery.

**Module 2: Profile Assembly** discovers cell type marker profiles through spatial colocalization analysis. For each marker pair, the module computes same-spot co-occurrence, expression correlation, adjacent-spot enrichment, and bivariate Moran's I. Significance-filtered relationships are assembled into a colocalization network, where connected components represent distinct lineages. Hierarchical clustering with dynamic tree cutting extracts marker profiles, which are validated against known cell type associations.

**Module 3: Spatial Deconvolution** performs two-pass optimization. Pass 1 estimates cell type proportions using protein markers as anchors, applying Laplacian smoothing to enforce spatial coherence across neighborhoods. Neighborhood-aware finetuning adjusts estimates based on local context. Pass 2 deconvolves bulk gene expression into cell type-specific layers using the estimated proportions as constraints, producing per-spot, per-cell-type expression matrices.

**Module 4: Program Discovery** applies non-negative matrix factorization (NMF) to each cell type-specific expression layer, discovering gene expression programs within each population. Moran's I validates that discovered programs exhibit spatial coherence—programs with high spatial autocorrelation represent biologically meaningful spatial organization rather than technical noise. Bivariate Moran's I between program pairs identifies co-localized and mutually exclusive relationships.

**Module 5: Cross-Sample Integration** aligns programs across multiple samples using Harmony-based batch correction in a shared latent space. Hierarchical clustering identifies conserved programs appearing across patients. Response analysis compares program enrichment between clinical groups (e.g., responders vs progressors), enabling discovery of clinically relevant signatures.

Throughout the pipeline, outputs are stored in standard formats compatible with scanpy, enabling seamless integration with established analysis tools (Figure 1B-C).

### 2.2 Modules 1-2: Automated Profile Discovery from Spatial Colocalization

A key advantage of same-slide proteomics is the ability to discover cell type profiles directly from the spatial data, rather than relying on predefined marker lists. CITEgeist's Modules 1-2 implement an automated discovery pipeline that identifies which protein markers are informative and groups them into biologically coherent cell type profiles (Figure 2A-B).

Module 1 applies three statistical gates to filter the ~30-50 antibodies in a typical CITE-seq panel. Moran's I spatial autocorrelation (I > 0.1, p < 0.05) identifies markers with organized spatial patterns—randomly distributed markers provide little information about tissue architecture. Gaussian mixture modeling fits a two-component distribution to each marker, and markers with signal-to-noise ratio below 1.5 are excluded. Kurtosis analysis (κ > 2.0) identifies markers with peaked distributions, which often indicate cell type-specific expression. Markers passing either the Moran's I gate or the kurtosis gate (combined with the GMM filter) are considered "interesting" and forwarded to Module 2.

Module 2 constructs a spatial colocalization network from pairwise marker relationships. For each marker pair, we compute: (1) same-spot co-occurrence via Jaccard similarity on binarized expression; (2) Pearson correlation on continuous values; (3) adjacent-spot enrichment, measuring whether one marker's expression predicts the other in neighboring spots; and (4) bivariate Moran's I, quantifying spatial co-patterning. Edges passing significance thresholds (permutation-based p < 0.05) are retained, weighted by combined evidence.

Connected components in the filtered graph represent distinct cellular lineages—epithelial markers cluster separately from immune markers, for example. Within each component, hierarchical clustering with dynamic tree cutting extracts individual profiles. The optimal number of profiles is determined by reconstruction accuracy: we select the partition that best explains the observed protein measurements while maintaining non-redundancy.

To validate the approach, we applied Modules 1-2 to Xenium single-cell data where ground truth cell types were available from RNA-based clustering (Figure 2C). The discovered protein profiles correctly grouped known markers: CD3E and CD8A clustered together (T cells), CD68 and CD163 clustered together (macrophages), and EPCAM and KRT markers clustered together (epithelial cells). This demonstrates that spatial colocalization provides sufficient signal to recover biologically meaningful cell type definitions without external reference data (Figure 2D).

### 2.3 Module 3: Protein-Anchored Deconvolution with Spatial Regularization

Given cell type profiles from Module 2, Module 3 estimates cell type proportions at each spatial location and deconvolves bulk gene expression into cell type-specific layers (Figure 3A). Unlike reference-based methods that transfer signatures learned from external atlases, CITEgeist uses the protein measurements from the same tissue section as direct anchors.

Pass 1 formulates proportion estimation as a quadratic programming problem. For each spot, we find non-negative proportions that minimize reconstruction error between observed protein measurements and the profile-weighted mixture. To enforce spatial coherence, we add a Laplacian regularization term that penalizes differences between neighboring spots. This prevents noisy, salt-and-pepper proportion estimates while allowing genuine boundaries between tissue regions.

Following global optimization, neighborhood-aware finetuning adjusts proportions based on local context. For each spot, we examine the 6-8 nearest neighbors and their assignments, refining estimates where local evidence suggests corrections. This two-stage approach balances global consistency with local precision.

Pass 2 uses the estimated proportions to deconvolve gene expression. At each spot, bulk expression is modeled as a weighted sum of cell type-specific expression, where weights are the proportions from Pass 1. We solve for per-cell-type expression values using non-negative least squares with global and local enrichment priors derived from proportion-weighted spatial patterns. The result is a set of cell type-specific expression layers: for each gene at each spot, we obtain estimated expression within each cell type population.

We benchmarked CITEgeist against established deconvolution methods using Xenium data aggregated into pseudo-Visium spots where single-cell ground truth was available (Figure 3B). Across 5 tissue regions totaling 7,054 spots, CITEgeist achieved Pearson correlation r = 0.60 with ground truth proportions, comparable to reference-based methods Cell2Location (r = 0.61) and RCTD (r = 0.62). Jensen-Shannon divergence (lower is better) was 0.355 for CITEgeist versus 0.335 for Cell2Location and 0.347 for RCTD. Methods requiring label transfer performed substantially worse: Tangram achieved r = 0.14 and Seurat r = 0.17.

Importantly, CITEgeist achieved this accuracy without requiring any external reference data—the protein measurements from the same tissue section provided sufficient information for competitive deconvolution. This demonstrates that same-slide proteomics can replace reference atlases for proportion estimation while avoiding potential artifacts from reference-sample mismatches.

### 2.4 Module 4: Discovery of Spatially Coherent Gene Expression Programs

A central goal of spatial transcriptomics is discovering gene expression programs that exhibit organized spatial patterns. Module 4 applies non-negative matrix factorization (NMF) to the cell type-specific expression layers from Module 3, discovering programs within each cell type population (Figure 4A). Crucially, we validate that discovered programs are spatially coherent using Moran's I, ensuring they represent genuine spatial biology rather than technical variation.

For each cell type, we extract the deconvolved expression matrix (spots × genes) weighted by that cell type's proportion. NMF factorizes this matrix into program loadings (W matrix, genes × programs) and program activities (H matrix, spots × programs). We discover K programs per cell type, where K is set based on data complexity (typically 5-10).

To distinguish biologically meaningful programs from noise, we compute Moran's I spatial autocorrelation on each program's activity pattern across tissue space. Programs with significant positive Moran's I (I > 0.2, p < 0.01) exhibit spatial clustering—cells running these programs tend to co-localize. Programs with near-zero Moran's I have random spatial distributions and are less likely to represent coordinated tissue-level processes.

Applied to 5 Xenium tissue regions with 7 cell types, Module 4 discovered 175 total programs (25 programs per region, 5 per cell type) (Figure 4B-C). Moran's I values ranged from 0.04 to 0.37, with 68% of programs exceeding the I > 0.15 threshold for moderate spatial coherence. The highest spatial coherence was observed in fibroblast programs (mean I = 0.28), consistent with the known spatial organization of stromal cells in tissue architecture. T cell programs showed more variable coherence (I range 0.08-0.26), reflecting heterogeneity in immune infiltration patterns.

Example high-coherence programs included a fibroblast program enriched for HLA-DRA, FCGR3A, and HLA-DRB1 (Moran's I = 0.37), suggesting antigen-presenting fibroblast subpopulations with organized spatial distribution. A CD4+ T cell program with VIM, CD3E, and PTEN (I = 0.26) identified spatially clustered activated T cells (Figure 4D).

Module 4 also computes bivariate Moran's I between all program pairs, identifying co-localized programs (positive I) that may represent coordinated multicellular activities and mutually exclusive programs (negative I) that may represent competing cellular states (Figure 4E). This relationship network provides a spatial interaction map complementing conventional pathway analysis.

### 2.5 Module 5: Cross-Sample Integration Identifies Response-Associated Programs

To discover clinically relevant biology, Module 5 integrates programs across multiple patient samples, identifying conserved signatures and response-associated differences (Figure 5A). We applied this to a cohort of 14 breast cancer samples from 6 patients undergoing neoadjuvant therapy, including 5 samples from treatment responders and 9 from progressors.

Integration proceeds in three stages. First, program signatures (top genes per program) from all samples are embedded in a shared latent space using principal component analysis. Second, Harmony batch correction aligns programs across samples, removing patient-specific technical variation while preserving biological signal. Third, hierarchical clustering identifies groups of similar programs, which we term "aligned programs"—programs that appear consistently across patients.

From 590 individual programs across 14 samples, Module 5 identified 73 aligned programs representing conserved transcriptional states (Figure 5B). Of these, 5 programs appeared in more than 50% of samples, indicating highly conserved biology, while 36 were sample-specific. The aligned programs spanned all 7 cell types, with representation proportional to cell type abundance.

Response analysis compared program enrichment between responder and progressor samples (Figure 5C). Three programs were significantly enriched in responders: a macrophage program with FABP4/HBA2/TNXB, a macrophage program with VIM/TMSB4X/S100A6, and a cancer luminal program with MGP/MT-CO3/FOS. Four programs were enriched in progressors: a dendritic cell program with FTL/FN1/TIMP1, a cancer basal program with mitochondrial genes, a cancer luminal program with MGP/HBA2/LTF, and an endothelial program with KLHL17/MT-ND1/GPR101.

To validate these associations and identify specific genes driving response differences, we performed differential expression analysis using PyDESeq2 on pseudo-bulk aggregated expression (Figure 5D). Comparing 5 responder samples to 9 progressor samples across 13,371 tested genes, we identified 127 significantly differentially expressed genes (adjusted p < 0.05). The vast majority (122 genes) were upregulated in progressors, including matrix metalloproteinases (MMP13, MMP3, ADAMTS4, ADAMTS15), the necroptosis effector MLKL, and immune modulators (ALOX5AP, CLEC5A). Only 5 genes were upregulated in responders, including TMEM38B and TRIM72.

Module 5 also identifies conserved spatial relationships—program pairs that maintain consistent co-localization or exclusion across samples (Figure 5E). Of 191 conserved relationships detected, 26 (13.6%) showed consistent co-localization, 6 (3.1%) showed consistent mutual exclusion, and the remainder were independent. This relationship network highlights multicellular organization patterns that persist across patients despite inter-individual variation.

### 2.6 Interpretable Outputs Enable Downstream Analysis

A key design principle of CITEgeist is generating outputs compatible with established bioinformatics tools. All intermediate and final results are stored in standard formats (AnnData, CSV, JSON), enabling seamless integration with the broader computational biology ecosystem (Figure 6A).

**Differential expression.** Cell type-specific expression layers from Module 3 can be directly analyzed using PyDESeq2 or other DE frameworks. By filtering spots based on cell type proportions, researchers can perform DE analysis within specific populations—for example, comparing macrophage-enriched spots between conditions. The Module 5 response analysis demonstrated this workflow, identifying 127 DE genes between responder and progressor samples using standard PyDESeq2.

**Pathway analysis.** Gene lists from DE analysis or program loadings from Module 4 integrate directly with GSEApy/Enrichr for pathway enrichment. In our cohort, progressor-enriched genes showed enrichment for extracellular matrix organization (GO:0030198), collagen catabolic process (GO:0030574), and the MSigDB Hallmark epithelial-mesenchymal transition signature.

**Cell-cell communication.** CITEgeist's spatial neighbor graphs and cell type proportions are compatible with tools like COMMOT and CellPhoneDB for inferring ligand-receptor interactions. Previous analysis of this cohort identified midkine (MDK) as a potential mediator of estrogen receptor (ER) signaling, with subsequent wet lab validation confirming the mechanism (Supplementary Note 1) (Figure 6B-C).

**Clustering and visualization.** Standard scanpy workflows (Leiden clustering, UMAP, rank_genes_groups) operate directly on CITEgeist outputs. The integration with Harmony enables joint visualization of programs across samples, as demonstrated in Module 5.

This interoperability is essential for practical adoption. Rather than requiring users to learn new analysis paradigms, CITEgeist extends existing workflows with spatially-aware preprocessing and cell type stratification. Researchers can apply familiar methods to the enhanced data structures CITEgeist provides.

---

## 3. Discussion

CITEgeist addresses a fundamental gap in spatial transcriptomics analysis: the mismatch between data that is inherently spatial and methods designed for non-spatial contexts. While established tools like Cell2Location, RCTD, and Tangram have advanced the field significantly, they share a common architecture—learning cell type signatures from external scRNA-seq references and projecting them onto spatial data. This reference-dependent paradigm has three limitations that CITEgeist's spatial-native, integration-first design overcomes.

**From spatial-aware to spatially-native.** Many existing methods incorporate spatial information as regularization or post-processing—coordinates constrain an otherwise non-spatial algorithm. CITEgeist inverts this relationship: spatial statistics drive discovery at every stage. Moran's I filters markers in Module 1, spatial colocalization defines profiles in Module 2, Laplacian smoothing regularizes proportions in Module 3, and spatial coherence validates programs in Module 4. This difference is more than philosophical. Reference-based methods can only discover patterns present in the reference; if a spatial cell state doesn't exist in dissociated single-cell data, it cannot be transferred. CITEgeist discovers patterns from the spatial data itself, using the protein layer as ground truth rather than an external atlas.

**Same-slide anchoring eliminates reference artifacts.** A practical advantage of same-slide proteomics is the elimination of batch effects between reference and query data. scRNA-seq atlases are collected from different individuals, processed with different protocols, and sequenced on different platforms than the spatial data being analyzed. These technical differences propagate through reference-based pipelines, potentially introducing systematic biases. By using protein measurements from the identical tissue section, CITEgeist avoids this source of error. Our benchmarking on Xenium data demonstrates that same-slide proteomics provides sufficient information for competitive deconvolution accuracy (Pearson r = 0.60 vs 0.61 for Cell2Location) without any reference data.

**Architectural readiness for emerging technologies.** The spatial multi-omics field is evolving rapidly. Visium HD achieves 2 μm resolution approaching single-cell, while expanded antibody panels (100+ targets) and spatial proteomics modalities (CODEX, IMC) provide increasingly rich protein information. CITEgeist's resolution-agnostic design—demonstrated on both spot-level Visium and single-cell Xenium—positions it for these developments. As protein panels expand, Module 2's colocalization-based profile discovery will leverage additional markers automatically. As resolution increases, the same modules apply with adjusted spatial neighborhood definitions. This contrasts with methods requiring substantial rearchitecting for new data types.

**Limitations and future directions.** CITEgeist's reliance on same-slide proteomics is both a strength and a limitation. Current spatial platforms support 30-50 antibodies, constraining the granularity of discoverable cell types. Rare populations without specific protein markers may be underrepresented. As antibody panels expand, this limitation will diminish, but near-term applications may benefit from hybrid approaches combining same-slide proteomics with targeted reference information for specific rare populations.

The current implementation assumes that protein markers have relatively consistent meaning across cell types—a marker positive for macrophages should indicate macrophage presence regardless of tissue context. While this holds for canonical markers (CD68, CD3, EPCAM), some markers show context-dependent expression. Future versions could incorporate tissue-specific priors or learn marker-context relationships from training data.

Computationally, Module 3's quadratic programming optimization scales well to typical Visium datasets (~5,000 spots) but may require approximations for very large datasets (Xenium with millions of cells). Current implementation uses checkpointing and can process such datasets in batches, but native scalability improvements would benefit large-scale applications.

**Conclusions.** CITEgeist provides a principled framework for spatial multi-omics analysis that treats spatial information as foundational rather than supplementary. By anchoring all discovery on same-slide protein measurements, it eliminates reference dependency while maintaining competitive accuracy. The five-module architecture progresses from raw data through cell type identification, proportion estimation, program discovery, and cross-sample integration—each stage producing interpretable outputs compatible with established tools. As spatial technologies continue advancing toward higher resolution and richer multi-modal measurements, frameworks designed from inception for integrated spatial analysis will become increasingly essential.

---

## 4. Methods

### 4.1 Data Availability

Spatial transcriptomics data from the 14-patient breast cancer cohort will be deposited in GEO upon publication. Xenium benchmarking data is available from 10x Genomics (Xenium Human Breast Cancer Dataset). Code is available at https://github.com/leeoesterreich/CITEgeist.

### 4.2 Module 1: Marker Interest Detection

For each protein marker, we compute three statistics. Moran's I spatial autocorrelation is calculated using the formula:

$$I = \frac{n}{\sum_i \sum_j w_{ij}} \frac{\sum_i \sum_j w_{ij}(x_i - \bar{x})(x_j - \bar{x})}{\sum_i (x_i - \bar{x})^2}$$

where $w_{ij}$ is the spatial weight between spots $i$ and $j$ (1 if neighbors, 0 otherwise), and $n$ is the number of spots. Significance is assessed via permutation testing (1000 permutations).

Gaussian mixture modeling fits a 2-component GMM to each marker's expression distribution. Signal-to-noise ratio is computed as the ratio of the higher component's mean to the lower component's mean.

Kurtosis is computed as the fourth standardized moment:

$$\kappa = \frac{E[(X - \mu)^4]}{\sigma^4}$$

Markers are considered "interesting" if they pass either: (1) Moran's I > 0.1 with p < 0.05, or (2) kurtosis > 2.0; combined with GMM SNR > 1.5.

### 4.3 Module 2: Profile Discovery

Pairwise marker relationships are computed for all marker pairs passing Module 1. Same-spot co-occurrence uses Jaccard similarity on binarized expression (threshold: 75th percentile). Adjacent-spot enrichment computes the average expression of marker B in neighbors of spots expressing marker A, compared to random expectation. Bivariate Moran's I uses the cross-product form:

$$I_{XY} = \frac{n}{\sum_i \sum_j w_{ij}} \frac{\sum_i \sum_j w_{ij}(x_i - \bar{x})(y_j - \bar{y})}{\sqrt{\sum_i (x_i - \bar{x})^2 \sum_j (y_j - \bar{y})^2}}$$

Edges with permutation p < 0.05 for at least two statistics are retained. Hierarchical clustering uses Ward's method with Euclidean distance on the edge weight matrix. Profile count is selected by maximizing reconstruction accuracy (R² between predicted and observed protein expression) while maintaining minimum profile size of 2 markers.

### 4.4 Module 3: Spatial Deconvolution

Pass 1 solves the quadratic program:

$$\min_{\pi} \|P\pi - y\|_2^2 + \lambda \pi^T L \pi$$

where $P$ is the profile matrix (markers × cell types), $\pi$ are proportions (constrained $\pi_i \geq 0$, $\sum_i \pi_i = 1$), $y$ is observed protein expression, and $L$ is the graph Laplacian for spatial regularization. Default $\lambda = 0.1$.

Neighborhood finetuning adjusts each spot's proportions based on consistency with neighbor estimates, using weighted averaging with distance decay.

Pass 2 solves for gene expression:

$$\min_{G} \|G \pi - e\|_2^2 + \lambda_g \|G - G_{prior}\|_2^2$$

where $G$ is cell type-specific expression (genes × cell types), $e$ is observed bulk expression, and $G_{prior}$ is computed from global expression patterns weighted by proportions.

### 4.5 Module 4: Program Discovery

NMF factorizes the weighted expression matrix $X = WH$ where $W$ (genes × programs) contains program loadings and $H$ (programs × spots) contains program activities. We use sklearn's NMF with Frobenius norm loss and coordinate descent solver.

For each program (column of $H$), we compute Moran's I on the activity values across spots. Programs with I > 0.15 and p < 0.01 are considered spatially coherent.

Bivariate Moran's I is computed between all program pairs to identify spatial relationships.

### 4.6 Module 5: Cross-Sample Integration

Programs from all samples are embedded using PCA (n_components=30) on concatenated gene loadings. Harmony batch correction removes sample-specific variation using default parameters (theta=2.0, max_iter=20).

Aligned programs are identified via hierarchical clustering (Ward linkage) on the corrected embedding, with cluster count selected by silhouette score optimization.

Response analysis computes mean program activity in responder vs progressor samples, with significance assessed by Mann-Whitney U test (Benjamini-Hochberg correction).

### 4.7 Differential Expression Analysis

Pseudo-bulk aggregation sums gene expression within each sample. PyDESeq2 (v0.4) performs differential expression with design formula ~condition (responder vs progressor). Genes with adjusted p < 0.05 are considered significant.

### 4.8 Benchmarking

Xenium single-cell data was aggregated into pseudo-Visium spots (55 μm diameter) to create ground truth. Cell type proportions were computed from single-cell annotations within each spot. Methods were evaluated on Pearson correlation, RMSE, MAE, and Jensen-Shannon divergence between predicted and true proportions.

---

## Acknowledgments

This work was supported by the Breast Cancer Research Foundation, the Nicole Meloche Foundation, the Magee-Women's Research Institute, Susan G. Komen, and NIH grants R01CA285697, P30CA047904, and T32CA082084. We thank the patients who contributed tissue samples for this research.

## Author Contributions

A.C.C. developed the CITEgeist framework, implemented all modules, and performed analyses. B.T.S. contributed to algorithm development and benchmarking. N.C. and P.F.M. provided clinical samples and interpretation. S.O., R.S., and A.V.L. supervised the project. A.C.C. and A.V.L. wrote the manuscript with input from all authors.

## Declaration of Interests

The authors declare no competing interests.

---

## Figure Legends

**Figure 1. CITEgeist Framework Overview.**
(A) Schematic of the five-module pipeline: Module 1 (Marker Interest Detection), Module 2 (Profile Assembly), Module 3 (Spatial Deconvolution), Module 4 (Program Discovery), Module 5 (Cross-Sample Integration). Arrows indicate data flow between modules.
(B) Spatial-native operations highlighted at each stage: Moran's I for marker detection, colocalization networks for profile discovery, Laplacian regularization for proportion estimation, spatial coherence validation for programs.
(C) Resolution flexibility: the same framework operates on spot-level (Visium) and single-cell (Xenium) data.

**Figure 2. Modules 1-2: Automated Profile Discovery.**
(A) Marker interest detection using three statistical gates: Moran's I spatial autocorrelation, Gaussian mixture modeling for signal-to-noise, and kurtosis analysis.
(B) Profile assembly from spatial colocalization: pairwise marker relationships → network construction → hierarchical clustering → profile extraction.
(C) Xenium single-cell demonstration showing discovered profiles from spatial colocalization analysis.
(D) Validation: discovered profiles correctly recover known cell type markers (CD3E/CD8A for T cells, CD68/CD163 for macrophages, EPCAM/KRT for epithelial).

**Figure 3. Module 3: Spatial Deconvolution Benchmarking.**
(A) Two-pass deconvolution schematic: Pass 1 estimates proportions from protein markers with Laplacian regularization; Pass 2 deconvolves gene expression using proportion constraints.
(B) Benchmarking on Xenium pseudo-Visium spots (n=7,054): comparison of CITEgeist, Cell2Location, RCTD, Tangram, and Seurat on Pearson correlation, Jensen-Shannon divergence, and RMSE.
(C) Example spatial visualization of estimated proportions across a tissue region.

**Figure 4. Module 4: Spatial Program Discovery.**
(A) NMF-based program discovery schematic: cell type-specific expression layers → matrix factorization → program loadings and activities → Moran's I validation.
(B) Program examples: top genes for selected programs across cell types with associated Moran's I values.
(C) Moran's I validation: scatter plot of program activity variance vs spatial coherence (I values), with threshold for biologically meaningful programs.
(D) Xenium single-cell programs demonstrating resolution-agnostic application.
(E) Bivariate relationships: heatmap of program-program spatial correlations identifying co-localized and mutually exclusive pairs.

**Figure 5. Module 5: Cross-Sample Integration.**
(A) Integration schematic: per-sample programs → Harmony alignment → aligned programs → response analysis.
(B) UMAP visualization of 590 programs from 14 samples, colored by cell type, showing successful integration.
(C) Responder vs progressor enriched programs: bar plot showing 3 responder-enriched and 4 progressor-enriched aligned programs.
(D) PyDESeq2 volcano plot: 127 significant genes (5 responder-up, 122 progressor-up) from pseudo-bulk differential expression analysis.
(E) Conserved spatial relationship network: edges represent program pairs with consistent co-localization (green) or exclusion (red) across samples.

**Figure 6. Interpretable Outputs and Interoperability.**
(A) Workflow diagram showing CITEgeist outputs feeding into standard tools: PyDESeq2 for differential expression, GSEApy for pathway analysis, COMMOT for cell-cell communication.
(B) Midkine (MDK) discovery summary: spatially resolved ligand-receptor signaling identified through COMMOT analysis of CITEgeist outputs.
(C) Wet lab validation: summary of experimental confirmation of computationally predicted mechanisms.

