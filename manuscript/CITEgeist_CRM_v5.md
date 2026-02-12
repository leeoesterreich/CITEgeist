# CITEgeist: Accurate deconvolution of spatial transcriptomics with same-slide proteomics reveals midkine as a secreted microenvironment modulator in ESR1 mutant breast cancer

Alexander Chih-Chieh Chang^1,2,\*^, Brent T. Schlegel^1,2,\*^, Neil Carleton^1,2^, Priscilla F. McAuliffe^1,3^, Hunter Waltermire^1,2^, Steffi Oesterreich^1,4^, Russell Schwartz^5,6^, and Adrian V. Lee^1,4,7,#^

^1^ Women's Cancer Research Center, UPMC Hillman Cancer Center, Magee-Womens Research Institute, Pittsburgh, PA, USA
^2^ University of Pittsburgh, School of Medicine, Pittsburgh, PA, USA
^3^ Department of Surgery, Division of Breast Surgical Oncology, University of Pittsburgh School of Medicine, Pittsburgh, PA, USA
^4^ Department of Pharmacology and Chemical Biology, University of Pittsburgh, Pittsburgh, PA, USA
^5^ Ray and Stephanie Lane Computational Biology Department, Carnegie Mellon University, Pittsburgh, PA, USA
^6^ Department of Biological Sciences, Carnegie Mellon University, Pittsburgh, PA, USA
^7^ Institute of Precision Medicine, Pittsburgh, PA, USA

\* These authors contributed equally to this manuscript.
\# Corresponding author: leeav@upmc.edu

---

## Summary

Spatial transcriptomics methods typically rely on external single-cell RNA-seq atlases for cell type deconvolution, introducing batch effects and limiting discovery to patterns present in dissociated reference data. Here we present CITEgeist, a framework that instead leverages same-slide proteomics to anchor deconvolution directly in the tissue being studied. CITEgeist automatically discovers cell type profiles from spatial protein colocalization, estimates cell type proportions via spatially-regularized optimization, deconvolves gene expression into cell type-specific layers, and integrates programs across patient samples. Applied to clinical specimens from a primary endocrine therapy trial for older women with ER+ breast cancer (NCT05914792), CITEgeist identified upregulated midkine signaling in an ESR1 D538G-mutant case, validated by ELISA and immunofluorescence. Cross-sample integration across 12 specimens revealed response-associated transcriptional programs distinguishing responders from progressors. Benchmarking on Xenium single-cell ground truth demonstrates accuracy comparable to reference-based methods without requiring external references.

---

## Motivation

Current spatial transcriptomics deconvolution methods depend on external single-cell RNA-seq reference atlases, which are costly to generate, may not match the tissue under study, and risk imposing non-spatial expression patterns onto spatial data. CITEgeist addresses this gap by using same-slide protein measurements as direct anchors for cell type identification and gene expression deconvolution, eliminating reference dependency entirely. This approach is particularly advantageous in heterogeneous cancer tissues, where atlas composition rarely reflects the complexity of post-treatment tumor microenvironments, and where spatial cell states may not exist in dissociated references.

---

## Highlights

- Same-slide proteomics replaces external atlases for spatial deconvolution
- Automated profile discovery from spatial protein colocalization patterns
- MDK signaling upregulation in ESR1-mutant breast cancer validated in vitro
- Cross-sample integration reveals endocrine therapy response signatures

---

## eTOC Blurb

Chang, Schlegel et al. present CITEgeist, a spatial transcriptomics framework that uses same-slide proteomics rather than external single-cell references for cell type deconvolution. Applied to breast cancer clinical specimens, CITEgeist identified midkine signaling upregulation in ESR1-mutant tumors, validated experimentally, and discovered endocrine therapy response-associated programs across patients.

---

## Introduction

Estrogen receptor-positive (ER+) breast cancer accounts for approximately 85% of new breast cancer diagnoses^1^, with surgery and anti-estrogen therapy established as first-line treatment^2^. Spatial transcriptomics offers a unique opportunity to dissect how these cancers evolve in response to therapy within their native tissue architecture. Technologies from spot-based platforms such as 10x Genomics Visium to single-cell resolution methods like Xenium and CosMx now enable measurement of thousands of transcripts while preserving spatial context^17^. Yet computational methods have largely adapted existing single-cell analysis frameworks rather than developing approaches that treat spatial information as foundational^18^.

The dominant paradigm for analyzing spatial transcriptomics data relies on transferring cell type annotations from external single-cell RNA sequencing (scRNA-seq) atlases. Methods such as Cell2Location^8^, RCTD, and Tangram^16^ learn cell type signatures from reference datasets and project these onto spatial measurements. While powerful, this approach has fundamental limitations^7^. Reference atlases derive from dissociated tissues, losing spatial context that may influence gene expression programs^5,6^. Atlas composition may not match the tissue being studied; a reference built from healthy tissue may poorly represent tumor heterogeneity^26^, and batch effects between technologies introduce systematic biases. Most critically, these methods cannot discover spatial patterns absent from the reference, effectively imposing non-spatial organization onto spatial data^7^.

An alternative leverages the protein layer available from same-slide measurements. CITE-seq technology, adapted for spatial transcriptomics platforms, enables simultaneous measurement of transcripts and surface proteins from identical tissue sections. Unlike external references, same-slide proteomics captures the actual cellular composition of the tissue being studied. Protein markers with well-characterized cell type associations (CD3 for T cells, CD68 for macrophages, EPCAM for epithelial cells) provide direct anchors for downstream analysis without requiring reference transfer. However, no existing tool fully integrates this proteomic information for spatial deconvolution.

We present CITEgeist (Cellular Indexing of Transcriptomes and Epitopes for Guided Exploration of Intrinsic Spatial Trends), a computational framework built on three principles. First, CITEgeist is spatially native: spatial statistics drive discovery at every stage, from Moran's I-based marker detection through Laplacian-regularized proportion estimation to spatially coherent program validation. This contrasts with methods that merely append coordinates to otherwise non-spatial algorithms. Second, CITEgeist takes an integration-first approach, using protein measurements to inform transcriptomic analysis at every module rather than analyzing modalities separately. Third, the same analytical framework operates across resolutions, from spot-level Visium to single-cell Xenium, with spatial neighborhood definitions that scale automatically.

CITEgeist comprises five integrated modules: automated marker detection and profile discovery (Modules 1-2), cell type proportion estimation and gene expression deconvolution (Module 3), spatially coherent program discovery (Module 4), and cross-sample integration (Module 5). We demonstrate CITEgeist on clinical specimens from a prospective trial evaluating primary endocrine therapy in older women with ER+ breast cancer^14,15^, where it identified upregulated midkine signaling in an ESR1-mutant case, validated through in vitro experiments. We further apply the full pipeline across 12 specimens from 6 patients to discover response-associated transcriptional programs, and benchmark against established methods on Xenium data with single-cell ground truth.

---

## Results

### 2.1 Framework Overview

CITEgeist is organized as a modular pipeline where each stage leverages spatial information and passes interpretable outputs to subsequent modules (Figure 1A). The framework accepts spatial transcriptomics data with matched protein measurements in standard AnnData format, supporting both spot-level (Visium) and single-cell (Xenium, CosMx) resolutions. Unlike existing deconvolution methods that require external scRNA-seq reference atlases, CITEgeist operates entirely from same-slide protein measurements, eliminating batch effects and reference selection bias.

Module 1 evaluates each protein marker for spatial organization using three complementary statistics: Moran's I spatial autocorrelation identifies markers with non-random spatial patterns, Gaussian mixture modeling separates signal from background to compute signal-to-noise ratios, and kurtosis analysis detects peaked distributions indicative of localized expression. Markers passing these gates are forwarded to Module 2, which discovers cell type marker profiles through spatial colocalization analysis. For each marker pair, the module computes same-spot co-occurrence, expression correlation, adjacent-spot enrichment, and bivariate Moran's I. Significance-filtered relationships are assembled into a colocalization network, where connected components represent distinct lineages and hierarchical clustering with dynamic tree cutting extracts individual profiles.

Module 3 performs two-pass optimization: Pass 1 estimates cell type proportions using protein markers as anchors with Laplacian smoothing to enforce spatial coherence, followed by neighborhood-aware finetuning. Pass 2 deconvolves bulk gene expression into cell type-specific layers using the estimated proportions as constraints. Module 4 applies non-negative matrix factorization to each cell type-specific expression layer, discovering gene expression programs whose spatial coherence is validated through Moran's I. Module 5 aligns programs across multiple samples using Harmony-based batch correction and identifies conserved transcriptional states and response-associated differences. Throughout the pipeline, outputs are stored in standard formats compatible with scanpy and interoperable with tools such as PyDESeq2, GSEApy, and COMMOT (Figure 1B-C).

### 2.2 Automated Profile Discovery from Spatial Colocalization

A key advantage of same-slide proteomics is the ability to discover cell type profiles directly from the spatial data rather than relying on predefined marker lists. Modules 1-2 implement an automated discovery pipeline that identifies informative protein markers and groups them into biologically coherent cell type profiles (Figure 2A-B).

Module 1 applies three statistical gates to filter the approximately 30-50 antibodies in a typical CITE-seq panel. Moran's I spatial autocorrelation (I > 0.1, p < 0.05) identifies markers with organized spatial patterns, as randomly distributed markers provide little information about tissue architecture. Gaussian mixture modeling fits a two-component distribution to each marker, excluding those with signal-to-noise ratio below 1.5. Kurtosis analysis (kappa > 2.0) identifies markers with peaked distributions often indicative of cell type-specific expression. Markers passing either the Moran's I gate or the kurtosis gate, combined with the GMM filter, are forwarded to Module 2.

Module 2 constructs a spatial colocalization network from pairwise marker relationships. For each marker pair, we compute same-spot co-occurrence via Jaccard similarity on binarized expression, Pearson correlation on continuous values, adjacent-spot enrichment measuring whether one marker's expression predicts the other in neighboring spots, and bivariate Moran's I quantifying spatial co-patterning. Edges passing significance thresholds (permutation-based p < 0.05) are retained, weighted by combined evidence. Connected components in the filtered graph represent distinct cellular lineages; epithelial markers cluster separately from immune markers, for example. Within each component, hierarchical clustering with dynamic tree cutting extracts individual profiles, and the optimal number of profiles is determined by reconstruction accuracy (R^2^ between predicted and observed protein expression) while maintaining non-redundancy.

To validate this approach, we applied Modules 1-2 to Xenium single-cell data where ground truth cell types were available from RNA-based clustering (Figure 2C). The discovered protein profiles correctly grouped known markers: CD3E and CD8A clustered together (T cells), CD68 and CD163 clustered together (macrophages), and EPCAM and KRT markers clustered together (epithelial cells), demonstrating that spatial colocalization provides sufficient signal to recover biologically meaningful cell type definitions without external reference data (Figure 2D).

### 2.3 Protein-Anchored Deconvolution Achieves Competitive Accuracy Without External References

To validate that Modules 1-2 recover biologically meaningful profiles, we applied them to Xenium renal cell carcinoma (RCC) pseudo-Visium data where ground truth cell types were available from RNA-based clustering. The automated discovery correctly identified 7 cell type profiles: B cells (CD20), CD4+ T cells (CD3E/CD4), CD8+ T cells (CD3E/CD8A), macrophages (CD68/CD163), endothelial cells (CD31), epithelial cells (PanCK), and fibroblasts (alphaSMA). Of these, 6 of 7 profiles matched known marker assignments exactly, with one profile recovered as a superset containing an additional marker (Figure 3A).

Given these profiles, Module 3 estimates cell type proportions at each spatial location and deconvolves bulk gene expression into cell type-specific layers. Pass 1 formulates proportion estimation as a quadratic programming problem: for each spot, we find non-negative proportions that minimize reconstruction error between observed protein measurements and the profile-weighted mixture. To enforce spatial coherence, a Laplacian regularization term penalizes differences between neighboring spots, preventing noisy salt-and-pepper estimates while allowing genuine boundaries between tissue regions. Following global optimization, neighborhood-aware finetuning adjusts proportions based on local context by examining the 6-8 nearest neighbors. Pass 2 uses the estimated proportions to deconvolve gene expression, modeling bulk expression at each spot as a weighted sum of cell type-specific expression and solving for per-cell-type expression values using non-negative least squares with global and local enrichment priors.

We rigorously evaluated CITEgeist against established deconvolution methods across simulated and real tissue data, focusing particularly on the heterogeneous cancer architectures where existing methods struggle most. On real tissue data, we aggregated Xenium breast cancer single-cell data into pseudo-Visium spots where ground truth cell type proportions could be computed from single-cell annotations. Across 5 tissue regions totaling 7,054 spots, CITEgeist achieved Pearson correlation r = 0.60 +/- 0.05, comparable to Cell2Location^8^ (r = 0.61 +/- 0.04) and RCTD (r = 0.62 +/- 0.03). CITEgeist achieved the lowest absolute error metrics among all methods (RMSE = 0.167 +/- 0.006, MAE = 0.115 +/- 0.005), compared to Cell2Location (RMSE = 0.179 +/- 0.017) and RCTD (RMSE = 0.177 +/- 0.004). Pairwise comparisons showed no significant differences between CITEgeist and the reference-based methods (Wilcoxon signed-rank test: vs Cell2Location p = 0.31, vs RCTD p = 0.19). Methods designed for label transfer performed substantially worse: Tangram^16^ (r = 0.14 +/- 0.08) and Seurat^10^ (r = 0.17 +/- 0.07), confirming task mismatch (Figure 3B). Per-spot predicted versus ground truth proportions showed overall concordance across all cell types (Figure 3E), and spatial comparison maps confirmed that CITEgeist faithfully recapitulated the spatial distributions of epithelial, macrophage, and T cell populations (Figure 3F).

Using simulated Visium-like datasets generated with scCube^12^ from an ER+ scRNA-seq atlas^11^, we tested performance across both highly segmented and mixed (cancer-like) tissue compositions. CITEgeist achieved Pearson correlation r = 0.95 and Jensen-Shannon divergence of 0.16, maintaining consistent performance across both tissue types (RMSE 0.08 in both conditions). In contrast, reference-based methods showed degraded performance on mixed tissues: Cell2Location's RMSE increased from 0.08 to 0.167, and Seurat's from 0.10 to 0.133 when moving from segmented to mixed architectures. To test sensitivity to reference quality, we compared performance using comprehensive (30,000-cell) versus realistic single-sample (8,000-cell) references. RCTD showed severe reference dependency, with RMSE increasing 4-fold from 0.05 to 0.21 with the smaller reference, while other methods were minimally affected. Tangram notably over-called T cells due to their known over-representation in dissociated single-cell datasets^17^, reporting average JSD of 0.56 in segmented data compared to 0.15 for other methods (Figure 3C).

For gene expression deconvolution, we compared CITEgeist against methods capable of producing per-cell-type expression estimates. On the simulated benchmarks, CITEgeist achieved RMSE of 0.293 in highly segmented tissue and 0.257 in mixed tissue, compared to Cell2Location (0.493 and 1.303, respectively) and Tangram (1.269 and 1.436), representing a 1.7- to 5.6-fold improvement depending on method and condition (Figure 3D). CITEgeist's advantage was most pronounced in mixed (cancer-like) tissue, where Cell2Location's gene expression RMSE increased 5.1-fold and Tangram's 5.6-fold relative to CITEgeist. CITEgeist's complete spatial coverage is critical for downstream spatial analyses and program discovery. Taken together, these results demonstrate that same-slide proteomics provides sufficient information for competitive deconvolution without requiring any external reference data, while offering robustness to the tissue heterogeneity and reference quality issues that challenge existing methods.

### 2.4 CITEgeist Identifies Upregulated Midkine Signaling in ESR1-Mutant Breast Cancer

After establishing CITEgeist's accuracy, we applied it to clinical samples from a prospective trial evaluating primary endocrine therapy (aromatase inhibitors) in older women with ER+ breast cancer who forego upfront surgery (NCT05914792)^14,15^. We focused on a case harboring an ESR1 D538G mutation, initially identified through pilot bulk RNA sequencing and confirmed by droplet digital PCR. This mutation was of particular interest given the clinical relevance of ESR1 mutations in endocrine resistance and our laboratory's prior research in this area^13^.

We have previously shown that ESR1-mutant breast cancers express basal cytokeratins^13^, and CITEgeist deconvolution on spatial transcriptomics from this case confirmed the presence of basal gene signatures in the cancer cell layer, visualized on the H&E tissue backdrop (Figure 4A). Notably, these signals were not clustered in a single location as one might expect from a single subclone but were instead dispersed across the tissue section (Figure 4B), despite exhibiting the expected gene mutation signatures from EstroGene^20^ (Figure 4C) and known pathway alterations from prior work^13^ (Figure 4D). As one of the first spatial analyses of ESR1-mutant breast cancer in the native breast microenvironment, we examined the signaling effects of these mutant cancer cells.

A key strength of CITEgeist is that its use of linear assignment preserves raw integer transcript counts, and all outputs are stored in standard AnnData format. This interoperability enabled us to directly apply COMMOT^27^, a signaling prediction tool designed for spatial data, to the deconvolved cell type layers. To validate the biological relevance of this approach, we first confirmed that CITEgeist-derived secretion signals correlated with established cell type-specific secretion profiles from the Human Protein Atlas (Spearman rho = 0.46, p = 2.86e-37). The COMMOT analysis then revealed increased MDK (midkine), PTN, and MIF signaling in the mutant regions (Figure 4F). Spatial mapping of the MDK-SDC4 ligand-receptor communication score confirmed that midkine signaling was spatially co-localized with D538G-high regions (Figure 4E). Midkine has recently been implicated as a driver of age-related changes and mammary tumorigenesis^21,22^, making this an immediately relevant finding in the context of our older patient population.

To validate these computational findings, we performed ELISA on conditioned media from MCF7 cell lines harboring the ESR1 D538G mutation. The mutant cells secreted approximately double the midkine in both estrogen-deprived and estradiol-supplemented conditions (p < 0.0001) (Figure 4G), consistent with the constitutive receptor activity conferred by D538G. Because midkine is known to have a strong pericellular effect, with a significant proportion retained near the cell surface in non-cancerous settings^23^, we also conducted immunofluorescence analysis on the same cell lines. In the estrogen-deprived environment mimicking a patient on endocrine therapy, the mutant cell lines showed approximately double the midkine at the cell membrane (p < 0.001) as well as elevated intracellular midkine (p < 0.01) (Figure 4H), suggesting that the D538G mutation increases both midkine production and retention at the cell surface.

The ESR1 D538G mutation, occurring in the ligand-binding domain of the estrogen receptor, confers constitutive activity and resistance to endocrine therapies. Our findings suggest that beyond its cell-autonomous effects, this mutation influences the broader tumor ecosystem through altered paracrine signaling. The increased MDK secretion observed in mutant cells may modulate immune cell recruitment, angiogenesis, and stromal remodeling, creating a microenvironment more conducive to tumor progression.

We repeated the ELISA and imaging analyses in T47D, a second breast cancer cell line, which did not recapitulate the MCF7 results (Supplemental Figure S4). This is consistent with previous work establishing that ESR1 mutation effects are context-specific^24,25^; MCF7 and T47D have different baseline ER expression levels and chromatin accessibility, leading to opposite downstream effects of the same mutation. The purpose of our wet lab validation is not to establish a universal mechanism across all breast cancer models, but to demonstrate that CITEgeist output can be directly input into external analytical tools and still yield findings confirmable in an isogenic experimental system, proving that CITEgeist deconvolution provides a robust foundation for translational follow-up.

### 2.5 Full Pipeline Application Reveals Endocrine Therapy Response Signatures Across a Patient Cohort

Having demonstrated CITEgeist's ability to generate testable biological hypotheses from individual cases, we applied the complete five-module pipeline across an expanded cohort to discover clinically relevant patterns at the population level.

Our cohort consisted of 12 specimens from 6 patients enrolled in a prospective clinical trial (NCT05914792)^14,15^ evaluating primary endocrine therapy (aromatase inhibitors) for women aged 70 and older with ER+/HER2- breast cancer who chose to forego upfront surgery. The study collected matched core biopsy and surgical specimens from each patient, yielding 12 specimens. Four patients were classified as responders based on imaging assessment and circulating tumor DNA (ctDNA) monitoring; these patients elected surgery around the 18-month mark for reasons unrelated to treatment failure, yielding 8 specimens. Two patients were classified as progressors, showing clinical or ctDNA-based evidence of disease progression around 36 months, at which point surgery was required, yielding 4 specimens. To preserve patient anonymity, all samples are labeled HCC22-088-P#-S#^15^.

Comparing cell type proportions from Module 3 between biopsy and surgical specimens revealed divergent remodeling trajectories by response status (Figure 5A). Progressors lost CD4 T cells (-30.7%) and gained cancer cells (+28.0% luminal in one patient, +11.7% basal in the other), consistent with immune evasion accompanying tumor expansion. Responders gained fibroblasts (+7.8% on average, with one patient showing a +31% fibroblast increase) and lost cancer cells (-5.6%), suggesting stromal remodeling and tumor regression under endocrine therapy.

Module 4 applied non-negative matrix factorization to the deconvolved cell type-specific expression layers from Module 3, discovering 490 gene expression programs across the cohort. We validated spatial coherence using Moran's I (Figure 5G): programs with significant positive spatial autocorrelation (I > 0.15, p < 0.01) exhibit spatial clustering, indicating biologically coordinated activity rather than technical noise. Across the cohort, 47% of discovered programs (231 of 490) exceeded the threshold for moderate spatial coherence (I > 0.15, p < 0.01). CD4 T cell programs showed the highest mean coherence (I = 0.33), followed by cancer and fibroblast programs (I = 0.25-0.26), while B cell and endothelial programs showed lower coherence (mean I = 0.08), reflecting their more dispersed spatial distributions.

Module 5 integrated these 490 programs across all 12 specimens, aligning them in a shared latent space using PCA followed by Harmony batch correction. Hierarchical clustering identified 65 aligned programs representing conserved transcriptional states (Figure 5B). Cancer_Basal programs were the most heterogeneous (13 aligned programs, mostly patient-specific), while aligned_032 (CD8 T cells) was perfectly conserved across all 12 samples. Module 5 also identified 163 conserved spatial relationships between program pairs; a network visualization (Figure 5D) reveals the multicellular organization architecture, with co-localized programs (bivariate Moran's I > 0.15), mutually exclusive programs (I < -0.15), and spatially independent programs.

Response analysis compared program enrichment between responder and progressor specimens, identifying 4 programs with significant differential enrichment (Figure 5C). One program was enriched in responders: a macrophage program (aligned_004) characterized by FABP4, HBA2, and TNXB. Three programs were enriched in progressors: a CD4 T cell program (aligned_008), a cancer luminal program (aligned_002), and a monocyte program (aligned_016).

The most striking finding emerged from bivariate spatial analysis of the conserved relationships, revealing a spatial organizational principle that fundamentally distinguishes responders from progressors and cannot be captured by bulk or single-cell RNA-seq. In the responding tumor P3-S2, fibroblast and CD4 T cell proportions were spatially co-localized (bivariate Moran's I = 0.358), with T cells concentrated in the same tissue regions as fibroblasts (Figure 5E). This co-localization suggests coordinated immune-stromal activity in responding tumors, where fibroblastic remodeling may facilitate T cell infiltration. In contrast, the progressing tumor P1-S2 showed spatial exclusion between cancer luminal cells and CD4 T cells (bivariate Moran's I = -0.194), with T cells absent from tumor-dense regions (Figure 5F). This exclusion pattern is consistent with immune evasion mechanisms in progressing tumors. Spatial proportion maps (Figure 5H) provide visual confirmation of these divergent architectures: responding tumors exhibit fibrotic scarring with reduced cancer burden, while progressing tumors show tumor expansion with immune exclusion. These spatial relationships represent precisely the type of microenvironment organization that is invisible to bulk RNA-seq (which destroys spatial information) and inaccessible to single-cell RNA-seq (which requires tissue dissociation).

To complement the spatial program analysis, we performed differential expression analysis using PyDESeq2^9,19^ on pseudo-bulk aggregated expression with design formula ~condition + timepoint. This identified 203 significantly differentially expressed genes (adjusted p < 0.05) between the 8 responder specimens (4 patients) and 4 progressor specimens (2 patients) across 13,350 tested genes. Of these, 120 genes were upregulated in responders, including NEDD9, TMC5, GP2, ZNF655, DACH1, and ACTG2, and 83 genes were upregulated in progressors, including VPREB3, COL4A4, CCNE1, FAM111B, and PLK1. Pathway enrichment of responder-upregulated genes highlighted epithelial-mesenchymal transition and TNF-alpha signaling via NF-kB, while progressor-upregulated genes were enriched for E2F targets and G2-M checkpoint pathways, consistent with proliferative advantage in progressing tumors (Supplementary Figure S2).

Critically, CITEgeist's outputs integrate directly with the standard computational biology ecosystem, enabling researchers to apply familiar methods to the enhanced data structures rather than learning new analysis paradigms. The cell type-specific expression layers from Module 3 serve as direct input for PyDESeq2^9^ differential expression, as demonstrated in the response analysis above. Gene lists from DE analysis or program loadings from Module 4 integrate with GSEApy for pathway enrichment. The spatial neighbor graphs and cell type proportions feed into COMMOT^27^ for cell-cell communication analysis, as demonstrated in the midkine case study. Standard scanpy workflows (Leiden clustering, UMAP, rank_genes_groups) operate directly on all CITEgeist outputs. This interoperability is not merely a convenience but a design requirement: by maintaining standard AnnData structures and preserving integer counts throughout the pipeline, CITEgeist ensures that its outputs satisfy the statistical assumptions underlying these downstream tools.

---

## Discussion

CITEgeist addresses a fundamental gap in spatial transcriptomics analysis: the mismatch between data that is inherently spatial and methods designed for non-spatial contexts. While established tools like Cell2Location^8^, RCTD, and Tangram^16^ have advanced the field, they share a common architecture of learning signatures from external scRNA-seq references and projecting them onto spatial data^18^. CITEgeist takes a different approach, using same-slide protein measurements as direct anchors at every analytical stage.

Many existing methods incorporate spatial information as regularization or post-processing, constraining an otherwise non-spatial algorithm with coordinates. CITEgeist inverts this relationship: spatial statistics drive discovery from the outset. Moran's I filters markers, spatial colocalization defines profiles, Laplacian smoothing regularizes proportions, and spatial coherence validates programs. This distinction has practical consequences. Reference-based methods can only discover patterns present in their reference; if a spatial cell state does not exist in dissociated single-cell data, it cannot be transferred. CITEgeist discovers patterns from the spatial data itself, using the protein layer as ground truth rather than an external atlas.

A practical advantage of this approach is the elimination of batch effects between reference and query data^3,4^. scRNA-seq atlases are collected from different individuals, processed with different protocols, and sequenced on different platforms than the spatial data being analyzed. These technical differences propagate through reference-based pipelines, potentially introducing systematic biases^7^. By using protein measurements from the identical tissue section, CITEgeist avoids this source of error. Our benchmarking on Xenium data demonstrates that same-slide proteomics provides sufficient information for competitive deconvolution accuracy without any reference data.

The midkine finding illustrates the translational potential of this framework. The identification of upregulated MDK signaling emerged from the analytical pipeline rather than being a predetermined hypothesis, demonstrating how computational tools designed with interpretability as a core principle can drive biological discovery. The ESR1 D538G mutation's influence on the broader tumor ecosystem through altered paracrine signaling, particularly via midkine, expands our understanding of how genomic alterations affect the tumor microenvironment beyond cell-autonomous effects^13^. CITEgeist's preservation of raw integer counts and standard AnnData output format enabled direct application of COMMOT^27^ for signaling analysis, followed by straightforward wet lab validation, a workflow that would not be possible with methods whose outputs are not compatible with existing analysis tools.

The spatial multi-omics field is evolving rapidly. Visium HD achieves 2-micrometer resolution approaching single-cell, while expanded antibody panels and spatial proteomics modalities (CODEX, IMC) provide increasingly rich protein information. CITEgeist's resolution-agnostic design, demonstrated on both spot-level Visium and single-cell Xenium, positions it for these developments. As protein panels expand, the colocalization-based profile discovery will leverage additional markers automatically; as resolution increases, the same modules apply with adjusted spatial neighborhood definitions.

---

## Limitations of the Study

CITEgeist's reliance on same-slide proteomics is both a strength and a limitation. Current spatial platforms support 30-50 antibodies, constraining the granularity of discoverable cell types; rare populations without specific protein markers may be underrepresented. As antibody panels expand, this limitation will diminish, but near-term applications may benefit from hybrid approaches combining same-slide proteomics with targeted reference information for specific rare populations.

Our clinical cohort is limited in size (n = 6 patients, with progressors from only 2 patients), and response-associated findings should be interpreted as hypothesis-generating rather than definitive. Although the differential expression analysis accounts for timepoint (biopsy vs surgical) as a covariate, the paired structure of multiple specimens per patient could not be modeled due to insufficient degrees of freedom. Validation in larger, independent cohorts is required to confirm the biological relevance of the response-associated programs reported here.

The T47D cell line did not recapitulate the midkine findings observed in MCF7, highlighting the context-dependency of ESR1 mutation effects across cell models^24,25^. Further work is needed to establish the full pathway supporting increased midkine secretion in ESR1-mutant breast tumors and to determine the generalizability of this mechanism.

Computationally, Module 3's quadratic programming optimization scales well to typical Visium datasets (~5,000 spots) but may require approximations for very large datasets such as Xenium with millions of cells. The current implementation processes such datasets in batches with checkpointing, but native scalability improvements would benefit large-scale applications.

---

## Acknowledgments

We thank Jingyang Qian for their guidance in utilizing the scCube method to simulate and benchmark spatially resolved Visium-like datasets from a single-cell reference. We also thank Jian Chen from the Lee-Oesterreich Lab for her help with RNA sequencing prep and ddPCR.

We acknowledge the assistance of the physicians and clinical support staff at the University of Pittsburgh Medical Center health system for their contributions to the clinical trial, with special thanks to Natera Inc. for their support with the clinical trial NCT05914792.

This research was partly supported by the University of Pittsburgh Center for Research Computing (RRID: SCR_022735) through the resources provided. We utilized the HTC cluster, supported by NIH award S10OD028483, and the H2P cluster, supported by NSF award OAC-2117681. Additionally, this work was funded by US NIH Grant R01HG010589, NIH Grant 5F30CA264963-03, and a Developmental Pilot Award from the UPMC Hillman Cancer Center, supported through NIH grant P30CA047904.

Histology sectioning was performed in the Pitt Biospecimen Core (RRID: SCR_025229), and the services and instruments used in this project were supported, in part, by the University of Pittsburgh and the Office of the Senior Vice Chancellor for Health Sciences. Work conducted at the UPMC Hillman Cancer Center Tissue and Research Pathology Services (TARPS) Shared Resource Facility was also supported, in part, by the University of Pittsburgh and the National Cancer Institute of the NIH under Award Number P30CA047904.

10x Visium library generation, Takara Smart Stranded Total RNA library generation, and Illumina sequencing were performed in the Health Sciences Sequencing Core (RRID: SCR_023116) at UPMC Children's Hospital, Rangos Research Center. These services and instruments were generously supported, in part, by the University of Pittsburgh, the Office of the Senior Vice Chancellor for Health Sciences, the Department of Pediatrics, the Institute for Precision Medicine, and the Richard K Mellon Foundation for Pediatric Research.

## Author Contributions

A.C.C. conceived the project, designed the computational framework, wrote the CITEgeist code, conducted the bioinformatic analyses, and prepared the initial manuscript draft. B.T.S. performed the CITEgeist benchmarking, simulation studies, and statistical analysis and contributed to discussing benchmarking and simulation strategies. N.C. facilitated the collection of clinical samples and coordinated the clinical trial. P.F.M. served as the clinical trial's principal investigator and performed breast surgery, collecting the samples. H.W. completed the wet lab validation of MDK signaling in ESR1-mutant cell lines. S.O. oversaw the overall project management and shaped the project's goals. R.S. contributed to manuscript revisions and provided expertise in mathematical modeling and computational methodologies. A.V.L. guided the project design and objectives. All authors reviewed and edited the manuscript.

## Declaration of Interests

The authors declare no competing interests.

---

## Ethics

The protocol for the related clinical trial was approved by the University of Pittsburgh IRB under STUDY21100091 and conducted through the UPMC Hillman Cancer Center under protocol 22-088.

---

## STAR Methods

### Resource Availability

**Lead Contact.** Further information and requests for resources should be directed to and will be fulfilled by the lead contact, Adrian V. Lee (leeav@upmc.edu).

**Materials Availability.** Establishment of MCF-7 and T47D mutant ESR1 cell lines was previously reported^13^. No new unique reagents were generated in this study.

**Data and Code Availability.** The CITEgeist package is available on GitHub at https://github.com/leeoesterreich/CITEgeist. A persistent version will be made available via FigShare at publication (doi: 10.6084/m9.figshare.28385675). The Wu et al. breast cancer atlas scRNA-seq data is available through GEO (accession: GSE176078). Visium datasets from the clinical trial will be deposited in GEO upon publication; reviewers can access the data via the link and private access token provided in the cover letter. Xenium benchmarking data is publicly available from 10x Genomics (Xenium Human Breast Cancer Dataset).

### Key Resources Table

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
|---|---|---|
| **Antibodies** | | |
| Anti-MDK (rabbit) | ThermoFisher | Cat# MA5-32538 |
| Anti-E-cadherin (mouse) | BD Biosciences | Cat# 610182 |
| Anti-rabbit Alexa Fluor 647 | Invitrogen | Cat# A32733 |
| Anti-mouse Alexa Fluor 555 | Invitrogen | Cat# A32727 |
| **Critical Commercial Assays** | | |
| Human MDK ELISA Kit | Invitrogen | Cat# EH319RB |
| **Experimental Models: Cell Lines** | | |
| MCF-7 (WT and D538G ESR1) | Li et al.^13^ | N/A |
| T47D (WT and D538G ESR1) | Li et al.^13^ | N/A |
| **Deposited Data** | | |
| Wu et al. breast cancer atlas | GEO | GSE176078 |
| CITEgeist code and simulations | FigShare | 10.6084/m9.figshare.28385675 |
| Xenium Human Breast Cancer | 10x Genomics | Publicly available |
| **Software and Algorithms** | | |
| CITEgeist | This paper | https://github.com/leeoesterreich/CITEgeist |
| Python | Python Software Foundation | v3.10 |
| Gurobi | Gurobi Optimization | v11.0.2 |
| scanpy | Wolf et al. | v1.10.4 |
| squidpy | Palla et al. | v1.6.2 |
| PyDESeq2 | Muzellec et al. | v0.4 |
| COMMOT | Cang et al.^27^ | N/A |
| GSEApy | Fang et al. | N/A |
| scCube | Qian et al.^12^ | v2.0.0 |
| CellProfiler | Carpenter et al. | N/A |
| ImageJ | Schneider et al. | N/A |

### Experimental Model and Study Participant Details

**Clinical Trial.** As part of a prospective, pragmatic, hybrid de-centralized non-randomized clinical trial (NCT05914792)^14^ evaluating ctDNA levels as a means to augment longitudinal monitoring of older patients with ER+ breast cancer receiving primary endocrine therapy, tumor tissue specimens were collected at baseline and at surgical intervention^15^. Patients aged 70 years and older with ER+/HER2- non-metastatic breast cancer signed informed consent (protocol approved by the University of Pittsburgh IRB under STUDY21100091, conducted through UPMC Hillman Cancer Center under protocol 22-088) and chose to forego upfront surgery in favor of primary endocrine therapy (aromatase inhibitors). A subset of six patients underwent surgery during the study follow-up, providing matched biopsy and surgical specimens. All samples are labeled HCC22-088-P#-S# to preserve patient anonymity.

**Cell Lines.** MCF-7 and T47D wild-type and CRISPR-edited D538G ESR1-mutant cell clones were evenly pooled and maintained in DMEM (ThermoFisher 11965118) or RPMI (ThermoFisher 11875119), respectively, supplemented with 10% fetal bovine serum and 1x penicillin-streptomycin (Sigma-Aldrich P4333) at 37C in a humidified incubator with 5% CO2. For hormone treatment experiments, cells were steroid-deprived using phenol red-free Improved Minimum Essential Medium (Richter's Mod.) (Fisher Scientific MT10026CV) with 10% charcoal-stripped serum. 17-beta-estradiol (Sigma-Aldrich E2758-1G) was used at 1 nM for +E2 conditions.

### Method Details

**Module 1: Marker Interest Detection.** For each protein marker, Moran's I spatial autocorrelation is calculated with significance assessed via permutation testing (1,000 permutations). Gaussian mixture modeling fits a 2-component GMM, and signal-to-noise ratio is computed as the ratio of the higher component's mean to the lower. Kurtosis is computed as the fourth standardized moment. Markers are considered interesting if they pass either Moran's I > 0.1 with p < 0.05, or kurtosis > 2.0, combined with GMM SNR > 1.5.

**Module 2: Profile Discovery.** Pairwise marker relationships are computed for all markers passing Module 1. Same-spot co-occurrence uses Jaccard similarity on binarized expression (threshold: 75th percentile). Adjacent-spot enrichment computes the average expression of marker B in neighbors of spots expressing marker A, compared to random expectation. Edges with permutation p < 0.05 for at least two statistics are retained. Hierarchical clustering uses Ward's method with Euclidean distance on the edge weight matrix. Profile count is selected by maximizing reconstruction accuracy (R^2^) while maintaining minimum profile size of 2 markers.

**Module 3: Spatial Deconvolution.** Pass 1 solves the quadratic program: minimize ||P*pi - y||^2^ + lambda * pi^T^ L pi, where P is the profile matrix, pi are non-negative proportions summing to 1, y is observed protein expression, and L is the graph Laplacian (default lambda = 0.1). Neighborhood finetuning adjusts proportions based on consistency with neighbor estimates using weighted averaging with distance decay. Pass 2 solves for gene expression: minimize ||G*pi - e||^2^ + lambda_g * ||G - G_prior||^2^, where G is cell type-specific expression, e is observed bulk expression, and G_prior is computed from proportion-weighted global expression patterns.

**Module 4: Program Discovery.** NMF factorizes the weighted expression matrix X = WH using sklearn's NMF with Frobenius norm loss and coordinate descent solver. Moran's I is computed on each program's activity pattern; programs with I > 0.15 and p < 0.01 are considered spatially coherent. Bivariate Moran's I is computed between all program pairs.

**Module 5: Cross-Sample Integration.** Programs are embedded using PCA (n_components = 30) on concatenated gene loadings. Harmony batch correction removes sample-specific variation (theta = 2.0, max_iter = 20). Aligned programs are identified via hierarchical clustering (Ward linkage) with cluster count selected by silhouette score optimization. Response analysis uses Mann-Whitney U test with Benjamini-Hochberg correction.

**Differential Expression Analysis.** Pseudo-bulk aggregation sums gene expression within each sample across 12 specimens (8 responder from 4 patients, 4 progressor from 2 patients). PyDESeq2 (v0.4)^9^ performs differential expression with design formula ~condition + timepoint, where condition is responder versus progressor and timepoint is biopsy (pre-treatment) versus surgical (post-treatment). The timepoint covariate accounts for systematic expression differences between pre- and post-treatment specimens. Genes with adjusted p < 0.05 are considered significant.

**Benchmarking.** Xenium single-cell data was aggregated into pseudo-Visium spots (55-micrometer diameter) to create ground truth. Cell type proportions were computed from single-cell annotations within each spot. Methods were evaluated on Pearson correlation, RMSE, MAE, and Jensen-Shannon divergence. Simulated spatial data was generated using scCube^12^ from a downsampled ER+ scRNA-seq atlas of 12,000 cells from 11 patients^11^.

**ELISA.** ELISA was performed using the Invitrogen Human MDK ELISA Kit (EH319RB) according to manufacturer instructions on conditioned media from MCF-7 and T47D WT and D538G cell lines.

**Immunofluorescence.** Cells were plated in 8-well chamber slides (FisherScientific 08-774-26) at 500 cells/well, steroid-deprived for 5 days with medium changed twice daily. For +E2 conditions, 1 nM E2 was added for 24 hours. Cells were fixed in 4% PFA, permeabilized in 0.2% Triton X-100, and blocked in 1% BSA, 10% goat serum, 0.1% Tween20 in PBS. Primary antibodies for anti-MDK (1:1,000) and anti-E-cadherin (1:200) were incubated overnight at 4C. Secondary antibodies (anti-rabbit Alexa Fluor 647, anti-mouse Alexa Fluor 555) were incubated for 1 hour at room temperature. Counterstaining was performed with Hoechst 33342. Cells were imaged on a Nikon AX confocal microscope with 10-image z-stacks at minimum 3 locations per sample at 60x oil magnification. Quantification was performed with CellProfiler and ImageJ.

### Quantification and Statistical Analysis

Spatial statistics (Moran's I, bivariate Moran's I) were assessed via permutation testing (1,000 permutations) with significance at p < 0.05. Benchmarking comparisons between methods used Wilcoxon signed-rank tests (n = 5 tissue regions). Differential expression used PyDESeq2 with Benjamini-Hochberg correction for multiple testing; genes with adjusted p < 0.05 were considered significant. ELISA comparisons used unpaired t-tests. Immunofluorescence quantification used Mann-Whitney U tests. Response analysis for program enrichment used Mann-Whitney U tests with Benjamini-Hochberg correction. All statistical analyses were performed in Python 3.10 using scipy and statsmodels.

---

## Supplementary Information

Supplementary File 1: Full mathematical derivation of the CITEgeist method.
Supplementary File 2: Tables reporting metrics of all benchmarked methods.
Supplementary File 3: Raw data from wet lab experiments (ELISA and immunofluorescence).
Supplementary File 4: Supplementary statistical and sequencing methods.
Supplementary Figure S1: ddPCR confirmation of ESR1 D538G mutation status.
Supplementary Figure S2: Differential expression volcano plot and pathway enrichment analysis.
Supplementary Figure S3: Runtime and computational performance benchmarking.
Supplementary Figure S4: T47D ELISA and immunofluorescence results.

---

## References

1. Huang, H., Wei, T., Zhang, A., et al. (2024). Trends in the incidence and survival of women with hormone receptor-positive breast cancer from 1990 to 2019: a large population-based analysis. Sci. Rep. *14*, 23690.

2. Patel, R., Klein, P., Tiersten, A., and Sparano, J.A. (2023). An emerging generation of endocrine therapies in breast cancer: a clinical perspective. Npj Breast Cancer *9*, 1-12.

3. Regan, C., and Preall, J. (2022). Practical Considerations for Single-Cell Genomics. Curr. Protoc. *2*, e498.

4. Chang, A.C.-C., Balic, M., Bartholow, T., et al. (2025). Hope for OTHERS (Our Tissue Helping Enhance Research & Science): research results from the University of Pittsburgh rapid autopsy program for breast cancer. Breast Cancer Res. *27*.

5. Wess, M., Andersen, M.K., Midtbust, E., et al. (2025). Spatial integration of multi-omics data from serial sections using the novel Multi-Omics Imaging Integration Toolset. GigaScience *14*, giaf035.

6. Zaidi, M., Fu, F., Cojocari, D., McKee, T.D., and Wouters, B.G. (2019). Quantitative Visualization of Hypoxia and Proliferation Gradients Within Histological Tissue Sections. Front. Bioeng. Biotechnol. *7*, 397.

7. Chen, J., Liu, W., Luo, T., et al. (2022). A comprehensive comparison on cell-type composition inference for spatial transcriptomics data. Brief. Bioinform. *23*, bbac245.

8. Kleshchevnikov, V., Shmatko, A., Dann, E., et al. (2022). Cell2location maps fine-grained cell types in spatial transcriptomics. Nat. Biotechnol. *40*, 661-671.

9. Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. *15*, 550.

10. Butler, A., Hoffman, P., Smibert, P., Papalexi, E., and Satija, R. (2018). Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. *36*, 411-420.

11. Wu, S.Z., Al-Eryani, G., Roden, D.L., et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. Nat. Genet. *53*, 1334-1347.

12. Qian, J., Bao, H., Shao, X., et al. (2024). Simulating multiple variability in spatially resolved transcriptomics with scCube. Nat. Commun. *15*, 5021.

13. Li, Z., McGinn, O., Wu, Y., et al. (2022). ESR1 mutant breast cancers show elevated basal cytokeratins and immune activation. Nat. Commun. *13*, 2011.

14. McAuliffe, P. (2025). Longitudinal ctDNA Monitoring in Older Women With ER+ Breast Cancer Who Forego Upfront Surgery in Favor of Primary Endocrine Therapy. clinicaltrials.gov. NCT05914792.

15. Carleton, N., Chang, A.C., Chen, F., et al. (2025). Longitudinal ctDNA Surveillance in Older Women with ER+ Breast Cancer to Facilitate Surgical De-Escalation: A Prospective, Hybrid-Decentralized Trial with Correlative Studies. medRxiv, 2025.08.23.25332468.

16. Biancalani, T., Scalia, G., Buffoni, L., et al. (2021). Deep learning and alignment of spatially resolved single-cell transcriptomes with Tangram. Nat. Methods *18*, 1352-1362.

17. Li, X., and Wang, C.-Y. (2021). From bulk, single-cell to spatial RNA sequencing. Int. J. Oral Sci. *13*, 36.

18. Li, H., Zhou, J., Li, Z., et al. (2023). A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics. Nat. Commun. *14*, 1548.

19. Anders, S., and Huber, W. (2010). Differential expression analysis for sequence count data. Genome Biol. *11*, R106.

20. Li, Z., Chen, F., Chen, L., et al. (2024). EstroGene2.0: A multi-omic database of response to estrogens, ER-modulators, and resistance to endocrine therapies in breast cancer. bioRxiv, 2024.06.28.601163.

21. Visvader, J.E. (2024). Midkine links aging with breast cancer--A new predictor of cancer risk. Cancer Cell *42*, 1815-1817.

22. Yan, P., Jimenez, E.R., Li, Z., et al. (2024). Midkine as a driver of age-related changes and increase in mammary tumorigenesis. Cancer Cell *42*, 1936-1954.e9.

23. Muramatsu, T. (2010). Midkine, a heparin-binding cytokine with multiple roles in development, repair and diseases. Proc. Jpn. Acad. Ser. B *86*, 410-425.

24. Arnesen, S., Blanchard, Z., Williams, M.M., et al. (2021). Estrogen receptor alpha mutations in breast cancer cells cause gene expression changes through constant activity and secondary effects. Cancer Res. *81*, 539-551.

25. Li, Z., Wu, Y., Yates, M.E., et al. (2022). Hotspot ESR1 Mutations Are Multimodal and Contextual Modulators of Breast Cancer Metastasis. Cancer Res. *82*, 1321-1339.

26. Guo, L., Kong, D., Liu, J., et al. (2023). Breast cancer heterogeneity and its implication in personalized precision therapy. Exp. Hematol. Oncol. *12*, 3.

27. Cang, Z., Zhao, Y., Almet, A.A., et al. (2023). Screening cell-cell communication in spatial transcriptomics via collective optimal transport. Nat. Methods *20*, 218-228.
