# Section 2.3: Module 3 - Protein-Anchored Deconvolution with Spatial Regularization

**Referenced Figures**: Figure 3 (A, B, C)

## Section Text

Given cell type profiles from Module 2, Module 3 estimates cell type proportions at each spatial location and deconvolves bulk gene expression into cell type-specific layers (Figure 3A). Unlike reference-based methods that transfer signatures learned from external atlases, CITEgeist uses the protein measurements from the same tissue section as direct anchors.

Pass 1 formulates proportion estimation as a quadratic programming problem. For each spot, we find non-negative proportions that minimize reconstruction error between observed protein measurements and the profile-weighted mixture. To enforce spatial coherence, we add a Laplacian regularization term that penalizes differences between neighboring spots. This prevents noisy, salt-and-pepper proportion estimates while allowing genuine boundaries between tissue regions.

Following global optimization, neighborhood-aware finetuning adjusts proportions based on local context. For each spot, we examine the 6-8 nearest neighbors and their assignments, refining estimates where local evidence suggests corrections. This two-stage approach balances global consistency with local precision.

Pass 2 uses the estimated proportions to deconvolve gene expression. At each spot, bulk expression is modeled as a weighted sum of cell type-specific expression, where weights are the proportions from Pass 1. We solve for per-cell-type expression values using non-negative least squares with global and local enrichment priors derived from proportion-weighted spatial patterns. The result is a set of cell type-specific expression layers: for each gene at each spot, we obtain estimated expression within each cell type population.

We benchmarked CITEgeist against established deconvolution methods using Xenium data aggregated into pseudo-Visium spots where single-cell ground truth was available (Figure 3B). Across 5 tissue regions totaling 7,054 spots, CITEgeist achieved Pearson correlation r = 0.60 with ground truth proportions, comparable to reference-based methods Cell2Location (r = 0.61) and RCTD (r = 0.62). Jensen-Shannon divergence (lower is better) was 0.355 for CITEgeist versus 0.335 for Cell2Location and 0.347 for RCTD. Methods requiring label transfer performed substantially worse: Tangram achieved r = 0.14 and Seurat r = 0.17.

Importantly, CITEgeist achieved this accuracy without requiring any external reference data—the protein measurements from the same tissue section provided sufficient information for competitive deconvolution. This demonstrates that same-slide proteomics can replace reference atlases for proportion estimation while avoiding potential artifacts from reference-sample mismatches.

## Figure 3 Legend (for reference)

(A) Two-pass deconvolution schematic. Pass 1 estimates cell type proportions from protein markers using quadratic programming with Laplacian regularization to enforce spatial coherence across neighborhoods. Neighborhood-aware finetuning adjusts estimates based on local context. Pass 2 deconvolves bulk gene expression into cell type-specific layers using the estimated proportions as constraints, producing per-spot, per-cell-type expression matrices.

(B) Benchmarking on Xenium pseudo-Visium spots. Single-cell Xenium data was aggregated into pseudo-Visium spots (55 μm diameter) to create ground truth cell type proportions. Methods compared: CITEgeist (same-slide proteomics, no reference), Cell2Location (reference-based), RCTD (reference-based), Tangram (label transfer), Seurat (label transfer). Metrics: Pearson correlation with ground truth (higher is better), Jensen-Shannon divergence (lower is better), RMSE (lower is better). CITEgeist achieves r = 0.60, comparable to reference-based methods (Cell2Location r = 0.61, RCTD r = 0.62), without requiring external reference data. n = 7,054 spots across 5 tissue regions.

(C) Example spatial visualization of estimated cell type proportions across a tissue region, showing smooth spatial patterns consistent with tissue architecture.
