# Section 2.1: Framework Overview

**Referenced Figures**: Figure 1 (A, B, C)

## Section Text

CITEgeist is organized as a modular pipeline where each stage leverages spatial information and passes interpretable outputs to subsequent modules (Figure 1A). The framework accepts spatial transcriptomics data with matched protein measurements in standard AnnData format, supporting both spot-level (Visium) and single-cell (Xenium, CosMx) resolutions.

**Module 1: Marker Interest Detection** evaluates each protein marker for spatial organization using three complementary statistics. Moran's I spatial autocorrelation identifies markers with non-random spatial patterns. Gaussian mixture modeling separates signal from background, computing signal-to-noise ratios. Kurtosis analysis detects peaked distributions indicative of localized expression. Markers passing these gates are forwarded to profile discovery.

**Module 2: Profile Assembly** discovers cell type marker profiles through spatial colocalization analysis. For each marker pair, the module computes same-spot co-occurrence, expression correlation, adjacent-spot enrichment, and bivariate Moran's I. Significance-filtered relationships are assembled into a colocalization network, where connected components represent distinct lineages. Hierarchical clustering with dynamic tree cutting extracts marker profiles, which are validated against known cell type associations.

**Module 3: Spatial Deconvolution** performs two-pass optimization. Pass 1 estimates cell type proportions using protein markers as anchors, applying Laplacian smoothing to enforce spatial coherence across neighborhoods. Neighborhood-aware finetuning adjusts estimates based on local context. Pass 2 deconvolves bulk gene expression into cell type-specific layers using the estimated proportions as constraints, producing per-spot, per-cell-type expression matrices.

**Module 4: Program Discovery** applies non-negative matrix factorization (NMF) to each cell type-specific expression layer, discovering gene expression programs within each population. Moran's I validates that discovered programs exhibit spatial coherence—programs with high spatial autocorrelation represent biologically meaningful spatial organization rather than technical noise. Bivariate Moran's I between program pairs identifies co-localized and mutually exclusive relationships.

**Module 5: Cross-Sample Integration** aligns programs across multiple samples using Harmony-based batch correction in a shared latent space. Hierarchical clustering identifies conserved programs appearing across patients. Response analysis compares program enrichment between clinical groups (e.g., responders vs progressors), enabling discovery of clinically relevant signatures.

Throughout the pipeline, outputs are stored in standard formats compatible with scanpy, enabling seamless integration with established analysis tools (Figure 1B-C).

## Figure 1 Legend (for reference)

(A) Schematic of the five-module pipeline: Module 1 (Marker Interest Detection), Module 2 (Profile Assembly), Module 3 (Spatial Deconvolution), Module 4 (Program Discovery), Module 5 (Cross-Sample Integration). Arrows indicate data flow between modules. Input data consists of spatial transcriptomics with matched protein measurements (CITE-seq antibody capture). Each module produces interpretable outputs stored in standard AnnData format.

(B) Spatial-native operations highlighted at each stage: Moran's I spatial autocorrelation for marker detection (Module 1), colocalization networks for profile discovery (Module 2), Laplacian regularization for proportion estimation (Module 3), spatial coherence validation for programs (Module 4), and conserved relationship analysis across samples (Module 5).

(C) Resolution flexibility: the same framework operates on spot-level data (Visium, ~55 μm diameter, 1-10 cells per spot) and single-cell resolution data (Xenium, CosMx). Module implementations automatically adapt to data resolution through spatial neighbor definitions that scale with technology.
