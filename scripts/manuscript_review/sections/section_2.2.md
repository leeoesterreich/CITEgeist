# Section 2.2: Modules 1-2 - Automated Profile Discovery from Spatial Colocalization

**Referenced Figures**: Figure 2 (A, B, C, D)

## Section Text

A key advantage of same-slide proteomics is the ability to discover cell type profiles directly from the spatial data, rather than relying on predefined marker lists. CITEgeist's Modules 1-2 implement an automated discovery pipeline that identifies which protein markers are informative and groups them into biologically coherent cell type profiles (Figure 2A-B).

Module 1 applies three statistical gates to filter the ~30-50 antibodies in a typical CITE-seq panel. Moran's I spatial autocorrelation (I > 0.1, p < 0.05) identifies markers with organized spatial patterns—randomly distributed markers provide little information about tissue architecture. Gaussian mixture modeling fits a two-component distribution to each marker, and markers with signal-to-noise ratio below 1.5 are excluded. Kurtosis analysis (κ > 2.0) identifies markers with peaked distributions, which often indicate cell type-specific expression. Markers passing either the Moran's I gate or the kurtosis gate (combined with the GMM filter) are considered "interesting" and forwarded to Module 2.

Module 2 constructs a spatial colocalization network from pairwise marker relationships. For each marker pair, we compute: (1) same-spot co-occurrence via Jaccard similarity on binarized expression; (2) Pearson correlation on continuous values; (3) adjacent-spot enrichment, measuring whether one marker's expression predicts the other in neighboring spots; and (4) bivariate Moran's I, quantifying spatial co-patterning. Edges passing significance thresholds (permutation-based p < 0.05) are retained, weighted by combined evidence.

Connected components in the filtered graph represent distinct cellular lineages—epithelial markers cluster separately from immune markers, for example. Within each component, hierarchical clustering with dynamic tree cutting extracts individual profiles. The optimal number of profiles is determined by reconstruction accuracy: we select the partition that best explains the observed protein measurements while maintaining non-redundancy.

To validate the approach, we applied Modules 1-2 to Xenium single-cell data where ground truth cell types were available from RNA-based clustering (Figure 2C). The discovered protein profiles correctly grouped known markers: CD3E and CD8A clustered together (T cells), CD68 and CD163 clustered together (macrophages), and EPCAM and KRT markers clustered together (epithelial cells). This demonstrates that spatial colocalization provides sufficient signal to recover biologically meaningful cell type definitions without external reference data (Figure 2D).

## Figure 2 Legend (for reference)

(A) Module 1: Marker interest detection using three statistical gates. Moran's I spatial autocorrelation identifies markers with organized spatial patterns (threshold I > 0.1, p < 0.05). Gaussian mixture modeling separates signal from background, computing signal-to-noise ratios (threshold SNR > 1.5). Kurtosis analysis detects peaked distributions indicative of localized expression (threshold κ > 2.0). Markers passing either the Moran's I gate or the kurtosis gate (combined with the GMM filter) are forwarded to Module 2.

(B) Module 2: Profile assembly from spatial colocalization. For each marker pair, we compute same-spot co-occurrence, expression correlation, adjacent-spot enrichment, and bivariate Moran's I. Significance-filtered relationships (permutation p < 0.05) form a colocalization network where connected components represent distinct cellular lineages. Hierarchical clustering with dynamic tree cutting extracts marker profiles, with optimal count determined by reconstruction accuracy.

(C) Xenium single-cell demonstration: Modules 1-2 applied to single-cell resolution data, showing that spatial colocalization analysis discovers biologically coherent profiles without reference data.

(D) Validation: discovered profiles correctly recover known cell type markers. CD3E and CD8A cluster together (T cells), CD68 and CD163 cluster together (macrophages), EPCAM and KRT markers cluster together (epithelial cells). This demonstrates that spatial colocalization provides sufficient signal for cell type definition.
