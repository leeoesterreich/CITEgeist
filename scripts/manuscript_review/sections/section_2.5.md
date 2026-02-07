# Section 2.5: Module 5 - Cross-Sample Integration Identifies Response-Associated Programs

**Referenced Figures**: Figure 5 (A, B, C, D, E)

## Section Text

To discover clinically relevant biology, Module 5 integrates programs across multiple patient samples, identifying conserved signatures and response-associated differences (Figure 5A). We applied this to a cohort of 14 breast cancer samples from 6 patients undergoing neoadjuvant therapy, including 5 samples from treatment responders and 9 from progressors.

Integration proceeds in three stages. First, program signatures (top genes per program) from all samples are embedded in a shared latent space using principal component analysis. Second, Harmony batch correction aligns programs across samples, removing patient-specific technical variation while preserving biological signal. Third, hierarchical clustering identifies groups of similar programs, which we term "aligned programs"—programs that appear consistently across patients.

From 590 individual programs across 14 samples, Module 5 identified 73 aligned programs representing conserved transcriptional states (Figure 5B). Of these, 5 programs appeared in more than 50% of samples, indicating highly conserved biology, while 36 were sample-specific. The aligned programs spanned all 7 cell types, with representation proportional to cell type abundance.

Response analysis compared program enrichment between responder and progressor samples (Figure 5C). Three programs were significantly enriched in responders: a macrophage program with FABP4/HBA2/TNXB, a macrophage program with VIM/TMSB4X/S100A6, and a cancer luminal program with MGP/MT-CO3/FOS. Four programs were enriched in progressors: a dendritic cell program with FTL/FN1/TIMP1, a cancer basal program with mitochondrial genes, a cancer luminal program with MGP/HBA2/LTF, and an endothelial program with KLHL17/MT-ND1/GPR101.

To validate these associations and identify specific genes driving response differences, we performed differential expression analysis using PyDESeq2 on pseudo-bulk aggregated expression (Figure 5D). Comparing 5 responder samples to 9 progressor samples across 13,371 tested genes, we identified 127 significantly differentially expressed genes (adjusted p < 0.05). The vast majority (122 genes) were upregulated in progressors, including matrix metalloproteinases (MMP13, MMP3, ADAMTS4, ADAMTS15), the necroptosis effector MLKL, and immune modulators (ALOX5AP, CLEC5A). Only 5 genes were upregulated in responders, including TMEM38B and TRIM72.

Module 5 also identifies conserved spatial relationships—program pairs that maintain consistent co-localization or exclusion across samples (Figure 5E). Of 191 conserved relationships detected, 26 (13.6%) showed consistent co-localization, 6 (3.1%) showed consistent mutual exclusion, and the remainder were independent. This relationship network highlights multicellular organization patterns that persist across patients despite inter-individual variation.

## Figure 5 Legend (for reference)

(A) Integration schematic. Programs from all samples (n = 14 breast cancer samples from 6 patients) are embedded in a shared latent space using PCA. Harmony batch correction aligns programs across samples, removing patient-specific technical variation. Hierarchical clustering identifies aligned programs—conserved transcriptional states appearing across patients. Response analysis compares program enrichment between responders (5 samples from patients P1, P5) and progressors (9 samples from patients P2-P4, P6).

(B) UMAP visualization of 590 individual programs from 14 samples, colored by cell type. Successful Harmony integration groups programs by biological identity rather than sample of origin. 73 aligned programs identified, with 5 appearing in >50% of samples.

(C) Response-associated programs. Bar plot showing 3 responder-enriched and 4 progressor-enriched aligned programs. Responder-enriched: macrophage programs with FABP4/HBA2/TNXB and VIM/TMSB4X/S100A6, cancer luminal program with MGP/MT-CO3/FOS. Progressor-enriched: dendritic cell program with FTL/FN1/TIMP1, cancer basal mitochondrial program, cancer luminal program with MGP/HBA2/LTF, endothelial program with KLHL17/MT-ND1/GPR101.

(D) PyDESeq2 differential expression volcano plot. Pseudo-bulk aggregation of deconvolved gene expression across 14 samples (5 responder, 9 progressor). 13,371 genes tested, 127 significant (adjusted p < 0.05). 122 genes upregulated in progressors (including MMP13, MMP3, ADAMTS4, ADAMTS15, MLKL, ALOX5AP, CLEC5A), 5 genes upregulated in responders (including TMEM38B, TRIM72). Top genes labeled.

(E) Conserved spatial relationship network. Of 191 conserved relationships detected across samples, 26 (13.6%) show consistent co-localization (green edges), 6 (3.1%) show consistent mutual exclusion (red edges), remainder are independent (gray). Node color indicates cell type. This network highlights multicellular organization patterns persisting across patients.
