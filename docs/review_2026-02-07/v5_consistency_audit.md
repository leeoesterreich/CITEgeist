# CITEgeist CRM v5 Manuscript -- Internal Consistency Audit

**Audit Date**: 2026-02-08
**File Audited**: `manuscript/CITEgeist_CRM_v5.md`
**Compared Against**: `manuscript/CITEgeist_Patterns_v4.md`, `manuscript/CITEgeist_Patterns_v4_Figures.md`, `docs/review_2026-02-07/verified_numbers.md`

**Severity Key**:
- **CRITICAL**: Factual error or internal contradiction that would mislead readers
- **WARNING**: Inconsistency between sections or between v4/v5 that needs resolution
- **INFO**: Stylistic, minor, or cosmetic issue

---

## A. Clinical Trial Details Consistency

### A1. Trial ID: NCT05914792

| Location | v5 Line | Text | Status |
|----------|---------|------|--------|
| Summary | L20 | "NCT05914792" | OK |
| Section 2.4 | L91 | "NCT05914792" | OK |
| Section 2.5 | L107 | "NCT05914792" | OK |
| STAR Methods (Clinical Trial) | L218 | "NCT05914792" | OK |
| Acknowledgments | L151 | "NCT05914792" | OK |

**PASS** -- Trial ID is consistent across all 5 occurrences.

### A2. Patient Population Description

| Location | v5 Line | Description |
|----------|---------|-------------|
| Summary | L20 | "primary endocrine therapy trial for older women with ER+ breast cancer" |
| Section 2.4 | L91 | "primary endocrine therapy (aromatase inhibitors) in older women with ER+ breast cancer who forego upfront surgery" |
| Section 2.5 | L107 | "primary endocrine therapy (aromatase inhibitors) for women aged 70 and older with ER+/HER2- breast cancer who chose to forego upfront surgery" |
| STAR Methods | L218 | "Patients aged 70 years and older with ER+/HER2- non-metastatic breast cancer...chose to forego upfront surgery in favor of primary endocrine therapy (aromatase inhibitors)" |

- **INFO (A2a)**: The Summary (L20) says "older women" but does not specify age >= 70 or HER2- status. The Methods (L218) adds "non-metastatic." These are acceptable levels of detail for Summary vs Methods, but for completeness, the Summary could say "ER+/HER2-".
- **INFO (A2b)**: v4 (L107 in v4) called these patients "undergoing neoadjuvant therapy." v5 correctly changed this to "primary endocrine therapy," which is clinically distinct. This is a corrected improvement, not an inconsistency.

### A3. Cohort Numbers (6 patients, 14 specimens)

| Location | v5 Line | Patients | Specimens | Responders | Progressors |
|----------|---------|----------|-----------|------------|-------------|
| Section 2.5 | L107 | 6 | 14 | 2 patients, 5 specimens | 4 patients, 9 specimens |
| Limitations | L139 | "n = 6 patients, with responders from only 2 patients" | -- | -- | -- |
| STAR Methods | L218 | "six patients" | -- | -- | -- |
| Fig 5A legend (v4 Figures) | L71 | 6 patients | 14 samples | 5 samples (P1, P5) | 9 samples (P2-P4, P6) |
| DE analysis | L113 | -- | -- | "5 responder" | "9 progressor specimens" |

**PASS** -- All patient/specimen numbers are internally consistent (6 patients, 14 specimens, 2 responders yielding 5 specimens, 4 progressors yielding 9 specimens).

- **WARNING (A3a)**: v4 Abstract (L15) said "14-patient breast cancer cohort" -- this was **14 samples from 6 patients**, not 14 patients. v5 correctly states "14 specimens from 6 patients" (L107). However, v4 Section 2.5 (L107) also said "14 breast cancer samples from 6 patients." The v4 abstract error has been fixed in v5. No residual issue in v5.

### A4. Response Criteria

| Location | v5 Line | Responder Criteria | Progressor Criteria |
|----------|---------|-------------------|---------------------|
| Section 2.5 | L107 | "imaging assessment and circulating tumor DNA (ctDNA) monitoring...elected surgery around the 18-month mark for reasons unrelated to treatment failure" | "clinical or ctDNA-based evidence of disease progression around 36 months, at which point surgery was required" |
| STAR Methods | L218 | Not explicitly restated | Not explicitly restated |

**PASS** -- Response criteria stated once in detail (Section 2.5 L107); Methods section describes the trial broadly. Consistent.

### A5. Sample Labeling Convention

| Location | v5 Line | Text |
|----------|---------|------|
| Section 2.5 | L107 | "HCC22-088-P#-S#" |
| STAR Methods | L218 | "HCC22-088-P#-S#" |

**PASS** -- Consistent in both occurrences.

### A6. Treatment Description

| Location | v5 Line | Text |
|----------|---------|------|
| Summary | L20 | "primary endocrine therapy trial" |
| Section 2.4 | L91 | "primary endocrine therapy (aromatase inhibitors)" |
| Section 2.5 | L107 | "primary endocrine therapy (aromatase inhibitors)" |
| STAR Methods | L218 | "primary endocrine therapy (aromatase inhibitors)" |

**PASS** -- Consistent; Summary omits "(aromatase inhibitors)" which is acceptable for brevity.

---

## B. Benchmarking Numbers

### B1. Xenium Proportion Deconvolution

**Verified against `docs/review_2026-02-07/verified_numbers.md`**

| Metric | v5 (L85) | v4 (L81-83) | Verified Source | Status |
|--------|----------|-------------|-----------------|--------|
| CITEgeist r | 0.60 +/- 0.05 | 0.60 +/- 0.05 | 0.5997 -> 0.60 | OK |
| Cell2Location r | 0.61 +/- 0.04 | 0.61 +/- 0.04 | 0.6115 -> 0.61 | OK |
| RCTD r | 0.62 +/- 0.03 | 0.62 +/- 0.03 | 0.6185 -> 0.62 | OK |
| Tangram r | 0.14 +/- 0.08 | 0.14 +/- 0.08 | 0.1432 -> 0.14 | OK |
| Seurat r | 0.17 +/- 0.07 | 0.17 +/- 0.07 | 0.1729 -> 0.17 | OK |
| CITEgeist RMSE | 0.167 +/- 0.006 | 0.167 +/- 0.006 | 0.167 +/- 0.006 | OK |
| CITEgeist MAE | 0.115 +/- 0.005 | 0.115 +/- 0.005 | 0.115 +/- 0.005 | OK |
| Cell2Location RMSE | 0.179 +/- 0.017 | 0.179 +/- 0.017 | 0.179 +/- 0.017 | OK |
| RCTD RMSE | 0.177 +/- 0.004 | 0.177 +/- 0.004 | 0.177 +/- 0.004 | OK |
| Spots | 7,054 | 7,054 | 7,054 | OK |
| Regions | 5 | 5 | 5 | OK |
| Wilcoxon p vs C2L | 0.31 | 0.31 | -- | OK |
| Wilcoxon p vs RCTD | 0.19 | 0.19 | -- | OK |

**PASS** -- All Xenium proportion benchmarking numbers match between v5, v4, and verified source data.

### B2. Number of Cell Types in Xenium Benchmarking

- **WARNING (B2a)**: The `verified_numbers.md` file header says "6 Cell Types, 5 Regions." However, v5 L87 says "mean Pearson r = 0.44 across **7 cell types**" for the GEX benchmarking. The v4 manuscript also says "7 cell types" (L87, L99 in v4). The v4 figure legend (Figure 4C, L59) says "7 cell types, 175 total programs."

  This discrepancy (6 vs 7 cell types) exists between the `verified_numbers.md` header and the manuscript text. The 6 in verified_numbers.md likely refers to proportion deconvolution cell types, while 7 in the manuscript may include an additional cell type for GEX evaluation. This needs clarification -- if the Xenium data truly has 7 ground truth cell types for GEX and 6 for proportions, the manuscript should state this explicitly. Currently, the manuscripts only say "7 cell types" in GEX context, which is internally consistent within the manuscript, but the proportion benchmarking section does not state a cell type count.

  **Severity: WARNING** -- Ambiguity about whether it is 6 or 7 cell types for different evaluations; if 7 is correct throughout, update verified_numbers.md.

### B3. JSD Values

v5 (L85) does not explicitly report JSD for the Xenium benchmarking in the Results text. v4 (L81) reported: CITEgeist JSD = 0.355 +/- 0.01, Cell2Location JSD = 0.335 +/- 0.03, RCTD JSD = 0.347 +/- 0.01.

- **INFO (B3a)**: v5 dropped JSD values from the Results prose, which is acceptable given the focus on r and RMSE. They are available in verified_numbers.md. No inconsistency.

### B4. GEX Benchmarking Numbers

| Metric | v5 (L87) | v4 (L87) | Status |
|--------|----------|----------|--------|
| CITEgeist GEX r | 0.44 | 0.44 | OK |
| Gene coverage | 100%, 395 genes | 100%, 395 genes | OK |
| scResolve r | 0.43 | 0.43 | OK |
| scResolve coverage | 32% of tissue locations | 32% spot coverage | OK |

**PASS** -- GEX benchmarking numbers consistent between v4 and v5.

### B5. Simulated Benchmarking Numbers

| Metric | v5 (L83) | v4 (L85-86) | Status |
|--------|----------|-------------|--------|
| CITEgeist simulated r | 0.95 | 0.95 | OK |
| CITEgeist simulated JSD | 0.16 | 0.16 | OK |
| CITEgeist simulated RMSE | 0.08 (both conditions) | Not stated in same way | OK |
| C2L RMSE degradation | 0.08 -> 0.167 | Not in v4 | NEW |
| Seurat RMSE degradation | 0.10 -> 0.133 | Not in v4 | NEW |
| RCTD reference sensitivity | 0.05 -> 0.21 (4-fold) | Not in v4 | NEW |
| Tangram JSD in segmented | 0.56 | Not in v4 | NEW |
| GEX simulated CITEgeist NRMSE | 0.04 | Not in v4 | NEW |
| GEX simulated C2L NRMSE | 0.16 | Not in v4 | NEW |
| GEX simulated Tangram NRMSE | 0.20 | Not in v4 | NEW |

- **INFO (B5a)**: v5 adds substantial new simulation benchmarking detail not in v4 (reference quality sensitivity, tissue type performance). These are new claims that cannot be cross-verified against v4 but are internally consistent within v5.

### B6. Simulated Data Description

| Metric | v5 (L83, L236) | v4 (not stated) | Status |
|--------|----------------|-----------------|--------|
| Simulation tool | scCube | scCube | OK |
| Reference atlas | Wu et al. ER+ scRNA-seq | Wu et al. | OK |
| Atlas size | "12,000 cells from 11 patients" | Not stated | -- |
| Comprehensive reference | 30,000 cells | Not in v4 | NEW |
| Reduced reference | 8,000 cells (single-sample) | Not in v4 | NEW |

- **WARNING (B6a)**: The Methods (L236) says "downsampled ER+ scRNA-seq atlas of 12,000 cells from 11 patients" but the Results (L83) says "comprehensive (30,000-cell) versus realistic single-sample (8,000-cell) references." The 12,000-cell figure in Methods appears to be the size used for CITEgeist's own simulation, while 30,000 and 8,000 are the sizes of the external references given to reference-based methods. However, this creates ambiguity -- the Methods should clarify that 12,000 cells refers to the scCube simulation input, while 30,000 and 8,000 are the reference atlas sizes for competing methods.

  **Severity: WARNING** -- Potential reader confusion about which cell count applies to which purpose.

---

## C. Module 4+5 Numbers

### C1. Program Counts (175 vs 590)

| Context | v5 Line | Number | Source |
|---------|---------|--------|--------|
| Module 4 on Xenium (validation) | Not explicitly in v5 | 175 | v4 Figure 4C legend (L59) |
| Module 4+5 full cohort | L111 | 590 individual programs | v5 Section 2.5 |
| Aligned programs | L111 | 73 | v5 Section 2.5 |

- **INFO (C1a)**: v5 does not mention "175 programs" anywhere in the main text. The 175 number appears only in the v4 figures document (Figure 4C legend). In v5, the Xenium program discovery validation is less explicitly numbered. The full cohort number of 590 is stated clearly. These are different contexts (Xenium validation = 175, full cohort = 590) and are not contradictory.

- **WARNING (C1b)**: v5 Results Section 2.5 (L109) describes Module 4 results on the full cohort with "68% of discovered programs exceeded the threshold for moderate spatial coherence." The v4 Figure 4C legend says "57% of discovered programs (100/175) exceed this threshold" for the Xenium validation. These percentages apply to different datasets (full cohort vs Xenium), but the v5 text does not clearly delineate this -- Section 2.5 describes the full cohort results, while the "68%" figure comes from the cohort context. This is consistent within v5, though readers familiar with v4 may wonder about the change from 57%.

### C2. Spatial Coherence Threshold Discrepancy (Module 4)

| Location | v5 Line | Threshold Used | Context |
|----------|---------|----------------|---------|
| Results 2.5 (cohort) | L109 | I > 0.2, p < 0.01 | Full cohort "68% exceeded threshold" |
| Methods (Module 4) | L230 | I > 0.15, p < 0.01 | Methods description |
| v4 Results (Xenium) | v4 L97 | I > 0.2, p < 0.01 | "significant positive Moran's I" |
| v4 Figure 4C | v4 Figs L59 | I > 0.15 | "Threshold line at I = 0.15; 57% exceed" |
| v4 Methods | v4 L209 | I > 0.15, p < 0.01 | Methods description |

- **CRITICAL (C2a)**: The v5 **Methods** (L230) says programs with "I > 0.15 and p < 0.01 are considered spatially coherent," but the v5 **Results** (L109) says "programs with significant positive spatial autocorrelation (I > 0.2, p < 0.01) exhibit spatial clustering" and "68% of discovered programs exceeded the threshold for moderate spatial coherence."

  **The Methods and Results use different Moran's I thresholds (0.15 vs 0.2).** If the threshold is 0.15, the 68% figure would apply to that cutoff. If the threshold is 0.2, a different (lower) percentage should apply. These cannot both be correct simultaneously -- either the threshold is 0.15 (update the Results text) or 0.2 (update the Methods text), and the percentage must match the actual threshold used.

  Additionally, the v4 Figure 4C legend says "57% (100/175) exceed [I > 0.15]" on Xenium data, while v5 says "68% exceeded the threshold" on the full cohort. If the full cohort uses I > 0.2 as stated in v5 Results, then 68% exceeding I > 0.2 would be a remarkably high rate; normally a stricter threshold produces a lower pass rate. Alternatively, if 68% is the proportion exceeding I > 0.15 on the full cohort, then the Results text should say I > 0.15, matching Methods.

  **Severity: CRITICAL** -- Methods/Results threshold mismatch (0.15 vs 0.2) makes the "68%" figure ambiguous.

### C3. Fibroblast Mean Moran's I

| Location | v5 Line | Value |
|----------|---------|-------|
| Section 2.5 (cohort) | L109 | mean I = 0.28 |
| v4 Section 2.4 / Figure 4B | v4 L99/Figs L58 | mean I = 0.28 |

**PASS** -- Fibroblast mean Moran's I = 0.28 is consistent across v4 and v5.

### C4. T Cell Moran's I Range

| Location | v5 Line | Value |
|----------|---------|-------|
| Section 2.5 | L109 | I range 0.08-0.26 |
| v4 Section 2.4 | v4 L99 | I range 0.08-0.26 |

**PASS** -- Consistent.

### C5. Module 5 Integration Numbers

| Metric | v5 (L111, L115) | v4 (L111, L117) | Figs (v4) | Status |
|--------|-----------------|-----------------|-----------|--------|
| Individual programs | 590 | 590 | 590 (Fig 5B) | OK |
| Aligned programs | 73 | 73 | 73 (Fig 5B) | OK |
| Programs in >50% samples | 5 | 5 | 5 (Fig 5B) | OK |
| Specimen-specific programs | 36 | 36 | -- | OK |
| Conserved relationships | 191 | 191 | 191 (Fig 5) | OK |
| Co-localization | 26 (13.6%) | 26 (13.6%) | 26 (13.6%) | OK |
| Mutual exclusion | 6 (3.1%) | 6 (3.1%) | 6 (3.1%) | OK |

**PASS** -- All Module 5 integration numbers match across v4, v5, and figure legends.

### C6. Response-Associated Programs

| Category | v5 (L113) | v4 (L113) | Status |
|----------|----------|----------|--------|
| Responder-enriched | 3 programs | 3 programs | OK |
| Progressor-enriched | 4 programs | 4 programs | OK |
| Macrophage (responder) | FABP4, HBA2, TNXB | FABP4/HBA2/TNXB | OK |
| Cancer luminal (responder) | MGP, MT-CO3, FOS | MGP/MT-CO3/FOS | OK |
| DC (progressor) | FTL, FN1, TIMP1 | FTL/FN1/TIMP1 | OK |

- **INFO (C6a)**: v5 mentions only 2 of the 3 responder-enriched programs by gene list (macrophage + cancer luminal), omitting the second macrophage program (VIM/TMSB4X/S100A6 from v4 L113). The v4 figure legend lists all three. v5 says "including" which is acceptable, not claiming to be exhaustive.

**PASS** within v5; one program name omitted from prose is acceptable with "including" language.

---

## D. Midkine Results

### D1. ELISA

| Claim | v5 (L97) | Status |
|-------|----------|--------|
| "approximately double the midkine" | OK -- stated for both E2-deprived and E2-supplemented | OK |
| p < 0.0001 | L97 | OK |
| Cell line | MCF7 (WT vs D538G) | OK |
| Conditions | "estrogen-deprived and estradiol-supplemented" | OK |

**PASS** -- ELISA claims consistent.

### D2. Immunofluorescence

| Claim | v5 (L97) | Status |
|-------|----------|--------|
| "approximately double the midkine at the cell membrane" | OK | OK |
| Membrane MDK p-value | p < 0.001 | OK |
| Intracellular MDK p-value | p < 0.01 | OK |
| Condition | "estrogen-deprived environment" | OK |

**PASS** -- IF claims consistent.

### D3. COMMOT Validation

| Claim | v5 (L95) | Status |
|-------|----------|--------|
| Spearman rho = 0.46 | OK | OK |
| p = 2.86e-37 | OK | OK |
| Comparison: "CITEgeist-derived secretion signals correlated with...Human Protein Atlas" | OK | OK |
| Signaling: MDK, PTN, MIF | OK | OK |

**PASS** -- COMMOT validation statistics consistent.

### D4. T47D Negative Result

| Claim | v5 (L101) | Status |
|-------|-----------|--------|
| "did not recapitulate the MCF7 results" | OK | OK |
| Supplemental Figure S3 reference | OK | OK |
| Context-specificity refs [24,25] | OK | OK |

**PASS** -- Negative result properly described.

---

## E. Differential Expression Analysis Numbers

| Metric | v5 (L113) | v4 (L115) | Figs (v4 Fig 5D) | Status |
|--------|----------|----------|-------------------|--------|
| Genes tested | 13,371 | 13,371 | 13,371 | OK |
| Significant (adj p < 0.05) | 127 | 127 | 127 | OK |
| Upregulated in progressors | 122 | 122 | 122 | OK |
| Upregulated in responders | 5 | 5 | 5 | OK |
| MMP13, MMP3 mentioned | Yes | Yes | Yes | OK |
| ADAMTS4, ADAMTS15 mentioned | Yes | Yes | Yes | OK |
| MLKL mentioned | Yes | Yes | Yes | OK |
| ALOX5AP, CLEC5A mentioned | Yes | Yes | Yes | OK |
| TMEM38B, TRIM72 mentioned | Yes | Yes | Yes | OK |

**PASS** -- All DE analysis numbers are consistent across v4, v5, and figure legends.

- **INFO (E1a)**: The sum 122 + 5 = 127 checks out mathematically.

---

## F. Methods-Results Cross-Check

### F1. Module 1 Thresholds

| Parameter | Methods (L224) | Results (L73) | Status |
|-----------|----------------|---------------|--------|
| Moran's I threshold | > 0.1, p < 0.05 | > 0.1, p < 0.05 | OK |
| GMM SNR | > 1.5 | "below 1.5" excluded | OK |
| Kurtosis | > 2.0 | kappa > 2.0 | OK |
| Gate logic | "either Moran's I or kurtosis, combined with GMM" | "either...or...combined with the GMM filter" | OK |

**PASS** -- Module 1 thresholds match between Methods and Results.

### F2. Module 3 Lambda

| Parameter | Methods (L228) | Results (L81) | Status |
|-----------|----------------|---------------|--------|
| Laplacian lambda | 0.1 | "Laplacian regularization term" (no value in Results) | OK |

**PASS** -- Lambda specified in Methods, referenced conceptually in Results.

### F3. Module 4 Moran's I Threshold

See **CRITICAL (C2a)** above. Methods says I > 0.15; Results Section 2.5 says I > 0.2.

**CRITICAL** -- Already flagged in C2a.

### F4. Module 2 Significance Threshold

| Parameter | Methods (L226) | Results (L75) | Status |
|-----------|----------------|---------------|--------|
| Edge retention | "permutation p < 0.05 for at least two statistics" | "permutation-based p < 0.05" | OK |

- **INFO (F4a)**: Results (L75) says "Edges passing significance thresholds (permutation-based p < 0.05) are retained, weighted by combined evidence," while Methods (L226) says "Edges with permutation p < 0.05 for at least two statistics are retained." The Methods adds the "at least two statistics" criterion that the Results omit. This is not contradictory (Results is less specific), but the Methods provides a more restrictive criterion. Readers relying on Results alone would have an incomplete picture.

  **Severity: INFO** -- Methods is more detailed than Results, which is normal for journal papers.

### F5. Statistical Tests

| Analysis | Methods Used (STAR Methods) | Reported in Results | Status |
|----------|-----------------------------|---------------------|--------|
| Benchmarking comparison | Wilcoxon signed-rank (L244) | Wilcoxon signed-rank (L85) | OK |
| DE analysis | PyDESeq2 BH correction (L234) | PyDESeq2 adj p < 0.05 (L113) | OK |
| ELISA | Unpaired t-tests (L244) | Not specified in Results | INFO |
| IF | Mann-Whitney U (L244) | Not specified in Results | INFO |
| Response programs | Mann-Whitney U + BH (L232, L244) | Not specified in Results | INFO |

- **INFO (F5a)**: Statistical test types for ELISA, IF, and response analysis are described in STAR Methods/Quantification but not explicitly named in Results text. This is standard journal practice.

**PASS** -- Statistical methods are consistent between the sections that describe them.

---

## G. Author List and Affiliations

### G1. Author Comparison (v4 vs v5)

**v4 (L3)**:
```
Alexander Chih-Chieh Chang^1,2^, Brent T. Schlegel^1^, Neil Carleton^1^,
Priscilla F. McAulife^1^, Steffi Oesterreich^1,2^, Russell Schwartz^2,3^,
Adrian V. Lee^1,2,*^
```

**v5 (L3)**:
```
Alexander Chih-Chieh Chang^1,2,*^, Brent T. Schlegel^1,2,*^, Neil Carleton^1,2^,
Priscilla F. McAuliffe^1,3^, Hunter Waltermire^1,2^, Steffi Oesterreich^1,4^,
Russell Schwartz^5,6^, and Adrian V. Lee^1,4,7,#^
```

Changes v4 -> v5:
1. **Hunter Waltermire added** (new author, position 5)
2. Chang and Schlegel now co-first authors (both marked with *)
3. Lee marked with # instead of * (corresponding)
4. Affiliation numbering completely restructured (7 affiliations in v5 vs 3 in v4)
5. McAuliffe spelling corrected from "McAulife" to "McAuliffe" (one 'f' vs two)

- **WARNING (G1a)**: The spelling of McAuliffe is "McAuliffe" (two f's) in v5 L3, but v4 L3 has "McAulife" (one f). The CLAUDE.md project instructions also have the old one-f spelling: "McAulife". Cross-reference with clinical trial records indicates "McAuliffe" is correct (two f's). v5 is correct.

- **INFO (G1b)**: Hunter Waltermire is new in v5. The Author Contributions (L161) credits H.W. with "completed the wet lab validation of MDK signaling in ESR1-mutant cell lines," which is consistent with the new midkine section added in v5.

### G2. Affiliation Numbering (v5)

| # | v5 Affiliation | Authors |
|---|----------------|---------|
| 1 | Women's Cancer Research Center, UPMC Hillman Cancer Center, Magee-Womens Research Institute | Chang, Schlegel, Carleton, McAuliffe, Waltermire, Oesterreich, Lee |
| 2 | University of Pittsburgh, School of Medicine | Chang, Schlegel, Carleton, Waltermire |
| 3 | Department of Surgery, Division of Breast Surgical Oncology, UPitt SOM | McAuliffe |
| 4 | Department of Pharmacology and Chemical Biology, UPitt | Oesterreich, Lee |
| 5 | Ray and Stephanie Lane Computational Biology Dept, CMU | Schwartz |
| 6 | Department of Biological Sciences, CMU | Schwartz |
| 7 | Institute of Precision Medicine | Lee |

- **WARNING (G2a)**: In v4, Schlegel had only affiliation ^1^. In v5, Schlegel has ^1,2^. Carleton went from ^1^ to ^1,2^. These are plausible updates (e.g., both are in the School of Medicine), but worth confirming they are intentional.

- **INFO (G2b)**: In v4, Schwartz had affiliations ^2,3^ (Dept of Computational and Systems Biology, UPitt + Computational Biology Dept, CMU). In v5, Schwartz has ^5,6^ (Lane Computational Biology Dept, CMU + Biological Sciences, CMU). The UPitt affiliation for Schwartz was dropped. This may reflect an actual change in appointment or a correction.

- **INFO (G2c)**: Oesterreich's affiliation changed from ^1,2^ (v4: WCRC + Comp Sys Bio) to ^1,4^ (v5: WCRC + Pharmacology). This is a significant affiliation change that should be verified as intentional.

- **INFO (G2d)**: Lee's affiliations changed from ^1,2^ (v4) to ^1,4,7^ (v5), adding Institute of Precision Medicine. This appears intentional.

**No numbering errors detected** -- All superscript numbers map to defined affiliations, and all authors have at least one affiliation.

---

## H. Reference Numbering

### H1. Sequential Check

References are numbered 1-27 in v5. Let me verify sequentiality:
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27.

**PASS** -- All 27 references present, no gaps, no duplicates.

### H2. Citation-Reference Matching

| Citation in Text | Expected Reference | Actual Reference (v5) | Status |
|------------------|--------------------|----------------------|--------|
| [8] = Cell2Location | Kleshchevnikov et al. 2022 | Ref 8: Kleshchevnikov et al. Cell2location. Nat Biotechnol | OK |
| [16] = Tangram | Biancalani et al. 2021 | Ref 16: Biancalani et al. Tangram. Nat Methods | OK |
| [10] = Seurat | Butler/Satija et al. | Ref 10: Butler et al. Integrating single-cell. Nat Biotechnol | OK |
| [9] = PyDESeq2/DESeq2 | Love et al. 2014 | Ref 9: Love et al. DESeq2. Genome Biol | OK |
| [11] = Wu atlas | Wu et al. 2021 | Ref 11: Wu et al. breast cancers atlas. Nat Genet | OK |
| [12] = scCube | Qian et al. 2024 | Ref 12: Qian et al. scCube. Nat Commun | OK |
| [13] = ESR1 mutant | Li et al. 2022 | Ref 13: Li et al. ESR1 mutant. Nat Commun | OK |
| [14] = Trial | McAuliffe 2025 | Ref 14: McAuliffe. clinicaltrials.gov NCT05914792 | OK |
| [15] = Trial medRxiv | Carleton et al. 2025 | Ref 15: Carleton et al. medRxiv | OK |
| [17] = spatial overview | Li & Wang 2021 | Ref 17: Li, Wang. bulk to spatial RNA seq | OK |
| [18] = benchmarking | Li et al. 2023 | Ref 18: Li et al. benchmarking deconvolution. Nat Commun | OK |
| [20] = EstroGene | Li et al. 2024 | Ref 20: Li et al. EstroGene2.0. bioRxiv | OK |
| [21] = Midkine/cancer | Visvader 2024 | Ref 21: Visvader. Midkine. Cancer Cell | OK |
| [22] = Midkine/aging | Yan et al. 2024 | Ref 22: Yan et al. Midkine driver. Cancer Cell | OK |
| [23] = Midkine biology | Muramatsu 2010 | Ref 23: Muramatsu. Midkine. Proc Jpn Acad | OK |
| [27] = COMMOT | Cang et al. 2023 | Ref 27: Cang et al. COMMOT. Nat Methods | OK |
| [24,25] = ESR1 context | Arnesen + Li 2022 | Refs 24,25: Arnesen et al., Li et al. Cancer Res | OK |

**PASS** -- All citations match their references.

### H3. References Without Text Citations

Checking whether any references are cited in text:

- Ref 1 (Huang et al.): Cited as ^1^ in L47 ("approximately 85%")
- Ref 2 (Patel et al.): Cited as ^2^ in L47 ("first-line treatment")
- Ref 3 (Regan & Preall): Cited as ^3^ in L127
- Ref 4 (Chang et al.): Cited as ^4^ in L127
- Ref 5 (Wess et al.): Cited as ^5^ in L49, L53
- Ref 6 (Zaidi et al.): Cited as ^6^ in L49
- Ref 7 (Chen et al.): Cited as ^7^ in L49, L125
- Ref 19 (Anders & Huber): Not cited in v5 text

- **WARNING (H3a)**: Reference 19 (Anders & Huber, 2010) is listed in the reference list but does not appear to be cited anywhere in the v5 text. In v4 it was cited as [9,19] in the Methods (L223). In v5, the Methods PyDESeq2 citation (L234) only uses ^9^. Reference 19 is an orphaned reference.

  **Severity: WARNING** -- Ref 19 is uncited. Either remove it from the reference list or add a citation (e.g., in the PyDESeq2 methods paragraph alongside ref 9, since DESeq2 builds on the Anders & Huber method).

- Ref 26 (Guo et al.): Cited as ^26^ in L49

**All other references are cited at least once.**

---

## I. Additional Cross-Document Consistency Issues

### I1. v4 vs v5 Title Change

- **v4**: "CITEgeist: A Spatially-Native Framework for Multi-Modal Integration and Program Discovery in Spatial Transcriptomics"
- **v5**: "CITEgeist: Accurate deconvolution of spatial transcriptomics with same-slide proteomics reveals midkine as a secreted microenvironment modulator in ESR1 mutant breast cancer"

**INFO** -- Title changed to reflect the biological finding (midkine). This is a significant revision but intentional.

### I2. v4 "immunotherapy" vs v5 "endocrine therapy"

- **v4 (L37)**: "identifying transcriptional programs that distinguish immunotherapy responders from progressors"
- **v5 (L107)**: "primary endocrine therapy (aromatase inhibitors)"

- This was a **CRITICAL error in v4** that has been correctly fixed in v5. The trial is for endocrine therapy, not immunotherapy. No residual issue in v5.

### I3. v4 "14-patient" vs v5 "14 specimens from 6 patients"

- **v4 Abstract (L15)**: "a 14-patient breast cancer cohort"
- **v5 Summary (L20)**: "14 specimens" (later clarified as "6 patients")
- **v5 Section 2.5 (L107)**: "14 specimens from 6 patients"

This was a **CRITICAL error in v4** (implying 14 patients when it was 6). Corrected in v5.

### I4. v4 Figure 6 Descriptions vs v5 Content

v4 had a separate Section 2.6 on interoperability with Figure 6 showing PyDESeq2, GSEApy, scanpy, COMMOT integration. v5 consolidates this into Section 2.5 (L117) and the Discussion. The Figure 6 content (midkine discovery-to-validation, interoperability demos) may need separate figure legends in v5's figure document, which is not yet created for v5.

**INFO** -- v5 appears to restructure figures (v4 had 6 figures; v5 appears to have 5 main figures based on references). The figure legends document has not been updated for v5.

### I5. Tangram RMSE Value

- **v5 Results (L85)**: Does not report Tangram's RMSE explicitly in the Results prose
- **v4 (L83)**: Reports Tangram RMSE = 8.33 (also in verified_numbers.md: 8.331 +/- 0.446)
- **v5 Figure 3B legend (from v4 Figs)**: Lists CITEgeist r = 0.60, but does not list Tangram RMSE

**INFO** -- Tangram RMSE omitted from v5 prose, which is fine since the r value is reported.

### I6. "Neoadjuvant therapy" terminology

- **v4 Section 2.5 (L107)**: "6 patients undergoing neoadjuvant therapy"
- **v5**: Uses "primary endocrine therapy" throughout

- **INFO**: "Primary endocrine therapy" is the correct clinical term for this trial (patients chose endocrine therapy *instead of* upfront surgery, which is distinct from neoadjuvant therapy given *before* planned surgery). v5 is correct; v4 was wrong.

---

## J. Figures Cross-Reference

The v5 manuscript references the following figures:

| Reference | v5 Location | Figure Content |
|-----------|-------------|----------------|
| Figure 1A | L63 | Pipeline overview |
| Figure 1B-C | L67 | Outputs/interoperability |
| Figure 2A-B | L72 | Profile discovery |
| Figure 2C | L77 | Xenium validation |
| Figure 2D | L77 | Marker clustering validation |
| Figure 3A | L81 | Deconvolution schematic |
| Figure 3B | L83 | Benchmarking (simulated) |
| Figure 4A | L93 | Basal signatures |
| Figure 4B | L93 | Spatial distribution |
| Figure 4C | L93 | EstroGene signatures |
| Figure 4D-E | L93 | Pathway alterations |
| Figure 4F | L95 | COMMOT signaling |
| Figure 4G | L97 | ELISA results |
| Figure 4H | L97 | IF results |
| Figure 5A | L109 | Module 4 programs |
| Figure 5B | L111 | Integration / 73 aligned programs |
| Figure 5C | L113 | Response programs |
| Figure 5D | L113 | DE volcano plot |
| Supplemental Figure S3 | L101 | T47D results |

- **WARNING (J1)**: v5 references Figures 1-5 with subpanels. v4 had Figures 1-6. The v4 figure legends document (`CITEgeist_Patterns_v4_Figures.md`) has legends for Figures 1-6. v5 appears to have consolidated: Figure 4 in v5 is the midkine case study (entirely new vs v4's Module 4 programs), and Figure 5 in v5 is the cohort/integration analysis (was Figures 4+5 in v4). A new v5 figures document with updated legends has not been created yet.

  **Severity: WARNING** -- The v4 figure legends document does not match v5 figure numbering. A v5 figure legends document should be created.

---

## Summary

### CRITICAL Issues (1)

1. **C2a**: Module 4 Moran's I threshold mismatch between Methods (I > 0.15) and Results (I > 0.2). The "68% exceeded threshold" claim cannot be evaluated without knowing which threshold was actually used.

### WARNING Issues (6)

1. **B2a**: Ambiguity about 6 vs 7 cell types between verified_numbers.md and manuscript
2. **B6a**: Methods says 12,000 cells for simulation; Results discusses 30,000 and 8,000 for reference sizes. Could confuse readers.
3. **G1a**: McAuliffe spelling corrected v4->v5 (good), but CLAUDE.md still has old spelling
4. **G2a**: Multiple authors gained new affiliations between v4 and v5 (verify intentional)
5. **H3a**: Reference 19 (Anders & Huber 2010) is uncited in v5 text
6. **J1**: v4 figure legends document does not match v5 figure numbering; v5 figure legends not yet created

### INFO Issues (15)

1. A2a: Summary omits age/HER2- specifics (normal for summary)
2. A2b: "neoadjuvant" -> "primary endocrine therapy" correction (improvement)
3. B3a: JSD values dropped from v5 Results prose (acceptable)
4. B5a: New simulation numbers in v5 not in v4 (new content)
5. C1a: 175 program count from v4 not explicitly in v5 (different context)
6. C1b: 57% (v4 Xenium) vs 68% (v5 cohort) are different datasets (OK but could confuse)
7. C6a: One responder program gene list omitted from v5 prose (acceptable with "including")
8. E1a: 122 + 5 = 127 (math checks out)
9. F4a: Methods more specific than Results on edge retention criteria (standard)
10. F5a: Statistical test names in STAR Methods but not Results prose (standard)
11. G1b: Hunter Waltermire addition justified by wet lab contribution
12. G2b-d: Affiliation changes for Schwartz, Oesterreich, Lee (verify intentional)
13. I1: Title change reflects biological finding addition
14. I4: Figure reorganization between v4 and v5
15. I5: Tangram RMSE omitted from v5 prose (acceptable)

---

## Recommended Actions

### Must Fix Before Submission
1. **Resolve Module 4 Moran's I threshold** (C2a): Decide whether I > 0.15 or I > 0.2 is the correct threshold, update both Methods and Results to match, and verify the "68%" figure against the chosen threshold.

### Should Fix
2. **Add citation for Reference 19** (H3a) or remove it from the reference list.
3. **Clarify simulation cell counts** (B6a): In Methods, distinguish between scCube simulation input (12,000 cells) and reference atlas sizes for competing methods (30,000 vs 8,000).
4. **Create v5 figure legends document** (J1): Update figure legends to match v5's reorganized figure numbering (Figures 1-5 + Supplementals).
5. **Verify affiliation changes** (G2a-d): Confirm that all affiliation number updates for Schlegel, Carleton, Schwartz, Oesterreich, and Lee are intentional.

### Minor / Can Wait
6. Update CLAUDE.md spelling of "McAulife" -> "McAuliffe" (G1a).
7. Clarify 6 vs 7 cell types context (B2a) in verified_numbers.md or manuscript.

---

*Audit performed on 2026-02-08 by automated review.*
*Files compared: CITEgeist_CRM_v5.md, CITEgeist_Patterns_v4.md, CITEgeist_Patterns_v4_Figures.md, verified_numbers.md*
