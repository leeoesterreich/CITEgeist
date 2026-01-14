# CITEgeist Pipeline Diagram

This document provides diagrams for the CITEgeist pipeline (Modules 1-5) intended as a rough draft for journal figures. Final versions should be created in BioRender or similar tools.

---

## Part 1: Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           RAW SPATIAL DATA                                  │
│         Antibody Capture (CITE-seq) + Gene Expression (Visium)             │
│                     N spots × M markers, N spots × G genes                  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 1: MARKER INTEREST DETECTION                                       │
│  ═══════════════════════════════════                                        │
│                                                                             │
│     ┌───────────┐     ┌───────────┐     ┌───────────┐                      │
│     │ Kurtosis  │     │    GMM    │     │ Moran's I │                      │
│     │  Gating   │     │    SNR    │     │  Spatial  │                      │
│     │  (>2.0)   │     │  (≥1.0)   │     │Autocorrel │                      │
│     └─────┬─────┘     └─────┬─────┘     └─────┬─────┘                      │
│           │                 │                 │                             │
│           └────────┐   ┌────┴────┐   ┌────────┘                             │
│                    ▼   ▼         ▼   ▼                                      │
│              (Kurtosis OR Moran's I) AND GMM                                │
│                            │                                                │
│                            ▼                                                │
│              interest_score = kurtosis × SNR × max(I, 0)                    │
│                                                                             │
│  Output: MarkerInterestResult (ranked interesting markers)                  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼ Interesting Markers (M')
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 2: PROFILE ASSEMBLY                                                 │
│  ═════════════════════════                                                  │
│                                                                             │
│  ┌─────────────────────┐   ┌─────────────────────┐   ┌─────────────────┐   │
│  │  2a: Colocalization │ → │  2b: Profile        │ → │  2c: Profile    │   │
│  │      Analysis       │   │      Discovery      │   │     Selection   │   │
│  │                     │   │                     │   │                 │   │
│  │ • Jaccard index     │   │ • FDR correction    │   │ • Greedy by     │   │
│  │ • Pearson/Spearman  │   │ • Top-k filtering   │   │   variance      │   │
│  │ • Cosine similarity │   │ • Graph components  │   │ • Dual target:  │   │
│  │ • Bivariate Moran's │   │ • Hierarchical      │   │   90% protein   │   │
│  │   I (40% weight)    │   │   clustering        │   │   90% spatial   │   │
│  └─────────────────────┘   └─────────────────────┘   └─────────────────┘   │
│                                                                             │
│  Output: ProfileSelectionResult (5-15 cell type marker profiles)           │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼ Cell Type Profiles
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 3: DECONVOLUTION                                                    │
│  ═══════════════════════                                                    │
│                                                                             │
│  ┌───────────────────────────────┐   ┌───────────────────────────────┐     │
│  │  Pass 1: Cell Proportions    │   │  Pass 2: Gene Expression      │     │
│  │  ─────────────────────────   │   │  ───────────────────────────  │     │
│  │                               │   │                               │     │
│  │  Input: Antibody data (N×M') │ → │  Input: GEX (N×G) + Props     │     │
│  │                               │   │                               │     │
│  │  Method: Gurobi QP            │   │  Method: Per-spot MIP         │     │
│  │  • Per-marker β optimization  │   │  • Global+local enrichment    │     │
│  │  • Laplacian smoothing        │   │  • Count-preserving           │     │
│  │  • Local finetuning           │   │  • Neighbor context           │     │
│  │                               │   │                               │     │
│  │  Output: Y (N×T proportions)  │   │  Output: Layers (N×G per T)   │     │
│  └───────────────────────────────┘   └───────────────────────────────┘     │
│                                                                             │
│  Output: Cell proportions + Deconvolved expression layers per cell type    │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼ Deconvolved Layers
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 4: ANCHORED PROGRAM DISCOVERY                                       │
│  ════════════════════════════════════                                       │
│                                                                             │
│  ┌───────────────────────────────┐   ┌───────────────────────────────┐     │
│  │  4a: Per-Cell-Type Programs  │   │  4b: Bivariate Relationships  │     │
│  │  ─────────────────────────   │   │  ─────────────────────────────│     │
│  │                               │   │                               │     │
│  │  For each cell type anchor:   │ → │  For all program pairs:       │     │
│  │  • Extract deconvolved layer  │   │  • Bivariate Moran's I        │     │
│  │  • Weighted NMF (K programs)  │   │  • Permutation test           │     │
│  │  • Spatial smoothing          │   │  • FDR correction             │     │
│  │  • Moran's I per program      │   │  • Classify: co-localized,    │     │
│  │  • Subpopulation detection    │   │    exclusive, independent     │     │
│  │                               │   │                               │     │
│  │  Output: W (G×K), H (K×N)     │   │  Output: Program network      │     │
│  └───────────────────────────────┘   └───────────────────────────────┘     │
│                                                                             │
│  Output: Spatial programs per cell type + within-sample relationships      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼ Programs + Relationships (per sample)
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 5: CROSS-SAMPLE INTEGRATION                                         │
│  ══════════════════════════════════                                         │
│                                                                             │
│  ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐   ┌──────────┐  │
│  │   Gene   │ → │  Harmony │ → │ Program  │ → │Conserved │ → │ Network  │  │
│  │Alignment │   │Correction│   │ Matching │   │Relations │   │ Building │  │
│  │          │   │          │   │          │   │          │   │          │  │
│  │ Union of │   │ PCA +    │   │ Cosine   │   │ Compare  │   │ Nodes:   │  │
│  │ all genes│   │ k-means +│   │similarity│   │bivariate │   │ programs │  │
│  │ Zero-pad │   │ iterative│   │clustering│   │ I across │   │ Edges:   │  │
│  │ W matrices│   │ correct  │   │          │   │ samples  │   │ conserved│  │
│  └──────────┘   └──────────┘   └──────────┘   └──────────┘   └──────────┘  │
│                                                                             │
│  Output: IntegrationResult (aligned programs + conserved relationships)    │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                            FINAL OUTPUT                                     │
│  ═══════════════════════════════════════                                    │
│                                                                             │
│  • Conserved gene expression programs across patients                       │
│  • Consensus gene signatures per program                                    │
│  • Spatial relationship network (co-localized vs exclusive programs)       │
│  • Conservation scores for cross-patient reproducibility                    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Part 2: Detailed Module Breakdowns

### Module 1: Marker Interest Detection

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  INPUT: Raw Antibody Expression Matrix X (N_spots × M_markers)              │
│         Spatial Coordinates (N_spots × 2)                                   │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
           ┌──────────────────────────┼──────────────────────────┐
           ▼                          ▼                          ▼
┌─────────────────────┐   ┌─────────────────────┐   ┌─────────────────────┐
│  TEST 1: KURTOSIS   │   │  TEST 2: GMM SNR    │   │  TEST 3: MORAN'S I  │
│  ─────────────────  │   │  ────────────────   │   │  ─────────────────  │
│                     │   │                     │   │                     │
│  Per marker:        │   │  Per marker:        │   │  Per marker:        │
│  κ = E[(X-μ)⁴]/σ⁴-3│   │  Fit 2-component    │   │  1. Spatial smooth  │
│                     │   │  GMM: N(μ_bg,σ_bg)  │   │     (k=6 neighbors) │
│  High κ (>2.0):     │   │      + N(μ_sig,σ_sig)│   │  2. Z-score         │
│  Peaked → Signal    │   │                     │   │  3. Compute I on    │
│                     │   │  SNR = (μ_sig-μ_bg) │   │     k=8 neighbors   │
│  Low κ (<2.0):      │   │        / σ_bg       │   │                     │
│  Flat → Noise       │   │                     │   │  I > threshold:     │
│                     │   │  SNR ≥ 1.0: Pass    │   │  Spatially clustered│
│                     │   │  SNR < 1.0: Fail    │   │                     │
└─────────┬───────────┘   └─────────┬───────────┘   └─────────┬───────────┘
          │                         │                         │
          │        Adaptive GMM on metric values              │
          │        learns thresholds from data                │
          │                         │                         │
          ▼                         ▼                         ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                           DECISION LOGIC                                    │
│  ────────────────────────────────────────────────────────────────────────   │
│                                                                             │
│                   ┌────────────┐                                            │
│                   │  Kurtosis  │                                            │
│                   │   Pass?    │                                            │
│                   └──────┬─────┘                                            │
│                          │                                                  │
│           ┌──────────────┼──────────────┐                                   │
│           │              │              │                                   │
│           ▼              ▼              ▼                                   │
│       ┌───────┐     ┌─────────┐    ┌───────┐                               │
│       │  Yes  │     │   OR    │    │  No   │                               │
│       └───┬───┘     └────┬────┘    └───┬───┘                               │
│           │              │             │                                    │
│           │              │             ▼                                    │
│           │              │      ┌────────────┐                              │
│           │              │      │ Moran's I  │                              │
│           │              │      │   Pass?    │                              │
│           │              │      └──────┬─────┘                              │
│           │              │             │                                    │
│           │         ┌────┴────┐        │                                    │
│           ▼         ▼         ▼        ▼                                    │
│       ┌─────────────────────────────────────┐                               │
│       │           AND GMM Pass?             │                               │
│       └──────────────────┬──────────────────┘                               │
│                          │                                                  │
│           ┌──────────────┴──────────────┐                                   │
│           ▼                             ▼                                   │
│    ┌─────────────┐               ┌─────────────┐                            │
│    │ INTERESTING │               │    BORING   │                            │
│    │   MARKER    │               │    MARKER   │                            │
│    └─────────────┘               └─────────────┘                            │
│                                                                             │
│  Combined Score: interest_score = kurtosis × GMM_SNR × max(Moran's_I, 0)   │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: MarkerInterestResult                                               │
│  • interesting_markers: List of markers passing gates                       │
│  • boring_markers: List of markers failing gates                            │
│  • Ranked by interest_score (descending)                                    │
│  • Per-marker: kurtosis, gmm_snr, morans_i, passed_* flags                 │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Module 2a: Marker Colocalization Analysis

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  INPUT: Interesting Markers from Module 1                                   │
│         Expression Matrix X (N_spots × M' interesting markers)              │
│         Spatial Coordinates                                                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
           ┌──────────────────────────┼──────────────────────────┐
           │                          │                          │
           ▼                          ▼                          ▼
┌─────────────────────┐   ┌─────────────────────┐   ┌─────────────────────┐
│  PREPROCESSING      │   │  For each marker    │   │  Build neighbor     │
│  ─────────────────  │   │  pair (i,j):        │   │  graph (k=6)        │
│                     │   │                     │   │                     │
│  1. k-NN graph      │   │  M'×(M'-1)/2 pairs  │   │  Spatial smoothing  │
│  2. Spatial smooth  │   │  ~parallel compute  │   │  per marker         │
│  3. Binarize @75th  │   │                     │   │                     │
└─────────────────────┘   └─────────────────────┘   └─────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│              FOUR COLOCALIZATION METRICS PER PAIR                           │
│  ────────────────────────────────────────────────────────────────────────   │
│                                                                             │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐             │
│  │ 1. CO-OCCURRENCE│  │ 2. CORRELATION  │  │ 3. SIMILARITY   │             │
│  │ (Binary)        │  │ (Continuous)    │  │ (Continuous)    │             │
│  │                 │  │                 │  │                 │             │
│  │ Jaccard index:  │  │ Pearson r:      │  │ Cosine:         │             │
│  │ J = |A∩B|/|A∪B| │  │ cov(A,B)/σ_Aσ_B│  │ A·B/(||A||||B||)│             │
│  │                 │  │                 │  │                 │             │
│  │ On binarized    │  │ Spearman ρ:     │  │ Scale-invariant │             │
│  │ expression      │  │ rank correlation│  │ pattern match   │             │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘             │
│                                                                             │
│  ┌─────────────────────────────────────────────────────────────┐           │
│  │ 4. BIVARIATE MORAN'S I (40% weight - KEY METRIC)           │           │
│  │                                                             │           │
│  │ I_AB = Σ_i Σ_j w_ij × (A_i - μ_A) × (B_j - μ_B)            │           │
│  │        ────────────────────────────────────────             │           │
│  │                    n × σ_A × σ_B                            │           │
│  │                                                             │           │
│  │ Measures: Does high A co-occur with neighboring high B?     │           │
│  │ Range: [-1, +1]  (+: co-localized, -: exclusive)           │           │
│  │                                                             │           │
│  │ Computed on spatially-smoothed data                         │           │
│  │ Permutation test (199 perms) for p-values                   │           │
│  │ Optional multi-scale: k ∈ [6, 12, 24, 48, 64]              │           │
│  └─────────────────────────────────────────────────────────────┘           │
│                                                                             │
│  COMBINED SCORE:                                                            │
│  colocalization_score = 0.3×Spearman + 0.3×Cosine + 0.4×Bivariate_I       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: ColocalizationResult                                               │
│  • pairs: List of MarkerPairColocalization (all metrics per pair)          │
│  • marker_specificity: Gini coefficient per marker (high = cell-specific)  │
│  • Ranked by colocalization_score                                           │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Module 2b: Profile Discovery

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  INPUT: ColocalizationResult from Module 2a                                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 1: FDR CORRECTION                                                     │
│  ────────────────────                                                       │
│                                                                             │
│  • Apply Benjamini-Hochberg on bivariate_morans p-values                   │
│  • FDR threshold α = 0.05                                                   │
│  • Keep edges with q-value < α                                              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 2: MUTUAL TOP-K SPARSIFICATION                                        │
│  ──────────────────────────────────                                         │
│                                                                             │
│         Marker A's top-k partners          Marker B's top-k partners       │
│         ┌───────────────────┐              ┌───────────────────┐           │
│         │ 1. B  (0.85)      │              │ 1. A  (0.85)      │           │
│         │ 2. C  (0.72)      │              │ 2. D  (0.68)      │           │
│         │ 3. D  (0.65)      │              │ 3. E  (0.55)      │           │
│         │ 4. E  (0.58)      │              │ 4. C  (0.52)      │           │
│         │ 5. F  (0.51)      │              │ 5. F  (0.48)      │           │
│         └───────────────────┘              └───────────────────┘           │
│                  │                                   │                      │
│                  └───────────┬───────────────────────┘                      │
│                              ▼                                              │
│                    A—B edge KEPT (both rank each other in top-k)           │
│                    A—C edge REMOVED (C doesn't rank A in top-k)            │
│                                                                             │
│  Prevents "hub collapse" where promiscuous markers connect everything      │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 3: GRAPH CONSTRUCTION & CONNECTED COMPONENTS                          │
│  ─────────────────────────────────────────────────                          │
│                                                                             │
│         Filtered edges form graph G                                         │
│                                                                             │
│              [A]───[B]           [E]───[F]                                  │
│               │     │                  │                                    │
│              [C]───[D]           [G]───[H]                                  │
│                                                                             │
│           Component 1          Component 2         [I]  [J]                 │
│           (4 markers)          (4 markers)       Singletons                 │
│                                                                             │
│  Each connected component = potential cell type lineage                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 4: HIERARCHICAL CLUSTERING + DYNAMIC TREE CUTTING                     │
│  ─────────────────────────────────────────────────────                      │
│                                                                             │
│  For each component:                                                        │
│                                                                             │
│  1. Build distance matrix: d_ij = 1 - colocalization_score_ij              │
│                                                                             │
│  2. Hierarchical clustering (average linkage / UPGMA)                       │
│                                                                             │
│  3. Dynamic tree cutting (gap-based):                                       │
│                                                                             │
│              Height                                                         │
│                │                                                            │
│           0.8 ─┼───────────────────────────┐                               │
│                │                           │ ← Large gap = lineage split   │
│           0.5 ─┼───────────┬───────────────┘                               │
│                │           │                                                │
│           0.3 ─┼───┬───────┤       ┌───────┐                               │
│                │   │       │       │       │                                │
│           0.0 ─┼───┴───────┴───────┴───────┴──                              │
│                │   A   B   C   D       E   F                                │
│                │                                                            │
│                └─────────────────────────────▶ Markers                      │
│                                                                             │
│           Cut at gap → Profile 1: [A,B,C,D], Profile 2: [E,F]              │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: ProfileDiscoveryResult                                             │
│  • profiles: List of discovered cell type marker profiles                   │
│  • lineage_dendrograms: Dendrogram per connected component                  │
│  • singletons: Markers forming their own profile                            │
│  • modularity: Quality metric (how well profiles explain structure)         │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Module 2c: Profile Selection

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  INPUT: Candidate profiles from Module 2b (may be 10-30 profiles)          │
│         Original marker expression data                                     │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  DUAL VARIANCE FRAMEWORK                                                    │
│  ═══════════════════════                                                    │
│                                                                             │
│  ┌─────────────────────────────────┐  ┌─────────────────────────────────┐  │
│  │  V_protein (Signal Variance)   │  │  V_spatial (Spatial Variance)   │  │
│  │  ─────────────────────────────  │  │  ───────────────────────────── │  │
│  │                                 │  │                                 │  │
│  │  Σ var(X[selected_markers])     │  │  Σ var(marker) × |Moran's I|   │  │
│  │  ─────────────────────────────  │  │  ──────────────────────────    │  │
│  │     Σ var(X[all_markers])       │  │      total_spatial_variance    │  │
│  │                                 │  │                                 │  │
│  │  Goal: Cover 90% of protein    │  │  Goal: Cover 90% of spatially  │  │
│  │        signal variance          │  │        structured variance     │  │
│  └─────────────────────────────────┘  └─────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  GREEDY SELECTION ALGORITHM                                                 │
│  ═════════════════════════                                                  │
│                                                                             │
│  selected = []                                                              │
│  while not (V_spatial ≥ 0.9 AND V_protein ≥ 0.9):                          │
│      for each remaining profile:                                            │
│          score = marginal_spatial_variance_gain(profile)                    │
│      best = argmax(score)                                                   │
│      selected.append(best)                                                  │
│      update V_spatial, V_protein                                            │
│                                                                             │
│  ┌───────────────────────────────────────────────────────────────────────┐ │
│  │  Step  │ Profile Added │ V_spatial │ V_protein │ Marginal Gain       │ │
│  │────────┼───────────────┼───────────┼───────────┼─────────────────────│ │
│  │   1    │ T-cell        │   0.25    │   0.20    │     0.25            │ │
│  │   2    │ Macrophage    │   0.48    │   0.42    │     0.23            │ │
│  │   3    │ B-cell        │   0.65    │   0.58    │     0.17            │ │
│  │   4    │ Epithelial    │   0.78    │   0.72    │     0.13            │ │
│  │   5    │ Fibroblast    │   0.88    │   0.85    │     0.10            │ │
│  │   6    │ Endothelial   │   0.92    │   0.91    │     0.04  ← STOP    │ │
│  └───────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: ProfileSelectionResult                                             │
│  • selected_profiles: Ordered list of 5-15 profiles                        │
│  • optimal_n: Number of profiles selected                                   │
│  • variance_explained: Cumulative variance at each step                     │
│  • stopping_reason: "targets_reached" / "max_profiles" / "no_remaining"    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Module 3: Deconvolution (Two-Pass Algorithm)

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 3 OVERVIEW                                                          │
│  ════════════════                                                           │
│                                                                             │
│                    ┌───────────────────────────────┐                        │
│                    │  Antibody Data (N × M)        │                        │
│                    │  Gene Expression (N × G)      │                        │
│                    │  Cell Type Profiles (T types) │                        │
│                    └───────────────┬───────────────┘                        │
│                                    │                                        │
│                                    ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  PASS 1: CELL PROPORTION ESTIMATION                                 │   │
│  │                                                                     │   │
│  │  Input:  Antibody S (N×M), Profile assignment (M→T)                │   │
│  │  Output: Proportions Y (N×T), Betas β (M scaling factors)          │   │
│  │                                                                     │   │
│  │  min Σ_i,m (S[i,m] - β[m] × Y[i, owner(m)])²                       │   │
│  │      + λ_reg × ||Y||²                                               │   │
│  │      + λ_lap × Σ_ij w_ij × ||Y[i] - Y[j]||²  (Laplacian smoothing) │   │
│  │                                                                     │   │
│  │  s.t.  0.9 ≤ Σ_j Y[i,j] ≤ 1.2  (soft sum-to-one)                   │   │
│  │        Y[i,j] ≥ 0                                                   │   │
│  │        β[m] ≥ 0                                                     │   │
│  │                                                                     │   │
│  │  Solver: Gurobi Quadratic Programming                               │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                    │                                        │
│                                    ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  LOCAL FINETUNING                                                   │   │
│  │                                                                     │   │
│  │  For each spot i:                                                   │   │
│  │    neighbors = get_neighbors(i, radius)                             │   │
│  │    local_optimize(S[neighbors], Y_global[neighbors])                │   │
│  │    Y_fine[i] = refined proportion                                   │   │
│  │                                                                     │   │
│  │  Captures local heterogeneity missed by global optimization         │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                    │                                        │
│                                    ▼                                        │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │  PASS 2: GENE EXPRESSION DECONVOLUTION                              │   │
│  │                                                                     │   │
│  │  Input:  GEX counts C (N×G), Proportions Y (N×T)                   │   │
│  │  Output: Deconvolved X (per cell type: N×G each)                   │   │
│  │                                                                     │   │
│  │  For each spot i:                                                   │   │
│  │    Variables: X[j,k] = gene k counts for cell type j               │   │
│  │                                                                     │   │
│  │    maximize: Σ_j,k enrichment[k,j] × Y[i,j] × X[j,k]               │   │
│  │              - λ_prior × deviation_from_prior(X)                    │   │
│  │                                                                     │   │
│  │    s.t.  Σ_j X[j,k] = C[i,k]  (counts sum to observed)             │   │
│  │          X[j,k] ≥ 0, integer                                        │   │
│  │                                                                     │   │
│  │  Solver: Gurobi Mixed Integer Programming (per spot)                │   │
│  │                                                                     │   │
│  │  enrichment = global_enrichment × (1-α) + local_enrichment × α     │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                    │                                        │
│                                    ▼                                        │
│                    ┌───────────────────────────────┐                        │
│                    │  OUTPUT:                      │                        │
│                    │  • Y: Proportions (N × T)     │                        │
│                    │  • Layers: T-cell_genes,      │                        │
│                    │    Macrophage_genes, etc.     │                        │
│                    │    (each N × G)               │                        │
│                    └───────────────────────────────┘                        │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Detailed Pass 1 (Proportions):**

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  PASS 1: CELL PROPORTION OPTIMIZATION                                       │
│  ════════════════════════════════════                                        │
│                                                                             │
│  Antibody Matrix S (N×M)              Profile Assignment (M→T)              │
│  ┌────────────────────┐               ┌────────────────────┐                │
│  │ Spot│CD3 │CD68│CD19│               │ CD3  → T-cell      │                │
│  │─────┼────┼────┼────│               │ CD4  → T-cell      │                │
│  │  1  │ 8.2│ 1.1│ 0.5│               │ CD68 → Macrophage  │                │
│  │  2  │ 2.1│ 6.5│ 0.8│               │ CD163→ Macrophage  │                │
│  │  3  │ 0.3│ 0.5│ 7.2│               │ CD19 → B-cell      │                │
│  │ ... │ ...│ ...│ ...│               │ ...                │                │
│  └────────────────────┘               └────────────────────┘                │
│           │                                    │                            │
│           └────────────────┬───────────────────┘                            │
│                            ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  Gurobi QP Problem Setup                                             │  │
│  │                                                                      │  │
│  │  Variables:                                                          │  │
│  │    Y[i,j] ∈ [0,1]: proportion of cell type j at spot i              │  │
│  │    β[m] ≥ 0: scaling factor for marker m                            │  │
│  │                                                                      │  │
│  │  Objective (minimize):                                               │  │
│  │    Σ_i Σ_m (S[i,m] - β[m] × Y[i, owner(m)])²   ← reconstruction    │  │
│  │    + λ_reg × Σ_i ||Y[i]||²                      ← regularization    │  │
│  │    + λ_lap × Σ_ij w_ij × ||Y[i] - Y[j]||²      ← spatial smooth    │  │
│  │                                                                      │  │
│  │  Constraints:                                                        │  │
│  │    0.9 ≤ Σ_j Y[i,j] ≤ 1.2  for all spots i (soft sum-to-one)       │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                            │                                                │
│                            ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  Output: Proportion Matrix Y (N × T)                                 │  │
│  │  ┌─────────────────────────────────┐                                │  │
│  │  │ Spot│T-cell│Macro │B-cell│Other│                                │  │
│  │  │─────┼──────┼──────┼──────┼─────│                                │  │
│  │  │  1  │ 0.72 │ 0.15 │ 0.08 │0.05 │  ← T-cell dominant            │  │
│  │  │  2  │ 0.18 │ 0.65 │ 0.12 │0.05 │  ← Macrophage dominant        │  │
│  │  │  3  │ 0.05 │ 0.08 │ 0.82 │0.05 │  ← B-cell dominant            │  │
│  │  │ ... │ ...  │ ...  │ ...  │ ... │                                │  │
│  │  └─────────────────────────────────┘                                │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Detailed Pass 2 (Gene Expression):**

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  PASS 2: GENE EXPRESSION DECONVOLUTION                                      │
│  ═════════════════════════════════════                                       │
│                                                                             │
│  Gene Expression C (N×G)              Proportions Y (N×T)                   │
│  ┌─────────────────────┐              ┌─────────────────────┐               │
│  │Spot│TRAC│CD68│MS4A1│              │Spot│T-cell│Macro│B-cell│              │
│  │────┼────┼────┼─────│              │────┼──────┼─────┼──────│              │
│  │ 1  │ 150│  20│   30│              │ 1  │ 0.72 │ 0.15│ 0.08 │              │
│  │ 2  │  40│ 180│   25│              │ 2  │ 0.18 │ 0.65│ 0.12 │              │
│  │ 3  │  10│  15│  175│              │ 3  │ 0.05 │ 0.08│ 0.82 │              │
│  └─────────────────────┘              └─────────────────────┘               │
│            │                                    │                            │
│            └────────────────┬───────────────────┘                            │
│                             ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐   │
│  │  For Spot 1: C = [TRAC:150, CD68:20, MS4A1:30]                       │   │
│  │              Y = [T-cell:0.72, Macro:0.15, B-cell:0.08]              │   │
│  │                                                                      │   │
│  │  Variables: X[j,k] = counts for cell type j, gene k                 │   │
│  │                                                                      │   │
│  │  Enrichment (gene k → cell type j preference):                      │   │
│  │  ┌─────────────────────────────────┐                                │   │
│  │  │ Gene │T-cell│Macro│B-cell│                                       │   │
│  │  │──────┼──────┼─────┼──────│                                       │   │
│  │  │ TRAC │ 0.95 │ 0.02│ 0.03 │  ← T-cell gene                       │   │
│  │  │ CD68 │ 0.05 │ 0.90│ 0.05 │  ← Macrophage gene                   │   │
│  │  │MS4A1 │ 0.02 │ 0.03│ 0.95 │  ← B-cell gene                       │   │
│  │  └─────────────────────────────────┘                                │   │
│  │                                                                      │   │
│  │  Objective: maximize Σ enrichment[k,j] × Y[j] × X[j,k]              │   │
│  │             (assign counts to most enriched cell type present)       │   │
│  │                                                                      │   │
│  │  Constraint: Σ_j X[j,k] = C[k]  (e.g., X[T,TRAC] + X[M,TRAC] + ... │   │
│  │              X[B,TRAC] = 150 for TRAC)                              │   │
│  └──────────────────────────────────────────────────────────────────────┘   │
│                             │                                                │
│                             ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐   │
│  │  Result for Spot 1:                                                  │   │
│  │  ┌─────────────────────────────────────┐                            │   │
│  │  │  Gene │ T-cell │ Macro │ B-cell │   │                            │   │
│  │  │ ──────┼────────┼───────┼────────│   │                            │   │
│  │  │ TRAC  │   145  │   3   │    2   │   │ ← Most to T-cell          │   │
│  │  │ CD68  │    2   │  16   │    2   │   │ ← Most to Macrophage      │   │
│  │  │ MS4A1 │    1   │   2   │   27   │   │ ← Most to B-cell          │   │
│  │  └─────────────────────────────────────┘                            │   │
│  └──────────────────────────────────────────────────────────────────────┘   │
│                             │                                                │
│                             ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐   │
│  │  OUTPUT: Deconvolved Layers (stored in AnnData)                      │   │
│  │                                                                      │   │
│  │  adata.layers['T-cell_genes_pass1'] = (N × G) T-cell expression     │   │
│  │  adata.layers['Macrophage_genes_pass1'] = (N × G) Macrophage expr   │   │
│  │  adata.layers['B-cell_genes_pass1'] = (N × G) B-cell expression     │   │
│  │  ...                                                                 │   │
│  │                                                                      │   │
│  │  Properties:                                                         │   │
│  │  • Count-preserving: Σ_j layers[j][i,k] = adata.X[i,k]             │   │
│  │  • Cell-type-specific: Each layer contains only that type's signal │   │
│  │  • Sparse: Many zeros (cell types don't express all genes)         │   │
│  └──────────────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Module 4: Anchored Program Discovery

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 4a: PER-CELL-TYPE PROGRAM DISCOVERY                                 │
│  ══════════════════════════════════════════                                  │
│                                                                             │
│  INPUT: Deconvolved expression layers from Module 3                         │
│         Cell type proportions Y (N × T)                                     │
│         One layer per cell type: (N × G)                                    │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  For each cell type anchor (e.g., T-cell, Macrophage, B-cell...):          │
│                                                                             │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  STEP 1: Extract Deconvolved Layer                                   │  │
│  │  ─────────────────────────────────                                   │  │
│  │                                                                      │  │
│  │  X_celltype = adata.layers['T-cell_genes_pass1']  (N × G)           │  │
│  │  weights = proportions['T-cell']                  (N × 1)           │  │
│  │                                                                      │  │
│  │  Filter to active spots: X_active = X[weights > threshold]          │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                            │                                                │
│                            ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  STEP 2: Weighted Non-Negative Matrix Factorization (NMF)           │  │
│  │  ─────────────────────────────────────────────────────────          │  │
│  │                                                                      │  │
│  │  X_weighted = X_celltype × diag(√weights)                           │  │
│  │                                                                      │  │
│  │  NMF decomposition (K programs):                                     │  │
│  │                                                                      │  │
│  │       X_weighted ≈ W × H                                             │  │
│  │                                                                      │  │
│  │       W: Gene loadings (G × K) - which genes define each program    │  │
│  │       H: Spot loadings (K × N) - how active each program per spot   │  │
│  │                                                                      │  │
│  │      ┌───────────────────────────────────────────────────────┐      │  │
│  │      │           X (N × G)                                   │      │  │
│  │      │   ┌───────────────────┐                               │      │  │
│  │      │   │     Gene 1...G    │                               │      │  │
│  │      │   │  ┌─────────────┐  │      W (G × K)    H (K × N)   │      │  │
│  │      │Spot│  │             │  │  =   ┌─────┐   × ┌─────────┐ │      │  │
│  │      │ 1 │  │   8.2  3.1  │  │      │     │     │  Prog 1 │ │      │  │
│  │      │ 2 │  │   2.1  6.5  │  │      │Gene │     │  Prog 2 │ │      │  │
│  │      │...│  │    ...      │  │      │loads│     │   ...   │ │      │  │
│  │      │   │  └─────────────┘  │      │     │     │  Prog K │ │      │  │
│  │      │   └───────────────────┘      └─────┘     └─────────┘ │      │  │
│  │      └───────────────────────────────────────────────────────┘      │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                            │                                                │
│                            ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  STEP 3: Spatial Smoothing (Optional Laplacian Regularization)      │  │
│  │  ───────────────────────────────────────────────────────────        │  │
│  │                                                                      │  │
│  │  H_smooth = (I + λ_spatial × L)⁻¹ × H                               │  │
│  │                                                                      │  │
│  │  Where L = Laplacian matrix from spatial neighbor graph             │  │
│  │  Encourages adjacent spots to have similar program activities       │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                            │                                                │
│                            ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  STEP 4: Program Characterization                                    │  │
│  │  ────────────────────────────                                        │  │
│  │                                                                      │  │
│  │  For each program k = 1...K:                                         │  │
│  │                                                                      │  │
│  │  • Top genes: argmax(W[:, k]) → top 50 genes                        │  │
│  │  • Variance explained: ||W[:,k] × H[k,:]||² / ||X||²                │  │
│  │  • Spatial Moran's I: autocorrelation of H[k,:] across spots        │  │
│  │    - I > 0: Program is spatially clustered                          │  │
│  │    - I ≈ 0: Program is randomly distributed                         │  │
│  │    - I < 0: Program is dispersed (rare)                             │  │
│  │  • Mean activity: mean(H[k,:])                                       │  │
│  │  • Active fraction: % spots with H[k,:] > median                    │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                            │                                                │
│                            ▼                                                │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  STEP 5: Subpopulation Detection                                     │  │
│  │  ───────────────────────────                                         │  │
│  │                                                                      │  │
│  │  Features = [H (program activities), spatial_coords]                 │  │
│  │  K-means clustering → spatially distinct subpopulations             │  │
│  │                                                                      │  │
│  │    Subpop 1: Spots 1-50 (left side, Program 1 dominant)             │  │
│  │    Subpop 2: Spots 51-100 (right side, Program 2 dominant)          │  │
│  │    ...                                                               │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: AnchoredProgramResult (per cell type)                              │
│  • anchor_name: "T-cell"                                                    │
│  • programs: List of SpatialProgram objects                                 │
│  • W: Gene loadings (G × K)                                                 │
│  • H: Spot activities (K × N)                                               │
│  • subpopulations: Spatially distinct groups                                │
└─────────────────────────────────────────────────────────────────────────────┘
```

**Module 4b: Bivariate Program Relationships:**

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 4b: BIVARIATE PROGRAM RELATIONSHIPS                                 │
│  ══════════════════════════════════════════                                  │
│                                                                             │
│  INPUT: H matrices from all cell type programs                              │
│         Spatial coordinates                                                 │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 1: Generate All Program Pairs                                         │
│  ──────────────────────────────────                                         │
│                                                                             │
│  Programs:                                                                  │
│    T-cell: [Prog_0, Prog_1, Prog_2]                                        │
│    Macrophage: [Prog_0, Prog_1]                                            │
│    B-cell: [Prog_0, Prog_1]                                                │
│                                                                             │
│  Pairs (all combinations):                                                  │
│    Same cell type: T-cell_0 ↔ T-cell_1, T-cell_0 ↔ T-cell_2, ...          │
│    Cross cell type: T-cell_0 ↔ Macro_0, T-cell_0 ↔ B-cell_0, ...          │
│                                                                             │
│  Total: n_programs × (n_programs - 1) / 2 pairs                            │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 2: Compute Bivariate Moran's I for Each Pair                          │
│  ─────────────────────────────────────────────────                          │
│                                                                             │
│  For pair (Program A, Program B):                                           │
│                                                                             │
│    H_A = activity vector for program A (N spots)                           │
│    H_B = activity vector for program B (N spots)                           │
│                                                                             │
│    Filter to co-active spots: where H_A > 0 AND H_B > 0                    │
│    Skip if < 20 co-active spots                                             │
│                                                                             │
│    Bivariate Moran's I:                                                     │
│                                                                             │
│    I_AB = mean_i( (H_A[i] - μ_A) × Σ_j w_ij × (H_B[j] - μ_B) )            │
│           ─────────────────────────────────────────────────────             │
│                            σ_A × σ_B                                        │
│                                                                             │
│    Interpretation:                                                          │
│    ┌─────────────────────────────────────────────────────────────────────┐ │
│    │  I_AB > 0.1:  CO-LOCALIZED                                          │ │
│    │               Programs peak together in same spatial regions        │ │
│    │               Example: Pro-inflammatory programs co-occurring       │ │
│    │                                                                     │ │
│    │  -0.1 < I < 0.1:  INDEPENDENT                                      │ │
│    │               No consistent spatial relationship                    │ │
│    │                                                                     │ │
│    │  I_AB < -0.1:  EXCLUSIVE                                           │ │
│    │               Programs avoid each other spatially                   │ │
│    │               Example: Cytotoxic vs regulatory T-cell programs     │ │
│    └─────────────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 3: Statistical Testing & FDR Correction                               │
│  ────────────────────────────────────────────                               │
│                                                                             │
│  Permutation test (199 permutations):                                       │
│    - Shuffle H_B, recompute I_AB                                           │
│    - p-value = (# shuffled I ≥ observed I) / n_perms                       │
│                                                                             │
│  FDR correction (Benjamini-Hochberg):                                       │
│    - Correct for multiple testing across all pairs                         │
│    - Keep pairs with q-value < 0.05                                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: BivariateProgramResult                                             │
│                                                                             │
│  • significant_pairs: FDR-significant program relationships                 │
│  • all_pairs: Complete pair information                                     │
│                                                                             │
│  Example network visualization:                                             │
│                                                                             │
│         [T-cell_Prog0]────────[Macro_Prog0]                                │
│              │      co-localized (I=0.35)                                   │
│              │                                                              │
│         exclusive                                                           │
│         (I=-0.22)                                                           │
│              │                                                              │
│              ▼                                                              │
│         [T-cell_Prog1]        [B-cell_Prog0]                               │
│                    └──────────────┘                                         │
│                     independent                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

### Module 5: Cross-Sample Integration

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  MODULE 5: CROSS-SAMPLE INTEGRATION                                         │
│  ══════════════════════════════════                                          │
│                                                                             │
│  INPUT: Module 4 + 4b results from multiple patient samples                 │
│         Sample 1: W matrices, H matrices, bivariate relationships          │
│         Sample 2: W matrices, H matrices, bivariate relationships          │
│         Sample 3: ...                                                       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 1: GENE SET ALIGNMENT                                                 │
│  ────────────────────────                                                   │
│                                                                             │
│  Problem: Different samples may have different gene sets                    │
│                                                                             │
│    Sample 1 genes: Gene_001...Gene_080 (80 genes)                          │
│    Sample 2 genes: Gene_010...Gene_100 (90 genes, overlapping)             │
│    Sample 3 genes: Gene_001...Gene_075 (75 genes)                          │
│                                                                             │
│  Solution: Create union and zero-pad                                        │
│                                                                             │
│    Union: Gene_001...Gene_100 (110 unique genes)                           │
│                                                                             │
│    W_sample1_aligned: (110 × K1) - padded with zeros                       │
│    W_sample2_aligned: (110 × K2) - padded with zeros                       │
│    W_sample3_aligned: (110 × K3) - padded with zeros                       │
│                                                                             │
│  Stack all programs: total_programs = Σ (K1 + K2 + K3) per cell type       │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 2: HARMONY-STYLE BATCH CORRECTION                                     │
│  ──────────────────────────────────────                                     │
│                                                                             │
│  Goal: Remove sample-specific batch effects while preserving biology       │
│                                                                             │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  1. PCA reduction: W_stacked (n_genes × n_programs) → Z (n_progs×30) │  │
│  │                                                                      │  │
│  │  2. L2 normalization of each program embedding                       │  │
│  │                                                                      │  │
│  │  3. K-means initialization: 50 clusters                              │  │
│  │                                                                      │  │
│  │  4. Iterative correction:                                            │  │
│  │     ┌─────────────────────────────────────────────────────────────┐ │  │
│  │     │  REPEAT until convergence:                                  │ │  │
│  │     │                                                             │ │  │
│  │     │    a. Compute global centroids C (50 × 30)                  │ │  │
│  │     │    b. For each sample s:                                    │ │  │
│  │     │       - Compute sample-specific centroids C_s               │ │  │
│  │     │       - Correction: Z_s += R_s × (C - C_s)                  │ │  │
│  │     │       - R_s = soft cluster assignments for sample s         │ │  │
│  │     │    c. Re-normalize embeddings                               │ │  │
│  │     │    d. Update cluster assignments                            │ │  │
│  │     │                                                             │ │  │
│  │     │  STOP when Δembeddings < 1e-4                               │ │  │
│  │     └─────────────────────────────────────────────────────────────┘ │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│  Result: Integrated embedding Z where similar programs cluster together    │
│          regardless of which sample they came from                          │
│                                                                             │
│  ┌────────────────────────────────────────┐                                │
│  │     BEFORE Harmony       AFTER Harmony │                                │
│  │                                        │                                │
│  │    ●○○                    ●○◆         │  ● = Sample 1                  │
│  │   ●●○ ◆                  ●○◆          │  ○ = Sample 2                  │
│  │    ●  ◆◆                  ●○◆         │  ◆ = Sample 3                  │
│  │                           ●○◆         │                                │
│  │  (separated by           (grouped by  │                                │
│  │   sample batch)           biology)    │                                │
│  └────────────────────────────────────────┘                                │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 3: PROGRAM MATCHING                                                   │
│  ────────────────────────                                                   │
│                                                                             │
│  1. Compute pairwise cosine similarity in integrated space                 │
│                                                                             │
│     Similarity(prog_i, prog_j) = Z_i · Z_j / (||Z_i|| × ||Z_j||)          │
│                                                                             │
│  2. Hierarchical clustering with distance threshold = 1 - sim_threshold    │
│                                                                             │
│  3. Group programs with similarity > 0.7                                    │
│                                                                             │
│     Aligned_001: [Sample1_T-cell_Prog0, Sample2_T-cell_Prog0,              │
│                   Sample3_T-cell_Prog0]                                     │
│                   → Conservation score = 3/3 = 1.0                          │
│                   → Consensus signature = mean(W vectors)                   │
│                   → Top genes from consensus                                │
│                                                                             │
│     Aligned_002: [Sample1_Macro_Prog1, Sample3_Macro_Prog0]                │
│                   → Conservation score = 2/3 = 0.67                         │
│                                                                             │
│     Aligned_003: [Sample2_B-cell_Prog0]  (sample-specific)                 │
│                   → Conservation score = 1/3 = 0.33                         │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 4: RELATIONSHIP CONSERVATION                                          │
│  ────────────────────────────────                                           │
│                                                                             │
│  Compare Module 4b relationships across samples:                            │
│                                                                             │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │  Sample 1: T-cell_Prog0 ↔ Macro_Prog0:  I = 0.35 (co-localized)     │  │
│  │  Sample 2: T-cell_Prog0 ↔ Macro_Prog0:  I = 0.32 (co-localized)     │  │
│  │  Sample 3: T-cell_Prog0 ↔ Macro_Prog0:  I = 0.29 (co-localized)     │  │
│  │                                                                      │  │
│  │  Map to aligned programs:                                            │  │
│  │    All three map to: Aligned_001 ↔ Aligned_004                       │  │
│  │                                                                      │  │
│  │  Conserved relationship:                                             │  │
│  │    Programs: Aligned_001 ↔ Aligned_004                               │  │
│  │    Mean I: 0.32                                                      │  │
│  │    Std I: 0.03 (very consistent!)                                   │  │
│  │    Conservation: 3/3 = 1.0                                          │  │
│  │    Type: "co-localized"                                              │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
│                                                                             │
│  Classification:                                                            │
│    Mean I > 0.1:  "co-localized" (programs peak together)                  │
│    Mean I < -0.1: "exclusive" (programs avoid each other)                  │
│    Otherwise:     "independent"                                             │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  STEP 5: SIMILARITY NETWORK CONSTRUCTION                                    │
│  ───────────────────────────────────────                                    │
│                                                                             │
│  Build NetworkX graph:                                                      │
│                                                                             │
│    Nodes: Aligned programs (with conservation ≥ 0.3)                       │
│      - Attributes: cell_type, conservation_score, top_genes                │
│                                                                             │
│    Edges: Conserved relationships (with conservation ≥ 0.3)               │
│      - Weight: mean_bivariate_I                                            │
│      - Type: co-localized / exclusive / independent                        │
│                                                                             │
│  ┌──────────────────────────────────────────────────────────────────────┐  │
│  │                                                                      │  │
│  │         [T-cell_Aligned_001]═══════════[Macro_Aligned_004]          │  │
│  │         (conserv=1.0)        co-loc     (conserv=1.0)               │  │
│  │              ║                I=0.32                                │  │
│  │              ║ exclusive                                            │  │
│  │              ║ I=-0.18                                              │  │
│  │              ▼                                                       │  │
│  │         [T-cell_Aligned_002]           [B-cell_Aligned_005]         │  │
│  │         (conserv=0.67)                 (conserv=0.67)               │  │
│  │                                                                      │  │
│  │  Legend: ═══ co-localized, ─── independent, ║ exclusive            │  │
│  │          Node size ~ conservation score                              │  │
│  │          Edge width ~ |bivariate I|                                  │  │
│  └──────────────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT: IntegrationResult                                                  │
│                                                                             │
│  • aligned_programs: List of AlignedProgram objects                        │
│    - program_id, cell_type, samples_present                                │
│    - consensus_signature (meta-analyzed gene loadings)                      │
│    - conservation_score (0-1)                                               │
│    - top_genes                                                              │
│                                                                             │
│  • conserved_relationships: List of ConservedRelationship objects          │
│    - program1_id, program2_id                                               │
│    - mean_bivariate_i, std_bivariate_i                                      │
│    - relationship_type                                                      │
│    - conservation_score                                                     │
│                                                                             │
│  • similarity_graph: NetworkX graph for visualization                       │
│                                                                             │
│  • integration_embedding: 30D embedding for all programs                    │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Part 3: Data Flow Summary Table

| Module | Step | Input Shape | Output Shape | Key Data Structure |
|--------|------|-------------|--------------|-------------------|
| **1** | Marker Interest | N×M (spots×markers) | M' (ranked markers) | `MarkerInterestResult` |
| **2a** | Colocalization | M'×M' pairs | M'×(M'-1)/2 scored pairs | `ColocalizationResult` |
| **2b** | Profile Discovery | Scored pairs | P profiles (lists of markers) | `ProfileDiscoveryResult` |
| **2c** | Profile Selection | P profiles | P' selected profiles (5-15) | `ProfileSelectionResult` |
| **3.1** | Proportions | N×M', M'→T mapping | N×T proportions, M' betas | `cell_prop_results.csv` |
| **3.2** | Expression | N×G counts, N×T props | T layers of N×G | AnnData layers |
| **4a** | NMF Programs | N×G per cell type | W(G×K), H(K×N) per type | `AnchoredProgramResult` |
| **4b** | Bivariate | All H matrices | Program pair relationships | `BivariateProgramResult` |
| **5** | Integration | S samples × programs | Aligned programs, network | `IntegrationResult` |

**Legend:**
- N = number of spots
- M = number of markers (antibodies)
- M' = number of interesting markers
- G = number of genes
- T = number of cell types
- P = number of discovered profiles
- K = number of programs per cell type
- S = number of samples

---

## Part 4: BioRender Suggestions

### Recommended Figure Layout

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  PANEL A: Pipeline Overview (simplified version of Part 1)                 │
│                                                                             │
│  Horizontal flow with 5 module boxes, minimal text                         │
│  Use icons: tissue section → filters → clustering → optimization → network│
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│  PANEL B: Input Data Schematic                                             │
│                                                                             │
│  Left: Tissue section with spots colored by expression                     │
│  Right: Matrices (antibody, gene expression)                               │
└─────────────────────────────────────────────────────────────────────────────┘

┌────────────────────────────────┐  ┌────────────────────────────────────────┐
│  PANEL C: Module 1-2           │  │  PANEL D: Module 3                     │
│  Marker selection & profiles   │  │  Deconvolution schematic               │
│                                │  │                                        │
│  Show: markers → colocalization│  │  Show: mixed spot → separated layers  │
│  graph → profiles              │  │  Per cell type expression             │
└────────────────────────────────┘  └────────────────────────────────────────┘

┌────────────────────────────────┐  ┌────────────────────────────────────────┐
│  PANEL E: Module 4             │  │  PANEL F: Module 5                     │
│  NMF & relationships           │  │  Cross-sample integration             │
│                                │  │                                        │
│  Show: W×H decomposition,      │  │  Show: multiple samples → aligned     │
│  program network               │  │  programs, conserved network          │
└────────────────────────────────┘  └────────────────────────────────────────┘
```

### Color Scheme Recommendations

| Element | Color | Hex Code |
|---------|-------|----------|
| T-cells | Blue | #3498DB |
| Macrophages | Orange | #E67E22 |
| B-cells | Purple | #9B59B6 |
| Epithelial | Green | #27AE60 |
| Fibroblast | Brown | #8B4513 |
| Unknown/Other | Gray | #95A5A6 |
| Co-localized edge | Green | #2ECC71 |
| Exclusive edge | Red | #E74C3C |
| Independent edge | Gray | #BDC3C7 |

### Key Visual Metaphors

1. **Tissue section**: Hexagonal spots on tissue background
2. **Heatmaps**: For expression matrices and colocalization
3. **Dendrograms**: For hierarchical clustering in Module 2b
4. **Network graphs**: For program relationships (4b, 5)
5. **Stacked matrices**: For NMF decomposition (W × H)
6. **Venn-style overlaps**: For cross-sample program conservation

### Typography

- Module headers: Bold, 14pt
- Data annotations: Regular, 10pt
- Equations: Italic, 11pt
- Use arrows (→) liberally for flow direction

---

## Appendix: Key Equations

### Module 1: Interest Score
```
interest_score = kurtosis × GMM_SNR × max(Moran's_I, 0)
```

### Module 2a: Colocalization Score
```
score = 0.3×Spearman + 0.3×Cosine + 0.4×Bivariate_Moran's_I
```

### Module 3 Pass 1: Proportion Optimization
```
min Σ_i,m (S[i,m] - β[m]×Y[i,owner(m)])² + λ_reg×||Y||² + λ_lap×Σ_ij w_ij×||Y[i]-Y[j]||²
```

### Module 4a: NMF
```
X ≈ W × H, where W ≥ 0, H ≥ 0
```

### Module 4b/5: Bivariate Moran's I
```
I_AB = Σ_i Σ_j w_ij × (A_i - μ_A) × (B_j - μ_B) / (n × σ_A × σ_B)
```

---

*Last updated: 2026-01-13*
*For use as rough draft - final figure to be created in BioRender*
