# CITEgeist — Automatic Multi-Order Antibody Profile Discovery (Adjusted Plan)

**Author:** A. Chang (edited with implementation notes)  
**Purpose:** Replace manual cell profile specification with statistically inferred antibody profiles (k=1–6)  
**Target:** Core model + CLI integration  
**Status:** Final design for coding agent

---

## Executive Summary
The revised plan delivers a fully unsupervised, threshold-free framework that automatically discovers the number and composition of antibody profiles. Key properties:
- Discovers the optimal number of profiles via BIC, requiring only 1–2 user-facing parameters.
- Supports mixed profile sizes (singles, doubles, triplets) in one loop.
- Handles shared markers through learned specificity weights (β) using EM, consistent with existing CITEgeist logic.
- Deterministic and reference-independent: all learning is driven by the antibody capture matrix.

---

## Alignment with Project Goals

### ✅ Strong Alignment
1. Addresses user burden by removing manual cell-type profile specification.
2. Remains reference-independent, learning directly from antibody data.
3. Keeps downstream Gurobi ILP compatibility (schema unchanged).
4. Aligns with CITEgeist’s EM/β approach already used in deconvolution.

---

## 1) Objectives & Scope
- Auto-discover single, double, and higher-order antibody profiles (default `max_k=3`; scalable to 6) with no manual thresholds.
- Only two user-facing parameters for discovery: `max_k` and `seed`.
- Output profiles in the existing dictionary schema so downstream Gurobi code remains unchanged.
- Deterministic, seeded, reproducible; no random sampling without explicit or default seed.

Non-goals: modifying Gurobi constraints, exhaustive combinatorics for k≥4, or introducing non-deterministic randomness.

---

## 2) File/Module Layout
- **New:** `CITEgeist/model/auto_profile_discovery.py` (export via `CITEgeist/model/__init__.py`).
- **New docs/tests:** `tests/test_profile_discovery.py`, doc updates here and brief note in `README.md` under development usage.
- **Touch:** `CITEgeist/model/citegeist_model.py` (new kwargs), CLI entrypoint (`examples/compute_sample.py` – current argparse path) for flags. If another canonical CLI exists, wire there too.

---

## 3) Data Interfaces & Contracts
- Input: `X` from `self.antibody_capture_adata.X` (can be dense or scipy sparse). Convert once with `np.asarray(X)` if not sparse, or `X.toarray()` when needed. Avoid copying twice.
- Marker names: `adata.var_names` list/Index.
- Preprocessing: use post-winsorization + CLR matrix, then z-score per marker (optionally robust median/MAD) prior to discovery.
- Output profiles:
  ```python
  {
    "EPCAM": {"Major": ["EPCAM"]},
    "CD3D_CD4": {"Major": ["CD3D", "CD4"]},
    "CD68_CD163_MRC1": {"Major": ["CD68", "CD163", "MRC1"]},
  }
  ```
- Merge behavior: when `profile_mode="auto"` or `--auto-profiles`, auto profiles replace manual ones; emit a log warning if both are provided. `profile_mode="manual"` (default) retains current behavior.

---

## 4) Math & Statistics: Unsupervised Profile Discovery

This framework discovers the optimal set of antibody profiles through iterative greedy search with EM refinement and automatic model selection via BIC.

### 4.1 Problem Formulation

The antibody data X (spots × markers) is modeled as:

    X ≈ Y @ A @ diag(β) + noise

Where:
- Y: (spots × K) — profile proportions per spot
- A: (K × markers) — binary profile definitions  
- β: (markers,) — marker specificity weights
- K: number of profiles (DISCOVERED, not specified)

### 4.2 Preprocessing

1. Input: X from `self.antibody_capture_adata.X` (after winsorization + CLR)  
2. Z-score standardization per marker:
   
       Z[s,i] = (X[s,i] - μ_i) / (σ_i + ε)
   
   Optional robust mode: use median/MAD instead of mean/std

### 4.3 Core Algorithm: Greedy Discovery with BIC Stopping

```python
def discover_profiles(Z, marker_names, max_k=3, seed=1234):
    """
    Iteratively discover profiles until BIC stops improving.
    
    Parameters (user-facing):
        Z: standardized antibody matrix (spots × markers)
        marker_names: list of marker names
        max_k: maximum profile size (default=3)
        seed: random seed for reproducibility (default=1234)
    
    Returns:
        profiles: dict of discovered profiles
        beta: dict of marker specificity weights
    """
    rng = np.random.default_rng(seed)
    n_spots, n_markers = Z.shape
    
    profiles = []           # List of marker index sets
    beta = np.ones(n_markers)
    prev_bic = np.inf
    
    while len(profiles) < 20:  # Safety cap
        
        # === STEP 1: Compute residual ===
        if profiles:
            Y, beta, _ = em_refine(Z, profiles, beta)
            A = build_profile_matrix(profiles, n_markers)
            residual = Z - Y @ A @ np.diag(beta)
        else:
            residual = Z
        
        # === STEP 2: Find best new profile ===
        best = find_best_candidate(residual, profiles, beta, max_k, rng)
        
        if best is None:
            break  # No significant candidates
        
        profiles.append(best)
        
        # === STEP 3: EM refinement ===
        Y, beta, log_lik = em_refine(Z, profiles, beta)
        
        # === STEP 4: BIC model selection ===
        n_params = len(profiles) * n_spots + n_markers
        bic = -2 * log_lik + n_params * np.log(n_spots)
        
        if bic > prev_bic:
            profiles.pop()  # Remove last, it made model worse
            break
        
        prev_bic = bic
    
    # Final refinement
    Y, beta, _ = em_refine(Z, profiles, beta)
    
    return format_output(profiles, beta, marker_names)
```

### 4.4 Finding the Best Candidate Profile

```python
def find_best_candidate(residual, existing, beta, max_k, rng, B=500):
    """
    Find the k-set that best explains residual variance.
    Must pass permutation significance test.
    """
    n_markers = residual.shape[1]
    existing_flat = set().union(*existing) if existing else set()
    
    best_score = -np.inf
    best_candidate = None
    
    # Generate candidates: all k-sets for k in 1..max_k
    # Skip markers already fully explained by existing profiles
    for k in range(1, max_k + 1):
        for candidate in combinations(range(n_markers), k):
            candidate = set(candidate)
            
            # Skip if this exact set already exists
            if candidate in existing:
                continue
            
            # Score: variance explained in residual
            score = score_candidate(residual, candidate, beta)
            
            # Permutation test
            null_scores = []
            for _ in range(B):
                perm_idx = rng.permutation(residual.shape[0])
                perm_residual = residual[perm_idx, :]
                null_scores.append(score_candidate(perm_residual, candidate, beta))
            
            p_value = (1 + np.sum(null_scores >= score)) / (1 + B)
            
            if p_value < 0.05 and score > best_score:
                best_score = score
                best_candidate = candidate
    
    return best_candidate


def score_candidate(Z, marker_set, beta):
    """
    Score = mean of β-weighted product of positive z-scores.
    High when all markers in set are jointly elevated.
    """
    markers = list(marker_set)
    # Product of max(z, 0) weighted by beta
    joint = np.ones(Z.shape[0])
    for m in markers:
        joint *= beta[m] * np.maximum(Z[:, m], 0)
    return joint.mean()
```

### 4.5 EM Refinement: Estimating Y and β

```python
def em_refine(Z, profiles, beta_init, max_iter=20, tol=1e-4):
    """
    Given profile definitions, estimate:
    - Y[s,k]: proportion of profile k in spot s  
    - β[i]: specificity weight for marker i
    
    Mirrors CITEgeist's existing EM for deconvolution.
    """
    n_spots, n_markers = Z.shape
    K = len(profiles)
    
    # Build binary profile matrix A
    A = build_profile_matrix(profiles, n_markers)
    
    beta = beta_init.copy()
    Y = np.full((n_spots, K), 1.0 / K)  # Uniform init
    
    prev_ll = -np.inf
    
    for _ in range(max_iter):
        
        # === E-step: Update Y given β ===
        # Soft assignment based on how well each profile explains each spot
        
        for k, profile in enumerate(profiles):
            markers = list(profile)
            # Likelihood: product of Gaussian pdfs for profile markers
            profile_fit = np.exp(-0.5 * np.sum(
                (Z[:, markers] - beta[markers]) ** 2, axis=1
            ))
            Y[:, k] = profile_fit
        
        # Normalize to proportions
        Y = Y / (Y.sum(axis=1, keepdims=True) + 1e-9)
        
        # === M-step: Update β given Y ===
        # β reflects marker specificity: high if marker discriminates profiles
        
        for i in range(n_markers):
            # Which profiles contain this marker?
            containing = [k for k in range(K) if A[k, i] > 0]
            
            if len(containing) == 0:
                beta[i] = 0.1  # Not in any profile
            elif len(containing) == 1:
                beta[i] = 1.0  # Unique to one profile — full weight
            else:
                # Shared marker: weight by concentration of Y when marker is high
                high_mask = Z[:, i] > 0
                if high_mask.sum() > 10:  # Need enough spots
                    Y_high = Y[high_mask][:, containing]
                    # Max proportion across profiles when marker is high
                    concentration = Y_high.max(axis=1).mean()
                    # Scale: if always one dominant profile → 1.0
                    #        if uniform across profiles → ~0.3
                    beta[i] = 0.3 + 0.7 * concentration
                else:
                    beta[i] = 0.5
        
        # === Log-likelihood for convergence ===
        expected = Y @ A @ np.diag(beta)
        ll = -0.5 * np.sum((Z - expected) ** 2)
        
        if abs(ll - prev_ll) < tol:
            break
        prev_ll = ll
    
    return Y, beta, ll
```

### 4.6 Output Formatting

```python
def format_output(profiles, beta, marker_names):
    """
    Convert internal representation to CITEgeist-compatible schema.
    """
    output = {}
    beta_out = {}
    
    for profile in profiles:
        markers = sorted([marker_names[i] for i in profile])
        name = "_".join(markers)
        output[name] = {"Major": markers}
    
    for i, name in enumerate(marker_names):
        beta_out[name] = float(beta[i])
    
    return output, beta_out
```

### 4.7 Helper Functions

```python
def build_profile_matrix(profiles, n_markers):
    """Build binary (K × markers) profile definition matrix."""
    K = len(profiles)
    A = np.zeros((K, n_markers))
    for k, profile in enumerate(profiles):
        for i in profile:
            A[k, i] = 1.0
    return A
```

### 4.8 Algorithm Properties

**Automatic K Discovery:**
- BIC penalizes complexity: -2*LL + n_params*log(n)
- Algorithm stops when adding profiles no longer improves BIC
- No need to specify number of profiles

**Mixed Profile Sizes:**
- Each iteration considers ALL k-sets (k=1 to max_k)
- Selects whichever best explains residual variance
- Naturally produces mixed sizes: {EPCAM}, {CD3D,CD4}, {CD68,CD163,MRC1}

**Shared Marker Handling:**
- β weights learned via EM
- Promiscuous markers (e.g., CD3D) get lower β
- Discriminating markers (e.g., CD4, CD8) retain high β
- Profiles scored with β-weighting, so specific profiles preferred

**Determinism:**
- Seeded RNG for permutation tests
- Deterministic tie-breaking: higher score, then lexicographic

### 4.9 Hyperparameter Summary

| Parameter | User-Facing | Default | Notes |
|-----------|-------------|---------|-------|
| max_k | Yes | 3 | Max markers per profile |
| seed | Yes | 1234 | Reproducibility |
| B | No | 500 | Permutations (internal) |
| max_iter | No | 20 | EM iterations (internal) |
| tol | No | 1e-4 | EM convergence (internal) |
| max_profiles | No | 20 | Safety cap (internal) |

**User API: 2 optional parameters**

---

## 5) Algorithm Steps (Implementation-Oriented)
1. Preprocess antibody matrix (winsorize + CLR already done) → z-score per marker (optionally robust median/MAD).
2. Initialize `profiles=[]`, `beta=ones`, `prev_bic=inf`, RNG with seed.
3. Loop (cap at `max_profiles`):
   - Compute residual using current profiles via `em_refine`; use raw Z if none exist.
   - Enumerate candidates (size 1..`max_k`), score with `score_candidate`, and require permutation significance.
   - If no candidate passes, stop.
   - Append best candidate, run `em_refine`, compute BIC; if BIC worsens, drop candidate and stop.
4. Final `em_refine`, then format outputs (profile dict + beta dict keyed by marker name).

---

## 6) API/Code Changes

```markdown
### New Module: CITEgeist/model/auto_profile_discovery.py

Exports:
- discover_profiles(X, marker_names, max_k=3, seed=1234) → (profiles, beta)
- em_refine(Z, profiles, beta) → (Y, beta, log_lik)  # For advanced users

### CitegeistModel Changes

class CitegeistModel:
    def __init__(
        self,
        ...,
        profile_mode: Literal["manual", "auto"] = "manual",
        max_profile_size: int = 3,
        rng_seed: int = 1234,
    ):
        self.profile_mode = profile_mode
        self.max_profile_size = max_profile_size
        self.rng_seed = rng_seed
    
    def _run_auto_profile_discovery(self):
        """Called after antibody preprocessing when profile_mode='auto'."""
        from .auto_profile_discovery import discover_profiles
        
        X = self.antibody_capture_adata.X
        marker_names = list(self.antibody_capture_adata.var_names)
        
        profiles, beta = discover_profiles(
            X, marker_names, 
            max_k=self.max_profile_size,
            seed=self.rng_seed
        )
        
        self.cell_type_profiles = profiles
        self.marker_beta = beta
        
        logger.info(f"Discovered {len(profiles)} profiles: {list(profiles.keys())}")

### CLI Changes (examples/compute_sample.py)

parser.add_argument("--auto-profiles", action="store_true",
                    help="Auto-discover antibody profiles (default: manual)")
parser.add_argument("--max-profile-size", type=int, default=3,
                    help="Max markers per profile for auto-discovery")
parser.add_argument("--seed", type=int, default=1234,
                    help="Random seed for reproducibility")

# Usage:
# python compute_sample.py --auto-profiles --max-profile-size 4
```

---

## 7) Determinism & Logging
- Seed handling: default to 1234 if none provided; record seed in logs and returned metadata.
- Tie-breakers: score desc, then lexicographic marker names.
- Log discovered profile count and names; warn when manual profiles are ignored due to auto mode.
- Log when no significant candidates are found or when BIC halts addition.

---

## 8) Complexity & Performance
- Candidate search explores all k-sets up to `max_k`; practical with `max_k<=3` and safety caps.
- Permutation test per candidate uses streaming counts to avoid large memory allocations.
- EM iterations bounded (`max_iter=20`) with simple Gaussian likelihood; reuse arrays when possible.
- Residual updates avoid redundant conversions; operate on dense NumPy arrays after a single sparse-to-dense conversion.

---

## 9) Failure Modes
- No significant candidate: stop gracefully and return profiles discovered so far (possibly zero).
- BIC worse after adding profile: revert last addition and stop.
- Shared marker ambiguity: β down-weights promiscuous markers via EM; still returns best-effort profiles.
- Limited markers (`n_markers < max_k`): automatically constrained by combinations; still deterministic.

---

## 10) Tests (`tests/test_profile_discovery.py`)

### Core Functionality

1. **test_single_dominant_population** — Input: one cell type with single marker (e.g., EPCAM+ cancer). Expected: discovers {EPCAM}.  
2. **test_double_profile** — Input: CD4 T cells with CD3D+CD4+ pattern. Expected: discovers {CD3D, CD4}.  
3. **test_mixed_profile_sizes** — Input: Cancer (EPCAM), T cells (CD3D+CD4), Macrophages (CD68+CD163+MRC1). Expected: returns 1 single, 1 double, 1 triplet.  
4. **test_shared_marker_discrimination** — Input: CD4 T cells (CD3D+CD4) and CD8 T cells (CD3D+CD8). Expected: both profiles discovered; β_CD3D < β_CD4 and β_CD3D < β_CD8.  
5. **test_bic_stops_overfitting** — Input: data with 3 true profiles + noise. Expected: discovers exactly 3 profiles, not more.  

### Edge Cases

6. **test_no_structure** — Random noise matrix → returns empty or minimal profiles; doesn't crash.  
7. **test_few_markers** — Only 2–3 markers → handles gracefully; `max_k` capped appropriately.  
8. **test_sparse_input** — scipy.sparse matrix → converts and processes correctly.  

### Determinism

9. **test_reproducibility** — Run twice with same seed → identical profiles and β values.  
10. **test_different_seeds** — Different seeds → same profiles (order may vary), similar β.  

### Integration

11. **test_output_schema_compatibility** — Output works with existing `run_cell_prop_model()`.  
12. **test_end_to_end_with_simulated_data** — Use existing fixtures; run full pipeline: auto-discovery → deconvolution → verify accuracy.  

---

## 11) Example: What Discovery Looks Like

```python
# Simulated tissue with 4 cell types
# Spots × Markers matrix after CLR normalization

>>> profiles, beta = discover_profiles(X, marker_names, max_k=3)

Iteration 1: Testing candidates...
  Best: {EPCAM}, score=0.82, p<0.001
  BIC: 15234.5
  
Iteration 2: Computing residual...
  Best: {CD3D, CD4}, score=0.45, p<0.001
  BIC: 12108.3 (improved)
  
Iteration 3: Computing residual...
  Best: {CD3D, CD8}, score=0.38, p<0.001
  BIC: 10892.1 (improved)
  
Iteration 4: Computing residual...
  Best: {CD68, CD163, MRC1}, score=0.31, p=0.003
  BIC: 10245.7 (improved)
  
Iteration 5: Computing residual...
  Best: {CD19}, score=0.12, p=0.08 (not significant)
  No valid candidates. Stopping.

>>> profiles
{
    "EPCAM": {"Major": ["EPCAM"]},
    "CD3D_CD4": {"Major": ["CD3D", "CD4"]},
    "CD3D_CD8": {"Major": ["CD3D", "CD8"]},
    "CD68_CD163_MRC1": {"Major": ["CD68", "CD163", "MRC1"]}
}

>>> beta
{
    "EPCAM": 1.0,    # Unique to cancer profile
    "CD3D": 0.58,    # Shared between T cell profiles
    "CD4": 1.0,      # Unique to CD4 T cells
    "CD8": 1.0,      # Unique to CD8 T cells
    "CD68": 0.85,    # Mostly specific to macrophages
    "CD163": 0.92,   # Specific to M2 macrophages
    "MRC1": 0.95,    # Specific to M2 macrophages
    "CD19": 0.1,     # Not in any profile
}
```

---

## 12) Summary of Changes from Original Plan

| Aspect | Original Plan | Revised Plan |
|--------|---------------|--------------|
| Threshold | Fixed 0.9 quantile | None — fully threshold-free |
| Scoring | Binary Jaccard | β-weighted product of z-scores |
| K discovery | Test all, return FDR-passing | Greedy + BIC model selection |
| Shared markers | Not handled | β-weighting via EM |
| Profile sizes | Separate phases per k | Mixed sizes in single loop |
| User params | 4–5 | 2 (max_k, seed) |

---

## 13) Conclusion
The revised framework is fully unsupervised, threshold-free, biologically grounded, and consistent with CITEgeist’s EM/β concepts while exposing only two optional user parameters. It is ready for implementation.

---

## 14) Execution Checklist
1. Add `CITEgeist/model/auto_profile_discovery.py` implementing the greedy EM+BIC workflow and helpers.  
2. Export helpers in `CITEgeist/model/__init__.py`.  
3. Wire `CitegeistModel` to use `profile_mode`, `max_profile_size`, and `rng_seed`; add `_run_auto_profile_discovery`.  
4. Update CLI (`examples/compute_sample.py`, plus any other official entrypoint) with `--auto-profiles`, `--max-profile-size`, and `--seed`.  
5. Add tests per Section 10.  
6. Update README docs to mention auto profiles and CLI flags.  
7. Confirm coverage and deterministic behavior end-to-end.  
