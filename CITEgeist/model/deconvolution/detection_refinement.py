"""Detection refinement: GEX-based detection and iterative sparsity tuning.

This module extends the protein-only GMM detection (detection.py) with three
complementary capabilities:

1. **GEX detection** (Functions 1-2): Identifies cell type presence using gene
   expression correlations with protein markers. Useful when protein signal is
   ambiguous or markers are shared across types.

2. **Detection fusion** (Function 3): Combines protein and GEX detection masks
   using adaptive logic based on marker exclusivity.

3. **Sparsity refinement** (Function 4): Post-Pass-1 iterative tightening of
   detection upper bounds based on estimated proportions. Suppresses types that
   the QP solver assigned negligible weight, and rescues gated types that the
   solver found significant.
"""
import logging
from typing import List, Optional

import numpy as np
from sklearn.mixture import GaussianMixture

from .detection import _compute_adaptive_threshold

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Function 1: Gene-type correlation matrix
# ---------------------------------------------------------------------------

def compute_gene_type_correlations(
    Y: np.ndarray,
    antibody_data: np.ndarray,
    antibody_names: List[str],
    cell_profile_dict: dict,
    type_names: List[str],
) -> np.ndarray:
    """Compute correlation between each gene and each cell type's protein signal.

    For each type, averages its Major marker protein levels across spots, then
    computes Pearson r with each gene's expression. Negative correlations are
    clipped to 0 (not useful for detection). Types with no matching markers
    receive a uniform row.

    Args:
        Y: (N, G) raw GEX count matrix.
        antibody_data: (N, M) raw antibody capture matrix.
        antibody_names: (M,) marker names matching antibody_data columns.
        cell_profile_dict: {type_name: {"Major": [marker_name, ...]}}.
        type_names: (T,) ordered type names.

    Returns:
        H: (T, G) correlation matrix, values in [0, 1].
    """
    N, G = Y.shape
    T = len(type_names)
    H = np.zeros((T, G), dtype=float)

    # Build name -> column index map
    name_to_idx = {name: i for i, name in enumerate(antibody_names)}

    # Identify genes with nonzero total expression
    gene_totals = Y.sum(axis=0)
    nz_mask = gene_totals > 0

    for t, tname in enumerate(type_names):
        markers = cell_profile_dict.get(tname, {}).get("Major", [])
        marker_indices = [name_to_idx[m] for m in markers if m in name_to_idx]

        if not marker_indices:
            # No matching markers -> uniform row for nonzero genes
            n_nz = nz_mask.sum()
            if n_nz > 0:
                H[t, nz_mask] = 1.0 / n_nz
            logger.debug(f"{tname}: no matching markers, using uniform row")
            continue

        # Mean protein signal for this type across its markers
        prot_mean = antibody_data[:, marker_indices].mean(axis=1)  # (N,)

        # Vectorized Pearson correlation with all nonzero genes
        mc = prot_mean - prot_mean.mean()
        ms = prot_mean.std()
        if ms < 1e-10:
            n_nz = nz_mask.sum()
            if n_nz > 0:
                H[t, nz_mask] = 1.0 / n_nz
            continue

        Y_nz = Y[:, nz_mask]
        Yc = Y_nz - Y_nz.mean(axis=0)
        Ys = Y_nz.std(axis=0)
        Ys[Ys < 1e-10] = 1.0
        gene_corrs = (mc @ Yc) / (N * ms * Ys)
        H[t, nz_mask] = gene_corrs

    # Clip to [0, 1]
    np.clip(H, 0.0, 1.0, out=H)
    return H


# ---------------------------------------------------------------------------
# Function 2: GEX-based cell type detection
# ---------------------------------------------------------------------------

def detect_cell_types_gex(
    Y: np.ndarray,
    H: np.ndarray,
    gene_names: List[str],
    type_names: List[str],
    k: int = 10,
    min_corr: float = 0.15,
    threshold: float = 0.5,
) -> np.ndarray:
    """Detect cell type presence per spot using gene expression signatures.

    For each type, selects the top-k most correlated AND specific genes,
    computes a composite log-expression score, and fits a 2-component GMM
    to separate signal from background.

    Args:
        Y: (N, G) raw GEX count matrix.
        H: (T, G) gene-type correlation matrix from compute_gene_type_correlations.
        gene_names: (G,) gene names (for logging).
        type_names: (T,) type names.
        k: Maximum number of genes to select per type.
        min_corr: Minimum correlation threshold for gene inclusion.
        threshold: Base GMM posterior threshold.

    Returns:
        detected: (N, T) boolean detection mask.
    """
    N, G = Y.shape
    T = len(type_names)
    detected = np.ones((N, T), dtype=bool)  # default all-True

    # Precompute argmax type per gene (specificity filter)
    gene_best_type = np.argmax(H, axis=0)  # (G,)

    for t in range(T):
        # Gene selection: rank by H[t,g], filter by min_corr and specificity
        corrs = H[t]
        candidates = []
        for g in np.argsort(-corrs):
            if corrs[g] < min_corr:
                break
            if gene_best_type[g] != t:
                continue
            candidates.append(g)
            if len(candidates) >= k:
                break

        if len(candidates) < 3:
            logger.debug(
                f"{type_names[t]}: only {len(candidates)} qualifying genes, "
                f"falling back to all-True"
            )
            detected[:, t] = True
            continue

        # Composite score: mean log1p expression of selected genes
        selected = np.array(candidates)
        score = np.log1p(Y[:, selected]).mean(axis=1)  # (N,)

        # Fit 2-component GMM
        gmm = GaussianMixture(
            n_components=2,
            covariance_type="full",
            random_state=42,
            n_init=3,
        )
        try:
            gmm.fit(score.reshape(-1, 1))
        except Exception as e:
            logger.warning(f"GMM fit failed for {type_names[t]}: {e}. All-True fallback.")
            continue

        # Signal cluster = higher mean
        signal_cluster = int(np.argmax(gmm.means_.ravel()))

        # Adaptive threshold
        effective_threshold = _compute_adaptive_threshold(
            gmm, signal_cluster, base_threshold=threshold,
        )

        posteriors = gmm.predict_proba(score.reshape(-1, 1))[:, signal_cluster]
        detected[:, t] = posteriors > effective_threshold

        n_det = detected[:, t].sum()
        logger.debug(
            f"{type_names[t]}: {n_det}/{N} detected via GEX "
            f"({len(candidates)} genes, thresh={effective_threshold:.2f})"
        )

    return detected


# ---------------------------------------------------------------------------
# Function 3: Fuse protein and GEX detection masks
# ---------------------------------------------------------------------------

def fuse_detection_masks(
    protein_detected: np.ndarray,
    gex_detected: np.ndarray,
    assignment_matrix: np.ndarray,
    mode: str = "adaptive",
) -> np.ndarray:
    """Fuse protein-based and GEX-based detection masks.

    Args:
        protein_detected: (N, T) bool from protein GMM detection.
        gex_detected: (N, T) bool from GEX GMM detection.
        assignment_matrix: (M, T) marker-type assignment matrix (used in adaptive mode).
        mode: Fusion strategy -- "union", "intersection", or "adaptive".

    Returns:
        fused: (N, T) boolean detection mask.
    """
    N, T = protein_detected.shape

    if mode == "union":
        fused = protein_detected | gex_detected
    elif mode == "intersection":
        fused = protein_detected & gex_detected
    elif mode == "adaptive":
        # Count exclusive markers per type
        # Exclusive = marker assigned to exactly 1 type (column sum == 1)
        marker_type_count = (assignment_matrix > 0).sum(axis=1)  # (M,)
        exclusive_per_type = np.zeros(T, dtype=int)
        for t in range(T):
            exclusive_per_type[t] = np.sum(
                (assignment_matrix[:, t] > 0) & (marker_type_count == 1)
            )

        fused = np.empty((N, T), dtype=bool)
        for t in range(T):
            if exclusive_per_type[t] >= 2:
                fused[:, t] = protein_detected[:, t] | gex_detected[:, t]
            else:
                fused[:, t] = protein_detected[:, t] & gex_detected[:, t]

        logger.debug(
            f"Adaptive fusion: exclusive markers per type = "
            f"{dict(zip(range(T), exclusive_per_type))}"
        )
    else:
        raise ValueError(f"Unknown mode: {mode!r}. Use 'union', 'intersection', or 'adaptive'.")

    # Rescue: ensure >= 2 types detected per spot
    for i in range(N):
        if fused[i].sum() < 2:
            fused[i] = True

    return fused


# ---------------------------------------------------------------------------
# Function 4: Iterative sparsity refinement from proportions
# ---------------------------------------------------------------------------

def refine_sparsity_from_proportions(
    Y: np.ndarray,
    sparsity_mask: np.ndarray,
    cellularity: Optional[np.ndarray] = None,
    suppress_threshold: float = 0.02,
    rescue_threshold: float = 0.08,
    detection_gate_ub: float = 0.01,
) -> np.ndarray:
    """Refine detection upper bounds using Pass 1 proportion estimates.

    Tightens bounds for types the solver assigned negligible weight (suppress),
    and opens bounds for gated types the solver found significant (rescue).

    Args:
        Y: (N, T) proportions from Pass 1.
        sparsity_mask: (N, T) current upper bounds (1.0 = ungated, <1.0 = gated).
        cellularity: (N,) optional nuclei counts per spot.
        suppress_threshold: Ungated types with proportion below this are suppressed.
        rescue_threshold: Gated types with proportion above this are rescued.
        detection_gate_ub: Upper bound applied to suppressed types.

    Returns:
        refined: (N, T) refined upper bounds (copy; input not modified).
    """
    refined = sparsity_mask.copy()
    N, T = refined.shape

    for i in range(N):
        for t in range(T):
            is_ungated = sparsity_mask[i, t] >= 1.0
            is_gated = sparsity_mask[i, t] < 1.0

            if is_ungated and Y[i, t] < suppress_threshold:
                # Suppress: tighten upper bound
                if cellularity is not None:
                    refined[i, t] = min(detection_gate_ub, 1.0 / cellularity[i])
                else:
                    refined[i, t] = detection_gate_ub
            elif is_gated and Y[i, t] > rescue_threshold:
                # Rescue: open upper bound
                refined[i, t] = 1.0

        # Rescue check: if <2 types with mask >= 1.0, revert this spot
        if (refined[i] >= 1.0).sum() < 2:
            refined[i] = sparsity_mask[i]

    return refined


# ---------------------------------------------------------------------------
# Function 5: Entropy-based marker weighting
# ---------------------------------------------------------------------------

def compute_marker_entropy_weights(
    marker_level_data: np.ndarray,
    marker_names: List[str],
    alpha: float = 1.0,
    eps: float = 1e-10,
    weight_floor: float = 0.1,
) -> np.ndarray:
    """Compute per-marker weights based on expression entropy.

    Markers with concentrated spatial expression (low entropy, high specificity)
    get high weights. Diffuse markers (high entropy) get low weights.
    Inspired by UCASpatial's Shannon entropy-based gene weighting.

    Args:
        marker_level_data: (N, M) marker signal matrix.
        marker_names: (M,) marker names for logging.
        alpha: Sharpness exponent. 1.0=linear, 2.0=quadratic penalty.
        eps: Numerical floor for log computation.
        weight_floor: Minimum weight (default 0.1).

    Returns:
        (M,) weights in [weight_floor, 1.0].
    """
    N, M = marker_level_data.shape
    log_N = np.log(max(N, 2))

    weights = np.ones(M, dtype=np.float64)
    for m in range(M):
        signal = np.maximum(marker_level_data[:, m], eps)
        p = signal / signal.sum()
        H = -np.sum(p * np.log(p))
        H_norm = H / log_N
        weights[m] = (1.0 - H_norm) ** alpha

    weights = np.maximum(weights, weight_floor)

    for m in range(M):
        logger.info("Entropy weight: %s = %.3f (H_norm=%.3f)", marker_names[m], weights[m], 1.0 - weights[m] ** (1.0 / max(alpha, 0.01)))

    return weights
