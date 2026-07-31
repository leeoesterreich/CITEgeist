"""Spatially Anchored Cell Expression (SACE) for per-cell GEX deconvolution.

Single-pass Poisson-multinomial allocation of spot-level GEX counts to
individual cells, guided by QP proportions and marker-derived gene profiles.

Grid search (1125 configs × 5 Xenium RCC regions) showed all internal
parameters are inert at single-pass: eps, n_0, and bandwidth affect only
the M-step diagnostic output, not the allocation itself.
These are hardcoded below.
"""

import logging
from typing import Dict, List, NamedTuple, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import KDTree

from ..assignment.cellularity_utils import round_counts_largest_remainder

logger = logging.getLogger(__name__)


class SaceInternals(NamedTuple):
    """Internal SACE arrays needed by project_sace_to_cells."""

    B: "np.ndarray"  # (N, T) — type mass per spot
    mu: "np.ndarray"  # (N, T, G) — gene profiles per spot-type
    K: "sp.csr_matrix"  # (N, N) — spatial kernel


# Internal constants — validated as inert via grid search (range < 0.007 r).
# With single-pass allocation, these only affect .obs diagnostic columns
# (shrinkage_alpha, n_eff_type), not the actual GEX predictions.
_EPS = 1e-8
_N_0 = 10.0  # Fixed smoothing constant for effective-count shrinkage
_BANDWIDTH = None  # auto-computed from median nearest-neighbor distance


def build_kernel_matrix(
    coords: np.ndarray,
    bandwidth: Optional[float] = None,
    truncate_at: float = 3.0,
) -> sp.csr_matrix:
    """Build row-normalized Gaussian spatial kernel matrix.

    Args:
        coords: (N, 2) spot spatial coordinates.
        bandwidth: Gaussian kernel bandwidth h. If None, uses
            median nearest-neighbor distance * 2.
        truncate_at: Truncate kernel at this many bandwidths (default 3.0).

    Returns:
        (N, N) sparse CSR matrix, row-normalized so each row sums to 1.
    """
    N = coords.shape[0]
    tree = KDTree(coords)

    if bandwidth is None:
        dists, _ = tree.query(coords, k=2)
        nn_dists = dists[:, 1]
        bandwidth = float(np.median(nn_dists)) * 2.0
        logger.info("Auto bandwidth: %.1f (median NN dist %.1f)", bandwidth, np.median(nn_dists))

    radius = bandwidth * truncate_at
    pairs = tree.query_pairs(r=radius, output_type="ndarray")

    if len(pairs) == 0:
        return sp.eye(N, format="csr")

    ii = pairs[:, 0]
    jj = pairs[:, 1]
    d = np.linalg.norm(coords[ii] - coords[jj], axis=1)
    w = np.exp(-0.5 * (d / bandwidth) ** 2)

    rows = np.concatenate([ii, jj, np.arange(N)])
    cols = np.concatenate([jj, ii, np.arange(N)])
    vals = np.concatenate([w, w, np.ones(N)])

    K = sp.csr_matrix((vals, (rows, cols)), shape=(N, N))

    row_sums = np.array(K.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    inv_sums = sp.diags(1.0 / row_sums)
    K = inv_sums @ K

    return K


def _e_step(
    Y: np.ndarray,
    B: np.ndarray,
    mu: np.ndarray,
) -> np.ndarray:
    """E-step: allocate spot counts to types proportionally.

    Args:
        Y: (N, G) spot counts.
        B: (N, T) total RNA mass per type per spot.
        mu: (N, T, G) local type profiles (normalized per (s,t)).

    Returns:
        (N, T, G) array. Invariant: sum over t equals Y exactly.
    """
    weights = B[:, :, None] * mu
    lam = weights.sum(axis=1)
    lam_safe = np.where(lam > 0, lam, 1.0)
    E_x = Y[:, None, :] * weights / lam_safe[:, None, :]
    E_x = np.where(lam[:, None, :] > 0, E_x, 0.0)
    return E_x


def _m_step_B(E_x: np.ndarray) -> np.ndarray:
    """M-step: update total RNA mass B_{st} = sum_g E[x_{s,t,g}]."""
    return E_x.sum(axis=2)


def _m_step_mu(
    E_x: np.ndarray,
    K: sp.csr_matrix,
    *,
    min_mass_threshold: float = 1.0,
    mu_prev: Optional[np.ndarray] = None,
) -> np.ndarray:
    """M-step: update spatially-adaptive type profiles with shrinkage.

    Uses Kish ESS for shrinkage: n_eff = (Sum K*B)^2 / Sum (K*B)^2.
    Rare types (total mass < threshold) are frozen at previous/global profile.

    Note: With single-pass allocation (production default), this output only
    populates per-cell diagnostic columns — it does not feed back into allocation.
    """
    N, T, G = E_x.shape

    # Global profiles
    global_num = E_x.sum(axis=0) + _EPS
    global_den = global_num.sum(axis=1, keepdims=True)
    mu_global = global_num / global_den

    # Local profiles via kernel smoothing
    mu_local = np.zeros((N, T, G))
    for t in range(T):
        smoothed = K @ E_x[:, t, :]
        smoothed_total = smoothed.sum(axis=1, keepdims=True)
        smoothed_total = np.where(smoothed_total > 0, smoothed_total, 1.0)
        mu_local[:, t, :] = (smoothed + _EPS) / (smoothed_total + G * _EPS)

    # Kish ESS shrinkage — element-wise (K*B)^2 in denominator
    B = E_x.sum(axis=2)
    alpha = np.zeros((N, T))
    frozen_types = set()
    total_mass_per_type = B.sum(axis=0)
    for t in range(T):
        if total_mass_per_type[t] < min_mass_threshold:
            alpha[:, t] = 0.0
            frozen_types.add(t)
            continue
        b_t = B[:, t]
        KB = np.array(K @ b_t).ravel()
        K_bt = K.multiply(b_t[np.newaxis, :])
        KB2 = np.array(K_bt.multiply(K_bt).sum(axis=1)).ravel()
        n_eff = np.where(KB2 > 0, KB**2 / KB2, 0.0)
        alpha[:, t] = n_eff / (n_eff + _N_0)

    # Blend local + global
    mu = alpha[:, :, None] * mu_local + (1 - alpha[:, :, None]) * mu_global[None, :, :]

    # Normalize
    mu_sum = mu.sum(axis=2, keepdims=True)
    mu_sum = np.where(mu_sum > 0, mu_sum, 1.0)
    mu = mu / mu_sum

    # Freeze rare types
    for t in frozen_types:
        if mu_prev is not None:
            mu[:, t, :] = mu_prev[:, t, :]
        else:
            mu[:, t, :] = mu_global[t, :]
        logger.debug("Type %d frozen (mass %.2f < %.2f)", t, total_mass_per_type[t], min_mass_threshold)

    return mu


def _poisson_log_likelihood(
    Y: np.ndarray,
    B: np.ndarray,
    mu: np.ndarray,
) -> float:
    """Observed-data Poisson log-likelihood (up to constant)."""
    lam = (B[:, :, None] * mu).sum(axis=1)
    lam_safe = np.where(lam > 0, lam, 1e-30)
    ll = (Y * np.log(lam_safe) - lam).sum()
    return float(ll)


def _marker_guided_init(
    Y: np.ndarray,
    antibody_data: np.ndarray,
    antibody_names: List[str],
    cell_profile_dict: dict,
    type_names: List[str],
) -> tuple:
    """Marker-guided initialization for SACE mu_global.

    For each cell type, computes spatial Pearson correlation between each
    gene's expression and the type's mean marker protein expression across
    spots. Genes that co-vary spatially with a type's protein markers get
    higher weight in that type's gene profile, producing genuinely distinct
    profiles (mean row corr ~0.42 vs ~0.999 for confounded init).

    Args:
        Y: (N, G) raw spot counts.
        antibody_data: (N, M) raw antibody capture matrix.
        antibody_names: List of antibody/marker names matching antibody_data columns.
        cell_profile_dict: Dict mapping type name to marker dict,
            e.g. {"TypeA": {"Major": ["marker1", "marker2"]}}.
        type_names: Ordered type names matching proportions columns.

    Returns:
        Tuple of (mu_global, diagnostics) where:
        - mu_global: (T, G) normalized type profiles, or None if degenerate
        - diagnostics: dict with status and mean_row_corr
    """
    _, G = Y.shape
    T = len(type_names)

    from ..deconvolution.detection_refinement import (  # pylint: disable=import-outside-toplevel
        compute_gene_type_correlations,
    )

    H = compute_gene_type_correlations(
        Y,
        antibody_data,
        antibody_names,
        cell_profile_dict,
        type_names,
    )

    # SACE needs 0.01 floor for profile init (shared helper clips to [0,1])
    gene_mask = Y.sum(axis=0) > 0
    H = np.maximum(H, 0.01)
    H[:, ~gene_mask] = _EPS

    # Check degeneracy (on nonzero genes only, matching original)
    H_nz = H[:, gene_mask]
    H_norm = H_nz / (np.linalg.norm(H_nz, axis=1, keepdims=True) + _EPS)
    corr_matrix = H_norm @ H_norm.T
    mask = ~np.eye(T, dtype=bool)
    mean_corr = float(corr_matrix[mask].mean())

    if np.any(np.isnan(H)) or mean_corr > 0.99:
        logger.warning("Marker-guided init degenerate (mean_corr=%.4f)", mean_corr)
        return None, {"status": "degenerate", "mean_row_corr": mean_corr}

    # Build mu_global: normalize nonzero genes, eps-fill zero genes
    mu_global = np.full((T, G), _EPS)
    H_nz_normalized = H_nz / (H_nz.sum(axis=1, keepdims=True) + _EPS)
    mu_global[:, gene_mask] = H_nz_normalized

    G_nz = int(gene_mask.sum())
    logger.info("Marker-guided init: mean_row_corr=%.4f, G_nz=%d", mean_corr, G_nz)
    return mu_global, {"status": "ok", "mean_row_corr": mean_corr, "G_nonzero": G_nz}


def _largest_remainder_round(values: np.ndarray, total: int) -> np.ndarray:
    """Round non-negative floats to integers summing to total."""
    return round_counts_largest_remainder(np.asarray(values, dtype=float), int(total))


def integerize_spot_type_gex(
    spot_type_gex: Dict[int, np.ndarray],
    spot_counts: np.ndarray,
) -> Dict[int, np.ndarray]:
    """Round continuous SACE allocations to integers conserving spot totals.

    For each spot s and gene g, applies largest-remainder rounding to the
    (T,) allocation vector so that allocations are non-negative integers
    summing to exactly spot_counts[s, g].

    Args:
        spot_type_gex: Dict[spot_idx -> (T, G) ndarray], continuous.
        spot_counts: (N, G) integer count matrix (raw observed counts).

    Returns:
        Dict[spot_idx -> (T, G) ndarray] with integer values.
    """
    if hasattr(spot_counts, "toarray"):
        spot_counts = spot_counts.toarray()
    spot_counts = np.asarray(spot_counts, dtype=int)

    result = {}
    for s_idx, alloc in spot_type_gex.items():
        T, G = alloc.shape
        int_alloc = np.zeros((T, G), dtype=int)
        for g in range(G):
            obs_total = int(spot_counts[s_idx, g])
            type_vals = alloc[:, g]
            int_alloc[:, g] = _largest_remainder_round(type_vals, obs_total)
        result[s_idx] = int_alloc
    return result


def run_sace_layers(
    spot_counts: np.ndarray,
    proportions: pd.DataFrame,
    spot_coords: np.ndarray,
    *,
    gene_names: List[str],
    spotwise_profiles_init: Optional[Dict[int, np.ndarray]] = None,
    min_mass_threshold: float = 1.0,
    antibody_data: Optional[np.ndarray] = None,
    antibody_names: Optional[List[str]] = None,
    cell_profile_dict: Optional[dict] = None,
    init_method: str = "marker",
) -> Tuple[Dict[int, np.ndarray], np.ndarray, "SaceInternals", Dict]:
    """Core SACE E-step: allocate spot counts to types without cell projection.

    This is the compute-intensive layer of SACE.  It produces spot-level
    allocations and all internal arrays required for downstream per-cell
    projection via ``_assemble_cell_adata``.

    Args:
        spot_counts: (N, G) raw count matrix (NOT normalized).
        proportions: (N, T) DataFrame of cell type proportions from Pass 1.
        spot_coords: (N, 2) spot spatial coordinates.
        gene_names: List of gene names matching columns of spot_counts.
        spotwise_profiles_init: Optional Pass 1 deconvolved profiles for initialization.
        min_mass_threshold: Mass threshold for type freeze in diagnostics.
        antibody_data: Optional (N, M) raw antibody capture matrix for
            marker-guided initialization.
        antibody_names: List of antibody names matching antibody_data columns.
        cell_profile_dict: Dict mapping type name to marker dict for
            marker-guided initialization, e.g. {"TypeA": {"Major": [...]}}.
        init_method: Initialization strategy for mu_global. "marker" uses
            antibody-guided init (current default). "prop" uses proportion-weighted
            GEX (confounded proportional init from validated QP proportions).

    Returns:
        4-tuple of:
        - spot_type_gex: Dict[spot_idx -> (T, G) array]
        - E_x: (N, T, G) ndarray of allocated counts
        - internals: SaceInternals(B, mu, K) for downstream cell projection
        - diagnostics: Dict with log_likelihood and init_method
    """
    # Handle sparse matrices (common when raw_counts layer is CSR)
    if hasattr(spot_counts, "toarray"):
        spot_counts = spot_counts.toarray()
    Y = np.asarray(spot_counts, dtype=float)
    N, G = Y.shape
    type_names = list(proportions.columns)
    T = len(type_names)
    props = proportions.values.astype(float)  # (N, T)

    logger.info("SACE: N=%d spots, T=%d types, G=%d genes", N, T, G)

    # --- Build kernel matrix (for diagnostic M-step only) ---
    K = build_kernel_matrix(spot_coords, bandwidth=_BANDWIDTH)

    # --- Initialize mu and B ---
    lib_sizes = Y.sum(axis=1)  # (N,)

    if spotwise_profiles_init is not None:
        # Aggregate from Pass 1 spot-type profiles
        all_profiles = np.zeros((N, T, G))
        for s_idx, profile in spotwise_profiles_init.items():
            if s_idx < N:
                all_profiles[s_idx] = profile  # (T, G)
        # Global: mean across spots weighted by proportions
        mu_global = np.zeros((T, G))
        for t in range(T):
            weighted = props[:, t : t + 1] * all_profiles[:, t, :]
            total_w = props[:, t].sum()
            if total_w > 0:
                mu_global[t] = weighted.sum(axis=0) / total_w
            else:
                mu_global[t] = 1.0 / G
        # Normalize global profiles
        mu_global = (mu_global + _EPS) / (mu_global + _EPS).sum(axis=1, keepdims=True)
    elif init_method == "prop":
        logger.info("Using proportion-weighted init (confounded proportional)")
        mu_global = np.zeros((T, G))
        for t in range(T):
            weighted = props[:, t : t + 1] * Y
            total_w = weighted.sum()
            if total_w > 0:
                mu_global[t] = weighted.sum(axis=0) / total_w
            else:
                mu_global[t] = 1.0 / G
        mu_global = (mu_global + _EPS) / (mu_global + _EPS).sum(axis=1, keepdims=True)
    else:
        mu_global = None
        # Try marker-guided init if antibody data is available
        if antibody_data is not None and cell_profile_dict is not None:
            mg_result = _marker_guided_init(
                Y,
                antibody_data,
                antibody_names or [],
                cell_profile_dict,
                type_names,
            )
            mu_global, mg_diag = mg_result
            if mu_global is not None:
                logger.info(
                    "Using marker-guided init (mean_row_corr=%.4f)",
                    mg_diag.get("mean_row_corr", 0),
                )

        if mu_global is None:
            # Fallback: confounded proportional init
            logger.info("Using confounded proportional init (no antibody data or marker-guided failed)")
            mu_global = np.zeros((T, G))
            for t in range(T):
                weighted = props[:, t : t + 1] * Y
                total_w = weighted.sum()
                if total_w > 0:
                    mu_global[t] = weighted.sum(axis=0) / total_w
                else:
                    mu_global[t] = 1.0 / G
            mu_global = (mu_global + _EPS) / (mu_global + _EPS).sum(axis=1, keepdims=True)

    # Initialize local mu = global everywhere
    mu = np.broadcast_to(mu_global[None, :, :], (N, T, G)).copy()

    # Initialize B from proportions and library sizes
    B = props * lib_sizes[:, None]  # (N, T)

    # --- Single-pass allocation ---
    E_x = _e_step(Y, B, mu)
    ll = _poisson_log_likelihood(Y, B, mu)

    # M-step: populates per-cell diagnostics (entropy, Kish ESS, shrinkage
    # alpha) but does NOT affect the allocation itself at single-pass.
    mu = _m_step_mu(E_x, K, min_mass_threshold=min_mass_threshold, mu_prev=mu)
    mu_sum = mu.sum(axis=2, keepdims=True)
    mu_sum = np.where(mu_sum > 0, mu_sum, 1.0)
    mu = mu / mu_sum
    B = _m_step_B(E_x)

    logger.info("SACE single-pass: LL=%.1f", ll)

    # --- Build outputs ---
    spot_type_gex = {}
    for s in range(N):
        spot_type_gex[s] = E_x[s]

    internals = SaceInternals(B=B, mu=mu, K=K)
    diagnostics = {
        "log_likelihood": ll,
        "init_method": init_method,
    }

    return spot_type_gex, E_x, internals, diagnostics


def run_sace(
    spot_counts: np.ndarray,
    proportions: pd.DataFrame,
    cell_assignments: Dict[str, str],
    cell_spot_map: pd.DataFrame,
    spot_coords: np.ndarray,
    *,
    gene_names: List[str],
    spotwise_profiles_init: Optional[Dict[int, np.ndarray]] = None,
    min_mass_threshold: float = 1.0,
    antibody_data: Optional[np.ndarray] = None,
    antibody_names: Optional[List[str]] = None,
    cell_profile_dict: Optional[dict] = None,
    init_method: str = "marker",
) -> Tuple[Dict[int, np.ndarray], sc.AnnData, Dict]:
    """Backward-compatible wrapper: run SACE layers + cell projection.

    Allocates spot-level GEX counts to individual cells via Poisson-multinomial
    E-step using QP proportions and marker-guided gene profiles. Grid search
    (1125 configs × 5 Xenium RCC regions) validated that internal parameters
    (eps, n_0, bandwidth) are inert at single-pass — they affect only diagnostic
    metadata, not predictions.

    Args:
        spot_counts: (N, G) raw count matrix (NOT normalized).
        proportions: (N, T) DataFrame of cell type proportions from Pass 1.
        cell_assignments: Dict[cell_id -> type_name].
        cell_spot_map: DataFrame with columns: cell_id, spot_barcode, spot_idx, x, y.
        spot_coords: (N, 2) spot spatial coordinates.
        gene_names: List of gene names matching columns of spot_counts.
        spotwise_profiles_init: Optional Pass 1 deconvolved profiles for initialization.
        min_mass_threshold: Mass threshold for type freeze in diagnostics.
        antibody_data: Optional (N, M) raw antibody capture matrix for
            marker-guided initialization.
        antibody_names: List of antibody names matching antibody_data columns.
        cell_profile_dict: Dict mapping type name to marker dict for
            marker-guided initialization, e.g. {"TypeA": {"Major": [...]}}.
        init_method: Initialization strategy for mu_global. "marker" uses
            antibody-guided init (current default). "prop" uses proportion-weighted
            GEX (confounded proportional init from validated QP proportions).

    Returns:
        Tuple of:
        - spot_type_gex: Dict[spot_idx -> (T, G) array]
        - cell_adata: AnnData (n_cells, G)
        - diagnostics: Dict with convergence info
    """
    spot_type_gex, E_x, internals, diagnostics = run_sace_layers(
        spot_counts=spot_counts,
        proportions=proportions,
        spot_coords=spot_coords,
        gene_names=gene_names,
        spotwise_profiles_init=spotwise_profiles_init,
        min_mass_threshold=min_mass_threshold,
        antibody_data=antibody_data,
        antibody_names=antibody_names,
        cell_profile_dict=cell_profile_dict,
        init_method=init_method,
    )

    type_names = list(proportions.columns)
    cell_adata = project_sace_to_cells(
        E_x,
        internals,
        cell_assignments,
        cell_spot_map,
        gene_names,
        type_names,
    )

    return spot_type_gex, cell_adata, diagnostics


def project_sace_to_cells(  # pylint: disable=too-many-positional-arguments
    E_x: np.ndarray,
    sace_internals: "SaceInternals",
    cell_assignments: Dict[str, str],
    cell_spot_map: pd.DataFrame,
    gene_names: List[str],
    type_names: List[str],
) -> sc.AnnData:
    """Project spot-level SACE layers to per-cell AnnData.

    Each cell receives E_x[s,t,:] / n_st (equal split among cells of the
    same type in a spot). Orphan counts (types with mass but no cells)
    are redistributed proportionally.

    Args:
        E_x: (N, T, G) final allocated counts from run_sace_layers.
        sace_internals: SaceInternals namedtuple from run_sace_layers.
        cell_assignments: Dict mapping cell_id -> type_name.
        cell_spot_map: DataFrame with columns: cell_id, spot_barcode,
            spot_idx, x, y.
        gene_names: List of gene names matching the G axis of E_x.
        type_names: List of cell type names matching the T axis of E_x.

    Returns:
        AnnData of shape (n_cells, G) with per-cell expression, obs
        metadata, and obsm["spatial"] coordinates.
    """
    N, T, G = E_x.shape
    B = sace_internals.B
    mu = sace_internals.mu
    K = sace_internals.K
    type_to_idx = {t: i for i, t in enumerate(type_names)}

    cell_ids = cell_spot_map["cell_id"].values
    spot_indices = cell_spot_map["spot_idx"].values.astype(int)
    n_cells = len(cell_ids)

    # Count cells per (spot, type) for equal split
    n_st = np.zeros((N, T), dtype=int)
    cell_type_indices = np.zeros(n_cells, dtype=int)
    for i, cid in enumerate(cell_ids):
        t_name = cell_assignments[cid]
        t_idx = type_to_idx[t_name]
        s_idx = spot_indices[i]
        n_st[s_idx, t_idx] += 1
        cell_type_indices[i] = t_idx

    # Per-cell expression: E[x_c,g] = E[x_{s,t,g}] / n_{st}
    # For count conservation, redistribute orphan type counts (types with
    # nonzero E_x but zero cells) proportionally to existing cells.
    X_cells = np.zeros((n_cells, G))

    # First pass: compute per-spot orphan counts
    orphan = np.zeros((N, G))
    for s in range(N):
        for t in range(T):
            if n_st[s, t] == 0 and E_x[s, t, :].sum() > 0:
                orphan[s] += E_x[s, t, :]

    # Compute total cells per spot for redistribution
    n_cells_per_spot = n_st.sum(axis=1)  # (N,)

    for i in range(n_cells):
        s = spot_indices[i]
        t = cell_type_indices[i]
        count = n_st[s, t]
        if count > 0:
            X_cells[i] = E_x[s, t, :] / count
        # Add share of orphan counts
        if n_cells_per_spot[s] > 0:
            X_cells[i] += orphan[s] / n_cells_per_spot[s]

    # Diagnostics per cell
    # Allocation entropy: H of type weights for this cell's spot
    weights = B[:, :, None] * mu  # (N, T, G)
    spot_weights = weights.sum(axis=2)  # (N, T) total weight per type
    spot_weight_sums = spot_weights.sum(axis=1, keepdims=True)
    spot_weight_sums = np.where(spot_weight_sums > 0, spot_weight_sums, 1.0)
    p_types = spot_weights / spot_weight_sums  # (N, T) normalized
    p_safe = np.where(p_types > 0, p_types, 1.0)
    entropy = -(p_types * np.log(p_safe)).sum(axis=1)  # (N,)

    # Kish ESS per (spot, type) — must match _m_step_mu formula
    n_eff_arr = np.zeros((N, T))
    for t in range(T):
        b_t = B[:, t]
        KB = np.array(K @ b_t).ravel()
        K_bt = K.multiply(b_t[np.newaxis, :])
        KB2 = np.array(K_bt.multiply(K_bt).sum(axis=1)).ravel()
        n_eff_arr[:, t] = np.where(KB2 > 0, KB**2 / KB2, 0.0)

    # Shrinkage alpha per (spot, type)
    alpha_arr = n_eff_arr / (n_eff_arr + _N_0)

    obs_data = {
        "cell_id": cell_ids,
        "spot_barcode": cell_spot_map["spot_barcode"].values,
        "cell_type": [cell_assignments[cid] for cid in cell_ids],
        "allocation_entropy": entropy[spot_indices],
        "n_eff_type": [n_eff_arr[spot_indices[i], cell_type_indices[i]] for i in range(n_cells)],
        "shrinkage_alpha": [alpha_arr[spot_indices[i], cell_type_indices[i]] for i in range(n_cells)],
        "B_st": [B[spot_indices[i], cell_type_indices[i]] for i in range(n_cells)],
    }
    obs = pd.DataFrame(obs_data)
    obs.index = pd.Index(obs["cell_id"])

    coords = cell_spot_map[["x", "y"]].values.astype(float)

    adata = sc.AnnData(
        X=X_cells,
        obs=obs,
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm["spatial"] = coords

    return adata


# Backward-compatible alias (used internally and by external callers).
_assemble_cell_adata = project_sace_to_cells
