"""Spatially-Adaptive Compositional EM (SACE) for per-cell GEX deconvolution.

Replaces heuristic per-cell GEX allocation with iterative Poisson-multinomial
EM that learns spatially-varying type profiles and between-type RNA mass.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial import KDTree

logger = logging.getLogger(__name__)


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
    n_0: float = 10.0,
    eps: float = 1e-6,
    min_mass_threshold: float = 1.0,
    mu_prev: Optional[np.ndarray] = None,
) -> np.ndarray:
    """M-step: update spatially-adaptive type profiles with shrinkage.

    Uses Kish ESS for shrinkage: n_eff = (Sum K*B)^2 / Sum (K*B)^2.
    Rare types (total mass < threshold) are frozen at previous/global profile.
    """
    N, T, G = E_x.shape

    # Global profiles
    global_num = E_x.sum(axis=0) + eps
    global_den = global_num.sum(axis=1, keepdims=True)
    mu_global = global_num / global_den

    # Local profiles via kernel smoothing
    mu_local = np.zeros((N, T, G))
    for t in range(T):
        smoothed = K @ E_x[:, t, :]
        smoothed_total = smoothed.sum(axis=1, keepdims=True)
        smoothed_total = np.where(smoothed_total > 0, smoothed_total, 1.0)
        mu_local[:, t, :] = (smoothed + eps) / (smoothed_total + G * eps)

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
        alpha[:, t] = n_eff / (n_eff + n_0)

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
    *,
    eps: float = 1e-10,
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
        eps: Numerical floor.

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
    H = np.maximum(H, eps)
    # Zero-expression genes stay at eps (not 0.01) to match original behavior
    H[:, ~gene_mask] = eps

    # Check degeneracy (on nonzero genes only, matching original)
    H_nz = H[:, gene_mask]
    H_norm = H_nz / (np.linalg.norm(H_nz, axis=1, keepdims=True) + eps)
    corr_matrix = H_norm @ H_norm.T
    mask = ~np.eye(T, dtype=bool)
    mean_corr = float(corr_matrix[mask].mean())

    if np.any(np.isnan(H)) or mean_corr > 0.99:
        logger.warning("Marker-guided init degenerate (mean_corr=%.4f)", mean_corr)
        return None, {"status": "degenerate", "mean_row_corr": mean_corr}

    # Build mu_global: normalize nonzero genes, eps-fill zero genes
    mu_global = np.full((T, G), eps)
    H_nz_normalized = H_nz / (H_nz.sum(axis=1, keepdims=True) + eps)
    mu_global[:, gene_mask] = H_nz_normalized

    G_nz = int(gene_mask.sum())
    logger.info("Marker-guided init: mean_row_corr=%.4f, G_nz=%d", mean_corr, G_nz)
    return mu_global, {"status": "ok", "mean_row_corr": mean_corr, "G_nonzero": G_nz}


def run_sace(
    spot_counts: np.ndarray,
    proportions: pd.DataFrame,
    cell_assignments: Dict[str, str],
    cell_spot_map: pd.DataFrame,
    spot_coords: np.ndarray,
    *,
    gene_names: List[str],
    spotwise_profiles_init: Optional[Dict[int, np.ndarray]] = None,
    n_0: float = 10.0,
    bandwidth: Optional[float] = None,
    max_iter: int = 10,
    tol: float = 1e-4,
    eps: float = 1e-6,
    min_mass_threshold: float = 1.0,
    damping_eta: float = 1.0,
    antibody_data: Optional[np.ndarray] = None,
    antibody_names: Optional[List[str]] = None,
    cell_profile_dict: Optional[dict] = None,
) -> Tuple[Dict[int, np.ndarray], sc.AnnData, Dict]:
    """Run Spatially-Adaptive Compositional EM for per-cell GEX.

    Args:
        spot_counts: (N, G) raw count matrix (NOT normalized).
        proportions: (N, T) DataFrame of cell type proportions from Pass 1.
        cell_assignments: Dict[cell_id -> type_name].
        cell_spot_map: DataFrame with columns: cell_id, spot_barcode, spot_idx, x, y.
        spot_coords: (N, 2) spot spatial coordinates.
        gene_names: List of gene names matching columns of spot_counts.
        spotwise_profiles_init: Optional Pass 1 deconvolved profiles for initialization.
        n_0: Shrinkage prior strength.
        bandwidth: Kernel bandwidth (None=auto).
        max_iter: Maximum EM iterations.
        tol: Convergence threshold.
        eps: Dirichlet floor.
        min_mass_threshold: Mass threshold for type freeze.
        damping_eta: Update damping (1.0 = no damping).
        antibody_data: Optional (N, M) raw antibody capture matrix for
            marker-guided initialization.
        antibody_names: List of antibody names matching antibody_data columns.
        cell_profile_dict: Dict mapping type name to marker dict for
            marker-guided initialization, e.g. {"TypeA": {"Major": [...]}}.

    Returns:
        Tuple of:
        - spot_type_gex: Dict[spot_idx -> (T, G) array]
        - cell_adata: AnnData (n_cells, G)
        - diagnostics: Dict with convergence info
    """
    # Handle sparse matrices (common when raw_counts layer is CSR)
    if hasattr(spot_counts, "toarray"):
        spot_counts = spot_counts.toarray()
    Y = np.asarray(spot_counts, dtype=float)
    N, G = Y.shape
    type_names = list(proportions.columns)
    T = len(type_names)
    props = proportions.values.astype(float)  # (N, T)

    logger.info("SACE: N=%d spots, T=%d types, G=%d genes, %d cells", N, T, G, len(cell_assignments))

    # --- Build kernel matrix ---
    K = build_kernel_matrix(spot_coords, bandwidth=bandwidth)

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
        mu_global = (mu_global + eps) / (mu_global + eps).sum(axis=1, keepdims=True)
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
                eps=eps,
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
            mu_global = (mu_global + eps) / (mu_global + eps).sum(axis=1, keepdims=True)

    # Initialize local mu = global everywhere
    mu = np.broadcast_to(mu_global[None, :, :], (N, T, G)).copy()

    # Initialize B from proportions and library sizes
    B = props * lib_sizes[:, None]  # (N, T)

    # --- EM loop ---
    ll_trace = []
    change_trace = []  # type: ignore[var-annotated]
    damping_activated = False
    E_x_prev = None

    for iteration in range(max_iter):
        # E-step
        E_x = _e_step(Y, B, mu)

        # Convergence check
        if E_x_prev is not None:
            max_change = np.abs(E_x - E_x_prev).max()
            rel_change = max_change / max(E_x.max(), 1.0)
            change_trace.append(float(rel_change))
            if rel_change < tol:
                logger.info("SACE converged at iteration %d (change=%.2e)", iteration, rel_change)
                break
        E_x_prev = E_x.copy()

        # Log-likelihood
        ll = _poisson_log_likelihood(Y, B, mu)
        ll_trace.append(ll)

        # Check for likelihood drop -> enable damping
        if len(ll_trace) >= 2 and ll_trace[-1] < ll_trace[-2]:
            if damping_eta >= 1.0:
                damping_eta = 0.5
                damping_activated = True
                logger.warning(
                    "SACE: log-likelihood dropped at iter %d, " "enabling damping (eta=%.2f)", iteration, damping_eta
                )

        # M-step: profiles
        mu_new = _m_step_mu(E_x, K, n_0=n_0, eps=eps, min_mass_threshold=min_mass_threshold, mu_prev=mu)

        # M-step: RNA mass
        B_new = _m_step_B(E_x)

        # Apply damping
        mu = (1 - damping_eta) * mu + damping_eta * mu_new
        B = (1 - damping_eta) * B + damping_eta * B_new

        # Re-normalize mu after damping
        mu_sum = mu.sum(axis=2, keepdims=True)
        mu_sum = np.where(mu_sum > 0, mu_sum, 1.0)
        mu = mu / mu_sum

        logger.info("SACE iter %d: LL=%.1f%s", iteration, ll, f" change={change_trace[-1]:.2e}" if change_trace else "")

    # Final E-step with converged parameters.
    # Skip when max_iter <= 1: the loop's E-step already used the init mu,
    # and the M-step overwrote mu with confounded profiles. Re-running E-step
    # here would discard the marker-guided init signal.
    if max_iter > 1:
        E_x = _e_step(Y, B, mu)

    # --- Build outputs ---
    # 1. Spot-type profiles: Dict[spot_idx -> (T, G)]
    spot_type_gex = {}
    for s in range(N):
        spot_type_gex[s] = E_x[s]  # (T, G)

    # 2. Per-cell AnnData
    cell_adata = _assemble_cell_adata(
        E_x,
        B,
        mu,
        K,
        n_0,
        cell_assignments,
        cell_spot_map,
        gene_names,
        type_names,
    )

    # 3. Diagnostics
    diagnostics = {
        "log_likelihood_trace": ll_trace,
        "allocation_change": change_trace,
        "damping_activated": damping_activated,
        "n_iterations": len(ll_trace),
        "final_bandwidth": float(bandwidth) if bandwidth else "auto",
    }

    return spot_type_gex, cell_adata, diagnostics


def _assemble_cell_adata(  # pylint: disable=too-many-positional-arguments
    E_x: np.ndarray,  # (N, T, G) final allocated counts
    B: np.ndarray,  # (N, T)
    mu: np.ndarray,  # (N, T, G)
    K: sp.csr_matrix,  # (N, N)
    n_0: float,
    cell_assignments: Dict[str, str],
    cell_spot_map: pd.DataFrame,
    gene_names: List[str],
    type_names: List[str],
) -> sc.AnnData:
    """Assemble per-cell AnnData from SACE results."""
    N, T, G = E_x.shape
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
    alpha_arr = n_eff_arr / (n_eff_arr + n_0)

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
