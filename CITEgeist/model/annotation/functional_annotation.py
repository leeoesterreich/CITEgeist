"""Module 3.5: Functional protein annotation for CITEgeist.

Post-deconvolution annotation step that runs after Module 3 QP proportions.
For each (cell type, functional marker) pair defined in the active mask, learns
NB emission rates and gates spots as functionally active based on
observed-to-expected ratio.

Generative model:
    S[i,m] ~ NB(mu_im, r_m)
    mu_im = s_i * (sum_t p_it * lambda_tm + b_m)

where:
    - p_it are fixed QP proportions from Module 3 (NOT re-optimized here)
    - s_i are fixed size factors (NOT optimized)
    - lambda_tm is the emission rate for type t and marker m (active pairs only)
    - b_m is per-marker background
    - r_m is per-marker NB dispersion
"""
import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Marker-type association table
# ---------------------------------------------------------------------------

DEFAULT_FUNCTIONAL_TABLE: Dict[str, Dict] = {
    # --- Immune checkpoint / exhaustion ---
    "PDCD1": {
        "function": "T cell exhaustion (PD-1)",
        "active_types": ["CD8_T_Cells", "CD4_T_Cells"],
    },
    "CD274": {
        "function": "Immune evasion (PD-L1)",
        "active_types": ["Epithelial", "Macrophages"],
    },
    # --- Proliferation (pan-type) ---
    "PCNA": {
        "function": "Proliferation",
        "active_types": [
            "Endothelial", "Fibroblasts", "B_Cells", "Macrophages", "Monocytes",
            "CD8_T_Cells", "CD4_T_Cells", "Epithelial", "Dendritic_Cells",
        ],
    },
    # --- Survival / apoptosis ---
    "BCL2": {
        "function": "Anti-apoptosis / survival",
        "active_types": ["B_Cells", "CD4_T_Cells", "Epithelial"],
    },
    # --- APC activation ---
    "CD40": {
        "function": "APC activation / co-stimulation",
        "active_types": ["B_Cells", "Dendritic_Cells", "Macrophages"],
    },
    # --- Memory / co-stimulation ---
    "CD27": {
        "function": "Memory / co-stimulation",
        "active_types": ["CD8_T_Cells", "CD4_T_Cells", "B_Cells"],
    },
    "CCR7": {
        "function": "Naive/central memory",
        "active_types": ["CD4_T_Cells", "CD8_T_Cells", "Dendritic_Cells"],
    },
    # --- Follicular homing ---
    "CXCR5": {
        "function": "Follicular homing (Tfh)",
        "active_types": ["CD4_T_Cells", "B_Cells"],
    },
    # --- B cell state ---
    "MS4A1": {
        "function": "B cell maturity (CD20)",
        "active_types": ["B_Cells"],
    },
    "CR2": {
        "function": "B cell activation (CD21)",
        "active_types": ["B_Cells"],
    },
    "PAX5": {
        "function": "B cell lineage TF",
        "active_types": ["B_Cells"],
    },
    # --- Myeloid activation ---
    "ITGAM": {
        "function": "Myeloid activation (CD11b)",
        "active_types": ["Macrophages", "Monocytes", "Dendritic_Cells"],
    },
    # --- Minor cell-typing markers repurposed as functional ---
    "FCGR3A": {
        "function": "NK/monocyte activation (CD16)",
        "active_types": ["Macrophages", "Monocytes"],
    },
    "VIM": {
        "function": "Mesenchymal / EMT (Vimentin)",
        "active_types": ["Fibroblasts", "Epithelial", "Endothelial"],
    },
    # --- Epithelial functional state markers ---
    "KRT5": {
        "function": "Basal epithelial state (Keratin 5)",
        "active_types": ["Epithelial"],
    },
    "SDC1": {
        "function": "Luminal / plasma cell marker (Syndecan-1 / CD138)",
        "active_types": ["Epithelial"],
    },
    # NOTE: CEACAM8 (CD66b) excluded — neutrophil marker, no matching cell type
    # NOTE: PTPRC (CD45) excluded — pan-immune, indistinguishable isoforms
}


# ---------------------------------------------------------------------------
# build_active_mask
# ---------------------------------------------------------------------------

def build_active_mask(
    functional_markers: List[str],
    cell_types: List[str],
    functional_table: Optional[Dict] = None,
) -> np.ndarray:
    """Build binary (T, M_func) active mask from the functional table.

    For each (cell_type, marker) pair, the mask is 1 if the marker is listed
    as active for that cell type in the functional table, 0 otherwise.

    Args:
        functional_markers: List of marker names of length M_func.
        cell_types: List of cell type names of length T.
        functional_table: Marker-type association table. Defaults to
            DEFAULT_FUNCTIONAL_TABLE.

    Returns:
        Binary numpy array of shape (T, M_func).
    """
    if functional_table is None:
        functional_table = DEFAULT_FUNCTIONAL_TABLE

    T = len(cell_types)
    M = len(functional_markers)
    mask = np.zeros((T, M), dtype=np.int32)

    type_to_idx = {ct: i for i, ct in enumerate(cell_types)}

    for m_idx, marker in enumerate(functional_markers):
        if marker not in functional_table:
            logger.warning("Marker '%s' not in functional_table; skipping.", marker)
            continue
        for ct in functional_table[marker].get("active_types", []):
            if ct in type_to_idx:
                mask[type_to_idx[ct], m_idx] = 1

    return mask


# ---------------------------------------------------------------------------
# learn_functional_emissions
# ---------------------------------------------------------------------------

def learn_functional_emissions(
    observed: np.ndarray,
    proportions: np.ndarray,
    active_mask: np.ndarray,
    size_factors: np.ndarray,
    max_iter: int = 200,
    lr: float = 0.01,
    early_stopping_patience: int = 20,
    holdout_fraction: float = 0.1,
    lambda_sigma: float = 2.0,
    device: str = "cpu",
    seed: int = 42,
) -> Dict:
    """Learn per-type NB emission rates for functional markers via PyTorch.

    Fixed quantities (NOT optimized):
        - size_factors (s_i): sequencing depth
        - proportions (p_it): QP results from Module 3

    Optimized quantities:
        - log_lambda[k]: one parameter per active (type, marker) pair
        - log_r[m]: per-marker NB dispersion
        - log_b[m]: per-marker background rate

    Args:
        observed: (N, M_func) raw integer count matrix.
        proportions: (N, T) cell-type proportions (frozen from QP).
        active_mask: (T, M_func) binary mask from build_active_mask().
        size_factors: (N,) fixed sequencing size factors.
        max_iter: Maximum optimization iterations.
        lr: Adam learning rate.
        early_stopping_patience: Stop if holdout NLL does not improve for
            this many consecutive iterations.
        holdout_fraction: Fraction of spots withheld for early stopping.
        lambda_sigma: Log-normal prior sigma for lambda parameters.
        device: PyTorch device string.
        seed: Random seed for holdout split reproducibility.

    Returns:
        Dict with keys:
            "lambda": (T, M_func) emission rate matrix (zeros for inactive pairs)
            "dispersion": (M_func,) per-marker NB dispersion
            "background": (M_func,) per-marker background
            "train_nll": list of training NLL values per iteration
            "holdout_nll": list of holdout NLL values per iteration
            "skipped_pairs": list of (type_idx, marker_idx) tuples with
                insufficient support (< min_spots where p > 0.05)
    """
    import torch  # pylint: disable=import-outside-toplevel

    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    N, M = observed.shape
    T = proportions.shape[1]

    assert active_mask.shape == (T, M), (
        f"active_mask shape {active_mask.shape} does not match (T={T}, M={M})"
    )
    assert size_factors.shape == (N,), (
        f"size_factors shape {size_factors.shape} does not match N={N}"
    )

    MIN_SUPPORT_SPOTS = 20
    MIN_DOMINANT_SPOTS = 20
    MIN_PROPORTION_THRESHOLD = 0.05

    dev = torch.device(device)

    # --- Identify active pairs and filter for support ---
    active_pairs = []      # list of (t_idx, m_idx) with sufficient support
    skipped_pairs = []     # list of (t_idx, m_idx) with insufficient support

    for t_idx in range(T):
        for m_idx in range(M):
            if active_mask[t_idx, m_idx] == 0:
                continue
            support = np.sum(proportions[:, t_idx] > MIN_PROPORTION_THRESHOLD)
            if support < MIN_SUPPORT_SPOTS:
                skipped_pairs.append((t_idx, m_idx))
                logger.debug(
                    "Skipping pair (t=%d, m=%d): only %d spots with p > %.2f",
                    t_idx, m_idx, support, MIN_PROPORTION_THRESHOLD,
                )
            else:
                active_pairs.append((t_idx, m_idx))

    n_active = len(active_pairs)
    logger.info(
        "learn_functional_emissions: %d active pairs, %d skipped (insufficient support)",
        n_active, len(skipped_pairs),
    )

    if n_active == 0:
        logger.warning("No active pairs with sufficient support; returning zeros.")
        return {
            "lambda": np.zeros((T, M), dtype=np.float32),
            "dispersion": np.ones(M, dtype=np.float32),
            "background": np.zeros(M, dtype=np.float32),
            "train_nll": [],
            "holdout_nll": [],
            "skipped_pairs": skipped_pairs,
        }

    # --- Type-specific initialization for lambda ---
    DOMINANT_THRESHOLD = 0.3
    lam_init = np.zeros(n_active, dtype=np.float64)

    for k, (t_idx, m_idx) in enumerate(active_pairs):
        dominant_mask = proportions[:, t_idx] > DOMINANT_THRESHOLD
        n_dominant = dominant_mask.sum()

        if n_dominant >= MIN_DOMINANT_SPOTS:
            s_dominant = size_factors[dominant_mask]
            p_dominant = proportions[dominant_mask, t_idx]
            obs_dominant = observed[dominant_mask, m_idx].astype(np.float64)
            # lambda_init = median(S[dominant,m] / s_i) / median(p[dominant,t])
            rate_dominant = obs_dominant / np.maximum(s_dominant, 1e-8)
            lam_init[k] = np.median(rate_dominant) / (np.median(p_dominant) + 1e-8)
        else:
            # Fallback: use all spots
            # lambda_init = median(S[:,m] / s_i) / mean(sum_t p[:,t] * active[t,m])
            rate_all = observed[:, m_idx].astype(np.float64) / np.maximum(size_factors, 1e-8)
            active_cols = [t2 for t2, m2 in active_pairs if m2 == m_idx]
            if len(active_cols) > 0:
                p_active_sum = proportions[:, active_cols].sum(axis=1)
            else:
                p_active_sum = proportions[:, t_idx]
            denom = np.mean(p_active_sum) + 1e-8
            lam_init[k] = np.median(rate_all) / denom

        # Clamp to reasonable range
        lam_init[k] = np.clip(lam_init[k], 1e-4, 1e4)

    # --- Holdout split ---
    n_holdout = max(1, int(N * holdout_fraction))
    holdout_idx = rng.choice(N, size=n_holdout, replace=False)
    train_mask = np.ones(N, dtype=bool)
    train_mask[holdout_idx] = False

    # --- Convert to torch tensors ---
    S = torch.tensor(observed.astype(np.float32), device=dev)
    p = torch.tensor(proportions.astype(np.float32), device=dev)
    s = torch.tensor(size_factors.astype(np.float32), device=dev)  # fixed, no grad

    train_mask_t = torch.tensor(train_mask, device=dev)
    holdout_mask_t = ~train_mask_t

    # Extract active pair indices as tensors
    t_idx_list = [t for t, m in active_pairs]
    m_idx_list = [m for t, m in active_pairs]
    t_active = torch.tensor(t_idx_list, dtype=torch.long, device=dev)
    m_active = torch.tensor(m_idx_list, dtype=torch.long, device=dev)

    # --- Parameters ---
    log_lam = torch.tensor(
        np.log(lam_init + 1e-8).astype(np.float32),
        device=dev,
        requires_grad=True,
    )
    log_r = torch.zeros(M, dtype=torch.float32, device=dev, requires_grad=True)
    log_b = torch.full(
        (M,),
        np.log(0.01),
        dtype=torch.float32,
        device=dev,
        requires_grad=True,
    )

    optimizer = torch.optim.Adam([log_lam, log_r, log_b], lr=lr)

    def _compute_nll_functional(spot_mask_t):
        """Compute NB NLL on spots selected by spot_mask_t."""
        lam_active_vals = torch.exp(log_lam)  # (n_active,)
        r = torch.exp(log_r)                  # (M,)
        b = torch.exp(log_b)                  # (M,)

        # Build full (T, M) lambda matrix (inactive pairs = 0)
        lam_full = torch.zeros(T, M, dtype=torch.float32, device=dev)
        lam_full[t_active, m_active] = lam_active_vals

        # mu[i, m] = s_i * (sum_t p[i,t] * lam[t,m] + b[m])
        # p @ lam_full -> (N, M)
        expected_rate = p @ lam_full  # (N, M)
        mu = s.unsqueeze(1) * (expected_rate + b.unsqueeze(0))  # (N, M)

        # NB log-pmf: lgamma(k+r) - lgamma(r) - lgamma(k+1)
        #              + r*log(r/(r+mu)) + k*log(mu/(r+mu))
        S_sel = S[spot_mask_t]
        mu_sel = torch.clamp(mu[spot_mask_t], min=1e-8)
        r_bc = r.unsqueeze(0).clamp(min=1e-3)  # broadcast to (n_sel, M)

        logp = (
            torch.lgamma(S_sel + r_bc)
            - torch.lgamma(r_bc)
            - torch.lgamma(S_sel + 1.0)
            + r_bc * torch.log(r_bc / (r_bc + mu_sel))
            + S_sel * torch.log(mu_sel / (r_bc + mu_sel))
        )
        nll = -logp.mean()
        return nll

    # --- Optimization loop ---
    train_nll_history = []
    holdout_nll_history = []
    best_holdout_nll = float("inf")
    patience_counter = 0
    best_lam = log_lam.detach().clone()
    best_r = log_r.detach().clone()
    best_b = log_b.detach().clone()

    lam_init_t = torch.tensor(lam_init.astype(np.float32), device=dev)

    for iteration in range(max_iter):
        optimizer.zero_grad()

        train_nll = _compute_nll_functional(train_mask_t)

        # Log-normal prior on lambda: log(lambda) ~ N(log(lam_init), sigma^2)
        lam_prior = (
            (log_lam - torch.log(lam_init_t + 1e-8)) ** 2
        ).sum() / (2.0 * lambda_sigma ** 2)

        loss = train_nll + lam_prior / max(train_mask.sum(), 1)
        loss.backward()
        optimizer.step()

        train_nll_history.append(float(train_nll.item()))

        with torch.no_grad():
            h_nll = _compute_nll_functional(holdout_mask_t).item()
        holdout_nll_history.append(h_nll)

        if h_nll < best_holdout_nll - 1e-5:
            best_holdout_nll = h_nll
            patience_counter = 0
            best_lam = log_lam.detach().clone()
            best_r = log_r.detach().clone()
            best_b = log_b.detach().clone()
        else:
            patience_counter += 1
            if patience_counter >= early_stopping_patience:
                logger.info(
                    "Early stopping at iteration %d (patience=%d).",
                    iteration, early_stopping_patience,
                )
                break

    # --- Reconstruct full (T, M) lambda matrix from best params ---
    lam_vals = torch.exp(best_lam).cpu().numpy()
    lam_full_np = np.zeros((T, M), dtype=np.float32)
    for k, (t_idx, m_idx) in enumerate(active_pairs):
        lam_full_np[t_idx, m_idx] = lam_vals[k]

    dispersion = torch.exp(best_r).cpu().numpy().astype(np.float32)
    background = torch.exp(best_b).cpu().numpy().astype(np.float32)

    return {
        "lambda": lam_full_np,
        "dispersion": dispersion,
        "background": background,
        "train_nll": train_nll_history,
        "holdout_nll": holdout_nll_history,
        "skipped_pairs": skipped_pairs,
    }


# ---------------------------------------------------------------------------
# gate_functional_markers
# ---------------------------------------------------------------------------

def gate_functional_markers(
    observed: np.ndarray,
    proportions: np.ndarray,
    lam: np.ndarray,
    background: np.ndarray,
    size_factors: np.ndarray,
    active_mask: np.ndarray,
    cell_types: List[str],
    functional_markers: List[str],
    gmm_min_proportion: float = 0.05,
    min_spots: int = 20,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
    """Gate spots as functionally active using observed-to-expected ratio GMM.

    For each active (type, marker) pair, computes:
        expected[i,m] = s_i * (sum_t p[i,t] * lambda[t,m] + b_m)
        ratio[i,m]    = S[i,m] / expected[i,m]

    Fits a 2-component GMM on log(ratio) for spots where p[i,t] > gmm_min_proportion.
    A spot is positive if:
        - GMM high-component posterior > 0.5  AND  ratio > 1.0

    Fallback (unimodal fit, component separation < 0.5 * pooled sigma):
        - Threshold at ratio > 1.0 only.

    Args:
        observed: (N, M_func) raw count matrix (same as NB input).
        proportions: (N, T) cell-type proportions.
        lam: (T, M_func) emission rate matrix from learn_functional_emissions().
        background: (M_func,) per-marker background from learn_functional_emissions().
        size_factors: (N,) sequencing depth.
        active_mask: (T, M_func) binary mask.
        cell_types: List of length T.
        functional_markers: List of length M_func.
        gmm_min_proportion: Minimum proportion threshold to include a spot in
            GMM fitting for a given type.
        min_spots: Minimum qualifying spots per pair to attempt GMM.

    Returns:
        Tuple of:
            intensity_df: (N, n_pairs) DataFrame with observed-to-expected ratios.
                Columns: "{type}_{marker}_ratio". Rows: spot indices.
            gates_df: (N, n_pairs) binary DataFrame.
                Columns: "{type}_{marker}_gate". Rows: spot indices.
            summary: Dict with per-pair gating statistics.
    """
    from sklearn.mixture import GaussianMixture  # pylint: disable=import-outside-toplevel

    N, M = observed.shape
    T = proportions.shape[1]

    # Compute expected counts: mu[i,m] = s_i * (p @ lam + b)[m]
    expected = size_factors[:, None] * (proportions @ lam + background[None, :])  # (N, M)
    expected = np.maximum(expected, 1e-8)

    ratio_cols = {}
    gate_cols = {}
    summary = {}

    for t_idx, ct in enumerate(cell_types):
        for m_idx, marker in enumerate(functional_markers):
            if active_mask[t_idx, m_idx] == 0:
                continue

            col_key = f"{ct}_{marker}"
            qualifying = proportions[:, t_idx] > gmm_min_proportion
            n_qual = qualifying.sum()

            ratio_all = observed[:, m_idx].astype(np.float64) / expected[:, m_idx]
            ratio_cols[f"{col_key}_ratio"] = ratio_all

            if n_qual < min_spots:
                # Insufficient data: fallback to ratio > 1
                gates = (ratio_all > 1.0).astype(np.float32)
                gate_cols[f"{col_key}_gate"] = gates
                summary[(ct, marker)] = {
                    "n_qualifying": int(n_qual),
                    "gating_method": "ratio_fallback_insufficient",
                    "n_positive": int(gates.sum()),
                    "fraction_positive": float(gates.mean()),
                }
                continue

            ratio_qual = ratio_all[qualifying]
            log_ratio_qual = np.log(np.maximum(ratio_qual, 1e-8)).reshape(-1, 1)

            gating_method = "ratio_fallback_unimodal"
            positive_mask = np.zeros(N, dtype=bool)

            try:
                gmm = GaussianMixture(
                    n_components=2,
                    covariance_type="full",
                    random_state=42,
                    n_init=3,
                    max_iter=100,
                )
                gmm.fit(log_ratio_qual)

                means = gmm.means_.flatten()
                stds = np.sqrt(gmm.covariances_.flatten())
                pooled_std = np.sqrt(
                    gmm.weights_[0] * stds[0] ** 2
                    + gmm.weights_[1] * stds[1] ** 2
                )
                component_separation = abs(means[0] - means[1])

                if component_separation < 0.5 * pooled_std:
                    # Unimodal distribution: use simple ratio > 1 threshold
                    gating_method = "ratio_fallback_unimodal"
                    qual_pos = qualifying & (ratio_all > 1.0)
                    positive_mask = qual_pos
                else:
                    # Bimodal: use GMM posteriors
                    # High component = higher mean
                    high_comp = int(np.argmax(means))
                    # Get posteriors for ALL spots (not just qualifying)
                    log_ratio_all = np.log(np.maximum(ratio_all, 1e-8)).reshape(-1, 1)
                    posteriors = gmm.predict_proba(log_ratio_all)  # (N, 2)
                    high_posterior = posteriors[:, high_comp]

                    positive_mask = (
                        qualifying
                        & (high_posterior > 0.5)
                        & (ratio_all > 1.0)
                    )
                    gating_method = "gmm"

            except Exception as exc:  # pylint: disable=broad-except
                logger.warning(
                    "GMM failed for (%s, %s): %s; using ratio > 1 fallback.",
                    ct, marker, exc,
                )
                gating_method = "ratio_fallback_exception"
                positive_mask = qualifying & (ratio_all > 1.0)

            gates = positive_mask.astype(np.float32)
            gate_cols[f"{col_key}_gate"] = gates
            summary[(ct, marker)] = {
                "n_qualifying": int(n_qual),
                "gating_method": gating_method,
                "n_positive": int(gates.sum()),
                "fraction_positive": float(gates.mean()),
            }

    intensity_df = pd.DataFrame(ratio_cols)
    gates_df = pd.DataFrame(gate_cols)

    return intensity_df, gates_df, summary


# ---------------------------------------------------------------------------
# compute_spatial_statistics
# ---------------------------------------------------------------------------

def compute_spatial_statistics(
    gates_df: pd.DataFrame,
    spot_coords: np.ndarray,
    active_pairs: List[Tuple[str, str]],
) -> Dict[Tuple[str, str], Dict]:
    """Compute Moran's I spatial autocorrelation on gate calls.

    For each (type, marker) pair in active_pairs, computes univariate
    Moran's I on the binary gate column using a k=6 neighbor graph.

    Args:
        gates_df: (N, n_pairs) binary gate DataFrame from gate_functional_markers().
            Column names follow "{type}_{marker}_gate" convention.
        spot_coords: (N, 2) spatial coordinates.
        active_pairs: List of (cell_type, marker) tuples to compute statistics for.

    Returns:
        Dict mapping (cell_type, marker) -> {"morans_i": float, "morans_p": float}.
        NaN values indicate insufficient data or computation failure.
    """
    try:
        import anndata  # pylint: disable=import-outside-toplevel
        import scipy.sparse as sp  # pylint: disable=import-outside-toplevel
        import squidpy as sq  # pylint: disable=import-outside-toplevel
        _has_squidpy = True
    except ImportError:
        _has_squidpy = False
        logger.warning(
            "squidpy not available; spatial statistics will return NaN."
        )

    results = {}
    N = spot_coords.shape[0]

    # Build spatial weight matrix once (shared across all pairs)
    W_norm = None
    W_raw = None  # original (unnormalized) for S0, S1, S2 computations
    if _has_squidpy and N >= 10:
        try:
            dummy_adata = anndata.AnnData(X=np.zeros((N, 1), dtype=np.float32))
            dummy_adata.obsm["spatial"] = spot_coords.astype(float)
            sq.gr.spatial_neighbors(dummy_adata, coord_type="generic", n_neighs=6)
            W_raw = dummy_adata.obsp["spatial_connectivities"]
            row_sums = np.asarray(W_raw.sum(axis=1)).ravel()
            row_sums[row_sums == 0] = 1.0
            W_norm = sp.diags(1.0 / row_sums) @ W_raw
        except Exception as exc:  # pylint: disable=broad-except
            logger.warning("Failed to build spatial weight matrix: %s", exc)
            W_norm = None
            W_raw = None

    # Precompute Moran's I weight sums if we have a valid weight matrix
    moran_params = None
    if W_norm is not None and W_raw is not None:
        try:
            S0 = float(np.asarray(W_raw.sum()))
            S1 = 0.5 * float(np.asarray((W_raw + W_raw.T).power(2).sum()))
            row_sums_vec = np.asarray(W_raw.sum(axis=1)).ravel()
            col_sums_vec = np.asarray(W_raw.sum(axis=0)).ravel()
            S2 = float(np.sum((row_sums_vec + col_sums_vec) ** 2))
            moran_params = {"S0": S0, "S1": S1, "S2": S2}
        except Exception as exc:  # pylint: disable=broad-except
            logger.warning("Failed to precompute Moran weight sums: %s", exc)
            moran_params = None

    from scipy.stats import norm as scipy_norm  # pylint: disable=import-outside-toplevel

    for ct, marker in active_pairs:
        col_name = f"{ct}_{marker}_gate"
        nan_result = {"morans_i": float("nan"), "morans_p": float("nan")}

        if col_name not in gates_df.columns:
            logger.debug("Column '%s' not in gates_df; skipping.", col_name)
            results[(ct, marker)] = nan_result
            continue

        x = gates_df[col_name].values.astype(float)
        n_positive = int((x > 0.5).sum())

        if n_positive < 5 or (N - n_positive) < 5:
            logger.debug(
                "Pair (%s, %s): too few positives (%d) or negatives for Moran's I.",
                ct, marker, n_positive,
            )
            results[(ct, marker)] = nan_result
            continue

        if W_norm is None or moran_params is None:
            results[(ct, marker)] = nan_result
            continue

        try:
            S0 = moran_params["S0"]
            S1 = moran_params["S1"]
            S2 = moran_params["S2"]

            if S0 == 0:
                results[(ct, marker)] = nan_result
                continue

            xz = (x - x.mean()) / (x.std() + 1e-12)
            Wxz = np.asarray(W_norm.dot(xz)).ravel()
            # Univariate Moran's I = (N / S0) * (z' W_norm z)
            I_val = float(N / S0 * xz.dot(Wxz))

            # Analytical variance under normality assumption
            # (Cliff & Ord 1981, standard formula)
            E_I = -1.0 / max(N - 1, 1)
            n = N
            b2 = float(np.mean(xz ** 4)) / (float(np.mean(xz ** 2)) ** 2 + 1e-12)
            var_I_numerator = (
                n * ((n ** 2 - 3 * n + 3) * S1 - n * S2 + 3 * S0 ** 2)
                - b2 * ((n ** 2 - n) * S1 - 2 * n * S2 + 6 * S0 ** 2)
            )
            var_I_denominator = (n - 1) * (n - 2) * (n - 3) * S0 ** 2
            if var_I_denominator > 0:
                var_I = var_I_numerator / var_I_denominator - E_I ** 2
            else:
                var_I = 1.0
            var_I = max(var_I, 1e-12)

            z_score = (I_val - E_I) / np.sqrt(var_I)
            p_val = float(2.0 * scipy_norm.sf(abs(z_score)))

            results[(ct, marker)] = {"morans_i": I_val, "morans_p": p_val}

        except Exception as exc:  # pylint: disable=broad-except
            logger.warning(
                "Moran's I computation failed for (%s, %s): %s",
                ct, marker, exc,
            )
            results[(ct, marker)] = nan_result

    return results


# ---------------------------------------------------------------------------
# gmm_gate_cells
# ---------------------------------------------------------------------------

def gmm_gate_cells(
    cell_protein: np.ndarray,
    cell_types: np.ndarray,
    type_names: List[str],
    marker_names: List[str],
    active_mask: np.ndarray,
    min_cells: int = 20,
    bimodality_threshold: float = 1.5,
    posterior_threshold: float = 0.5,
    min_high_component_log_mean: float = 1.0,
) -> Tuple[pd.DataFrame, Dict]:
    """Gate per-cell deconvolved protein levels via 2-component GMM.

    For each active (type, marker) pair, pools all cells of that type and fits
    a 2-component GMM. Cells in the higher component are called positive.

    Mirrors the spot-level gating pattern in gate_functional_markers() but
    operates on per-cell values from SACE protein allocation.

    Three-stage gating:
    1. Bimodality check: separation >= bimodality_threshold * pooled_std.
    2. Absolute signal check: high component log-mean >= min_high_component_log_mean.
       This suppresses spurious bimodality from SACE averaging leaking near-zero
       signal to negative cells (e.g., sparse markers like PD-L1 on Epithelial
       cells where high-comp log-mean ≈ -0.5, raw ≈ 0.6 vs well-calibrated pairs
       like T cell PD-1 where high-comp log-mean ≈ 3.0, raw ≈ 20).
    3. Posterior threshold: cell assigned to high component only if posterior > threshold.

    Args:
        cell_protein: (N_cells, M) deconvolved protein counts per cell.
        cell_types: (N_cells,) string array of assigned cell types.
        type_names: Ordered list of T type names.
        marker_names: Ordered list of M marker names.
        active_mask: (T, M) binary mask from build_active_mask().
        min_cells: Minimum cells of the active type to attempt GMM.
        bimodality_threshold: GMM fires only when component separation
            >= bimodality_threshold * pooled_std. Raise to suppress
            over-calling for sparse markers (default 1.5).
        posterior_threshold: Minimum posterior probability for the high
            component to call a cell positive (default 0.5).
        min_high_component_log_mean: Minimum log-mean for the high GMM
            component to call any cells positive. Prevents near-zero SACE
            allocation leakage from being called positive. Default 1.0
            (≈ e^1.0 ≈ 2.7 in raw units), calibrated on Xenium RCC Region 0
            where bad pairs had high-comp log-mean ≤ -0.49 and good pairs ≥ 2.96.

    Returns:
        Tuple of:
            gates_df: (N_cells, n_active_pairs) DataFrame.
                Columns: "func_{marker}_{type}_gate". Values: 0 or 1.
                Cells whose type doesn't match the pair get 0.
            summary: Dict mapping (type_name, marker_name) -> {
                n_cells, gating_method, n_positive, fraction_positive,
                bimodality_threshold, posterior_threshold,
                min_high_component_log_mean, high_component_log_mean}.
    """
    from sklearn.mixture import GaussianMixture  # pylint: disable=import-outside-toplevel

    N = cell_protein.shape[0]
    gate_cols = {}
    summary = {}

    for t_idx, ct in enumerate(type_names):
        for m_idx, marker in enumerate(marker_names):
            if active_mask[t_idx, m_idx] == 0:
                continue

            col_name = f"func_{marker}_{ct}_gate"
            gates = np.zeros(N, dtype=np.int32)
            type_mask = cell_types == ct
            n_type = int(type_mask.sum())

            if n_type < min_cells:
                gate_cols[col_name] = gates
                summary[(ct, marker)] = {
                    "n_cells": n_type,
                    "gating_method": "insufficient_cells",
                    "n_positive": 0,
                    "fraction_positive": 0.0,
                }
                continue

            values = cell_protein[type_mask, m_idx]
            n_nan = int(np.isnan(values).sum())
            if n_nan > 0:
                logger.warning(
                    "gmm_gate_cells: %d NaN values in cell_protein for (%s, %s); "
                    "replacing with 0 before GMM fit.", n_nan, ct, marker,
                )
                values = np.where(np.isnan(values), 0.0, values)
            log_values = np.log(np.maximum(values, 1e-8)).reshape(-1, 1)

            gating_method = "threshold_fallback_unimodal"
            high_log_mean = None

            try:
                gmm = GaussianMixture(
                    n_components=2,
                    covariance_type="full",
                    random_state=42,
                    n_init=3,
                    max_iter=100,
                )
                gmm.fit(log_values)

                means = gmm.means_.flatten()
                stds = np.sqrt(gmm.covariances_.flatten())
                pooled_std = np.sqrt(
                    gmm.weights_[0] * stds[0] ** 2
                    + gmm.weights_[1] * stds[1] ** 2
                )
                separation = abs(means[0] - means[1])

                if separation < bimodality_threshold * pooled_std:
                    gating_method = "threshold_fallback_unimodal"
                else:
                    high_comp = int(np.argmax(means))
                    high_log_mean = float(means[high_comp])
                    if high_log_mean < min_high_component_log_mean:
                        # High component has near-zero absolute signal — SACE
                        # averaging leaked small amounts to negative cells, creating
                        # a spurious bimodal distribution. Suppress the gate.
                        gating_method = "threshold_fallback_weak_signal"
                    else:
                        posteriors = gmm.predict_proba(log_values)
                        positive = posteriors[:, high_comp] > posterior_threshold
                        gates[type_mask] = positive.astype(np.int32)
                        gating_method = "gmm_bimodal"

            except Exception as exc:  # pylint: disable=broad-except
                logger.warning(
                    "GMM failed for (%s, %s): %s. Using fallback.", ct, marker, exc
                )
                gating_method = "threshold_fallback_error"

            gate_cols[col_name] = gates
            n_pos = int(gates[type_mask].sum()) if type_mask.any() else 0
            summary[(ct, marker)] = {
                "n_cells": n_type,
                "gating_method": gating_method,
                "n_positive": n_pos,
                "fraction_positive": n_pos / n_type if n_type > 0 else 0.0,
                "bimodality_threshold": bimodality_threshold,
                "posterior_threshold": posterior_threshold,
                "min_high_component_log_mean": min_high_component_log_mean,
                "high_component_log_mean": high_log_mean,
            }

    gates_df = pd.DataFrame(gate_cols)
    return gates_df, summary
