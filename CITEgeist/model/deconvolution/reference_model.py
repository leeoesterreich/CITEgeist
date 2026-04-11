"""NB reference model for proportion refinement using scRNA-seq marker genes.

Trains NB emission profiles from annotated scRNA-seq, then refines
cuOPT QP proportions via MAP inference with Dirichlet protein prior.
"""

import logging
import pickle
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)


@dataclass
class ReferenceProfile:
    """Frozen NB reference profiles learned from scRNA-seq.

    Attributes:
        mu: (T, G_marker) mean expression per type per marker gene (raw count scale).
        alpha: (G_marker,) per-gene NB overdispersion. Var(Y) = mu + alpha * mu^2.
        gene_names: Marker gene names (length G_marker).
        type_names: Cell type names (length T).
        de_results: Full DE results from scanpy.tl.rank_genes_groups.
        n_cells_per_type: Number of cells per type in the reference.
    """

    mu: np.ndarray
    alpha: np.ndarray
    gene_names: List[str]
    type_names: List[str]
    de_results: pd.DataFrame
    n_cells_per_type: Dict[str, int]

    def save(self, path: str) -> None:
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path: str) -> "ReferenceProfile":
        with open(path, "rb") as f:
            return pickle.load(f)


def train_reference(
    reference_adata: sc.AnnData,
    *,
    cell_type_col: str = "cell_type",
    spatial_genes: Optional[List[str]] = None,
    n_markers_per_type: int = 20,
    de_method: str = "wilcoxon",
    type_mapping: Optional[Dict[str, str]] = None,
) -> ReferenceProfile:
    """Train NB reference profiles from annotated scRNA-seq.

    Args:
        reference_adata: Annotated scRNA-seq AnnData with raw counts in .X.
        cell_type_col: Column in .obs with cell type annotations.
        spatial_genes: Gene names present in spatial data. If provided,
            reference is filtered to this gene universe before DE.
        n_markers_per_type: Number of top DE marker genes per type.
        de_method: scanpy DE method ("wilcoxon" or "t-test").
        type_mapping: Optional dict mapping reference type names to
            CITEgeist type names. Multiple ref types can map to one
            CITEgeist type (cells are pooled).

    Returns:
        Frozen ReferenceProfile with NB emission parameters.
    """
    adata = reference_adata.copy()

    # --- Apply type mapping (pool subtypes) ---
    if type_mapping is not None:
        raw_labels = adata.obs[cell_type_col].values.copy()
        mapped_labels = np.array([type_mapping.get(str(lbl), str(lbl)) for lbl in raw_labels])
        adata.obs[cell_type_col] = mapped_labels

    # --- Filter to spatial gene universe ---
    if spatial_genes is not None:
        shared = [g for g in adata.var_names if g in set(spatial_genes)]
        if len(shared) == 0:
            raise ValueError("No shared genes between reference and spatial data.")
        logger.info(
            "Filtered to %d shared genes (ref=%d, spatial=%d)",
            len(shared),
            adata.n_vars,
            len(spatial_genes),
        )
        adata = adata[:, shared].copy()

    # --- Differential expression ---
    # Store raw counts for mu estimation, but run DE on log-normalized data
    raw_X = adata.X.copy() if not hasattr(adata.X, "toarray") else adata.X.toarray().copy()
    adata_de = adata.copy()
    sc.pp.normalize_total(adata_de, target_sum=1e4)
    sc.pp.log1p(adata_de)
    # Request 3x markers to allow deduplication headroom
    sc.tl.rank_genes_groups(adata_de, groupby=cell_type_col, method=de_method, n_genes=n_markers_per_type * 3)
    de_results = sc.get.rank_genes_groups_df(adata_de, group=None)

    # --- Select top marker genes (deduplicated) ---
    type_names_sorted = sorted(adata.obs[cell_type_col].unique())
    selected_genes = []
    seen_genes = set()
    for tname in type_names_sorted:
        type_de = de_results[de_results["group"] == tname].sort_values("scores", ascending=False)
        count = 0
        for _, row in type_de.iterrows():
            if count >= n_markers_per_type:
                break
            gene = row["names"]
            if gene not in seen_genes:
                selected_genes.append(gene)
                seen_genes.add(gene)
                count += 1

    if len(selected_genes) == 0:
        raise ValueError("No marker genes selected after DE filtering.")

    logger.info(
        "Selected %d marker genes across %d types",
        len(selected_genes),
        len(type_names_sorted),
    )

    # --- Compute mu (raw count scale) and alpha (overdispersion) ---
    # Debug: verify all selected genes are in adata
    adata_genes = set(adata.var_names)
    missing_genes = [g for g in selected_genes if g not in adata_genes]
    if missing_genes:
        logger.warning(
            "%d/%d selected marker genes NOT in filtered adata (removing): %s",
            len(missing_genes),
            len(selected_genes),
            missing_genes[:5],
        )
        selected_genes = [g for g in selected_genes if g in adata_genes]
    logger.info("Subsetting adata (%d genes) to %d marker genes", adata.n_vars, len(selected_genes))
    marker_adata = adata[:, selected_genes]
    # Use raw counts for mu/alpha estimation (NOT the log-normalized data)
    marker_idx_in_adata = [list(adata.var_names).index(g) for g in selected_genes]
    if hasattr(raw_X, "toarray"):
        raw_X = raw_X.toarray() if not isinstance(raw_X, np.ndarray) else raw_X
    X = np.asarray(raw_X, dtype=np.float64)[:, marker_idx_in_adata]

    n_types = len(type_names_sorted)
    n_markers = len(selected_genes)
    mu = np.zeros((n_types, n_markers), dtype=np.float64)
    alpha_per_type = np.zeros((n_types, n_markers), dtype=np.float64)
    n_cells_per_type = {}

    for i, tname in enumerate(type_names_sorted):
        mask = (marker_adata.obs[cell_type_col] == tname).values
        X_type = X[mask]
        n_cells_per_type[tname] = int(mask.sum())

        mu[i] = X_type.mean(axis=0)
        var_type = X_type.var(axis=0)

        # Method of moments on raw counts: alpha = (var - mean) / mean^2
        mean_safe = np.maximum(mu[i], 1e-8)
        alpha_raw = (var_type - mu[i]) / (mean_safe**2)
        alpha_per_type[i] = np.maximum(alpha_raw, 0.01)

    # Pool alpha across types via median (robust to outlier types)
    alpha = np.median(alpha_per_type, axis=0)
    alpha = np.maximum(alpha, 0.01)

    return ReferenceProfile(
        mu=mu,
        alpha=alpha,
        gene_names=selected_genes,
        type_names=type_names_sorted,
        de_results=de_results,
        n_cells_per_type=n_cells_per_type,
    )


def refine_proportions_nb(
    spot_counts: np.ndarray,
    pi_protein: np.ndarray,
    reference: ReferenceProfile,
    *,
    kappa: float = 10.0,
    epsilon: float = 1.0,
    lr: float = 0.1,
    max_iter: int = 200,
    convergence_tol: float = 1e-4,
    device: str = "cpu",
) -> np.ndarray:
    """Refine protein QP proportions using NB MAP inference on marker genes.

    Args:
        spot_counts: (N, G_marker) observed counts for marker genes.
        pi_protein: (N, T) protein QP proportions.
        reference: Frozen ReferenceProfile with mu and alpha.
        kappa: Dirichlet prior concentration. Higher = more protein trust.
        epsilon: Floor for Dirichlet concentration (prevents zero issues).
        lr: Adam learning rate.
        max_iter: Maximum optimization iterations.
        convergence_tol: Early stop if max|delta_pi| < tol.
        device: "cpu" or "cuda".

    Returns:
        (N, T) refined proportions on the simplex.
    """
    import torch  # pylint: disable=import-outside-toplevel

    _, _ = spot_counts.shape
    _ = pi_protein.shape[1]

    Y = torch.tensor(spot_counts, dtype=torch.float64, device=device)
    mu = torch.tensor(reference.mu, dtype=torch.float64, device=device)
    alpha = torch.tensor(reference.alpha, dtype=torch.float64, device=device)
    pi_prot = torch.tensor(pi_protein, dtype=torch.float64, device=device)

    # Dirichlet concentration: alpha_t = epsilon + kappa * pi_protein_t
    # epsilon=1.0 ensures all concentrations >= 1 (uniform-like, not corner-seeking)
    dir_alpha = epsilon + kappa * pi_prot

    # Precompute observed marker totals (fixed)
    observed_marker_total = Y.sum(dim=1, keepdim=True)  # (N, 1)

    # Initialize logits from protein proportions
    z = torch.log(pi_prot + 0.01).clone().detach().requires_grad_(True)

    optimizer = torch.optim.Adam([z], lr=lr)

    # NB parameters: total_count = 1/alpha
    inv_alpha = 1.0 / alpha

    prev_pi = torch.softmax(z, dim=1).detach()

    for iteration in range(max_iter):
        optimizer.zero_grad()

        pi = torch.softmax(z, dim=1)

        # Update library size from CURRENT pi each iteration (not frozen)
        expected_total = (pi @ mu).sum(dim=1, keepdim=True)  # (N, 1)
        l_s = observed_marker_total / (expected_total.detach() + 1e-8)

        # Expected counts: l_s * sum_t(pi_st * mu_tg)
        expected = l_s * (pi @ mu)
        expected = torch.clamp(expected, min=1e-8)

        # NB log-likelihood (mean-dispersion parameterization)
        log_prob = (
            torch.lgamma(Y + inv_alpha)
            - torch.lgamma(inv_alpha)
            - torch.lgamma(Y + 1.0)
            + inv_alpha * torch.log(1.0 / (1.0 + alpha * expected))
            + Y * torch.log(alpha * expected / (1.0 + alpha * expected) + 1e-30)
        )
        nb_ll = log_prob.sum()

        # Dirichlet log-prior: sum_t (alpha_t - 1) * log(pi_t)
        dir_lp = ((dir_alpha - 1.0) * torch.log(pi + 1e-30)).sum()

        # MAP objective (maximize)
        loss = -(nb_ll + dir_lp)
        loss.backward()
        optimizer.step()

        # Convergence check
        with torch.no_grad():
            curr_pi = torch.softmax(z, dim=1)
            max_delta = (curr_pi - prev_pi).abs().max().item()
            if max_delta < convergence_tol:
                logger.info(
                    "NB refinement converged at iteration %d (max_delta=%.2e)",
                    iteration,
                    max_delta,
                )
                break
            prev_pi = curr_pi.clone()

    with torch.no_grad():
        pi_refined = torch.softmax(z, dim=1).cpu().numpy()

    return pi_refined
