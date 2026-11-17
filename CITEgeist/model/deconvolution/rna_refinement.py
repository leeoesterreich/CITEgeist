"""Reference-free RNA refinement of Phase 1 (protein) cell-type proportions.

Model inputs are only protein-derived proportions ``p0``, raw spot RNA counts ``Y`` and
nuclei counts ``n``. No scRNA reference is used.

Estimand: an RNA-mass fraction exported as a cell fraction under an equal-RNA-yield
assumption (every cell type contributes the same amount of RNA per cell).

Algorithm (per pass): learn per-type expression profiles by weighted NNLS of per-cell
counts ``Y/n`` on the current proportions, then refit each spot's proportions by
multiplicative EM on a multinomial (or Poisson) mixture with a Dirichlet-style prior of
pseudo-counts ``c * L_s * p0``. The prior always stays on the original ``p0``.
"""
from __future__ import annotations

import re
from typing import Optional

import numpy as np
from scipy.optimize import nnls

_EXCLUDE_GENE = re.compile(r"^MT-|^RP[LS]|^HB[AB]|^MALAT1$")
_MIN_GENES = 50
_HIGH_COSINE = 0.95
_MAX_COND = 1e8


def validate_raw_counts(Y: np.ndarray) -> None:
    """Raise ValueError unless Y is finite, non-negative and integer-like raw counts."""
    Y = np.asarray(Y)
    if Y.size == 0:
        raise ValueError("Expected raw counts, got an empty matrix")
    if not np.isfinite(Y).all():
        raise ValueError("Expected raw counts: matrix has non-finite values")
    if (Y < 0).any():
        raise ValueError("Expected raw counts: matrix has negative values")
    if np.abs(Y - np.round(Y)).max() > 1e-6:
        raise ValueError(
            "Expected raw counts (integers); got non-integer values, e.g. "
            "log-normalized or scaled data. Pass the raw count layer."
        )


def select_genes(
    Y: np.ndarray,
    genes: list,
    *,
    n_top: int = 2000,
    min_detect_frac: float = 0.05,
) -> np.ndarray:
    """Indices of genes used for refinement.

    Drop zero-total genes, keep genes detected in >= ``min_detect_frac`` of spots,
    exclude MT-/RPL/RPS/HBA/HBB/MALAT1, then take the top ``n_top`` highly variable
    genes (scanpy seurat_v3 on raw counts). If <= ``n_top`` remain, all are kept.
    Raises ValueError if fewer than 50 genes remain.
    """
    Y = np.asarray(Y, dtype=float)
    if Y.shape[1] != len(genes):
        raise ValueError("genes length does not match Y columns")
    detect = (Y > 0).mean(axis=0)
    keep = (Y.sum(axis=0) > 0) & (detect >= min_detect_frac)
    keep &= np.array([_EXCLUDE_GENE.search(str(g)) is None for g in genes])
    idx = np.flatnonzero(keep)
    if idx.size < _MIN_GENES:
        raise ValueError(f"Only {idx.size} genes remain after filtering; need at least {_MIN_GENES} genes")
    if idx.size <= n_top:
        return idx
    import anndata
    import scanpy as sc

    ad = anndata.AnnData(Y[:, idx].astype(np.float32))
    sc.pp.highly_variable_genes(ad, n_top_genes=n_top, flavor="seurat_v3")
    return np.sort(idx[np.flatnonzero(ad.var["highly_variable"].to_numpy())])


def learn_profiles(
    Y: np.ndarray,
    n: np.ndarray,
    P: np.ndarray,
    *,
    anchor_thr: float = 0.0,
    non_anchor_weight: float = 0.1,
) -> np.ndarray:
    """Per-gene weighted NNLS of per-cell counts Y/n on proportions P. Returns (K, G), unnormalized.

    Weights are n_s, times ``non_anchor_weight`` for spots whose max proportion is below
    ``anchor_thr`` (only when ``anchor_thr`` > 0). Spots with n_s <= 0 are excluded.
    """
    ok = n > 0
    Y, n, P = Y[ok], n[ok], P[ok]
    w = n * np.where(P.max(1) >= anchor_thr, 1.0, non_anchor_weight) if anchor_thr > 0 else n.copy()
    sw = np.sqrt(w)[:, None]
    A = P * sw
    Z = (Y / n[:, None]) * sw
    mu = np.zeros((P.shape[1], Y.shape[1]))
    for g in range(Y.shape[1]):
        mu[:, g], _ = nnls(A, Z[:, g])
    return mu


def _floor_profiles(mu: np.ndarray) -> np.ndarray:
    return mu + 1e-6 * mu.max()


def em_refine(
    Y: np.ndarray,
    phi: np.ndarray,
    p0: np.ndarray,
    *,
    c: float,
    tol: float = 1e-5,
    max_iter: int = 150,
    model: str = "multinomial",
) -> tuple:
    """Multiplicative EM for per-spot proportions with prior pseudo-counts ``c * L_s * p0``.

    ``model="multinomial"``: ``phi`` rows are normalized profiles (sum to 1).
    ``model="poisson"``: ``phi`` is the unnormalized per-cell rate profile; the update is
    the Poisson-mixture EM with each type's expected count divided by its total rate, so
    yields are accounted for (equal yields reduce to the multinomial update).
    EM always starts from p0 (for every c, including c == 0), as in the spike.
    For the Poisson model with c > 0 the update is a heuristic fixed-point iteration (the
    prior is added in q-space to yield-scaled counts), not a guaranteed monotone EM.
    Returns (p, {"iters", "converged"}).
    """
    L = Y.sum(1, keepdims=True)
    alpha = c * L * p0
    q = p0.copy()
    q = np.clip(q, 1e-6, None)
    q /= q.sum(1, keepdims=True)
    scale = None
    if model == "poisson":
        tot = phi.sum(1)
        scale = tot.mean() / np.maximum(tot, 1e-300)
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        pi = np.clip(q @ phi, 1e-12, None)
        R = q * ((Y / pi) @ phi.T)
        if scale is not None:
            R = R * scale
        denom = R.sum(1, keepdims=True) if scale is not None else L
        qn = (R + alpha) / (denom + alpha.sum(1, keepdims=True))
        qn = np.clip(qn, 1e-9, None)
        qn /= qn.sum(1, keepdims=True)
        done = np.abs(qn - q).max() < tol
        q = qn
        if done:
            converged = True
            break
    return q, {"iters": it, "converged": converged}


def _profile_diagnostics(mu: np.ndarray) -> tuple:
    """rank, condition number (of row-normalized profiles), cosine matrix."""
    rn = np.linalg.norm(mu, axis=1, keepdims=True)
    mn = mu / np.maximum(rn, 1e-300)
    s = np.linalg.svd(mn, compute_uv=False)
    rank = int(np.linalg.matrix_rank(mn))
    cond = float(s.max() / s.min()) if s.min() > 0 else float("inf")
    return rank, cond, mn @ mn.T


def _bootstrap_cv(Y, n, P, n_bootstrap, anchor_thr) -> np.ndarray:
    rng = np.random.default_rng(0)
    S = Y.shape[0]
    mus = np.stack(
        [
            learn_profiles(Y[i], n[i], P[i], anchor_thr=anchor_thr)
            for i in (rng.integers(0, S, S) for _ in range(n_bootstrap))
        ]
    )
    mean, sd = mus.mean(0), mus.std(0)
    cv = np.where(mean > 0, sd / np.maximum(mean, 1e-300), np.nan)
    return np.array([np.nanmedian(r) if np.isfinite(r).any() else np.nan for r in cv])


class _DropType(Exception):
    def __init__(self, k):
        self.k = k


def refine_proportions(
    Y: np.ndarray,
    n: np.ndarray,
    p0: np.ndarray,
    *,
    c: float = 0.1,
    n_passes: int = 3,
    anchor_thr: float = 0.0,
    model: str = "multinomial",
    n_bootstrap: int = 0,
    gene_mask: Optional[np.ndarray] = None,
    non_anchor_weight: float = 0.1,
    max_iter: int = 150,
) -> tuple:
    """Refine protein proportions ``p0`` (S x K) using raw RNA counts ``Y`` (S x G) and nuclei ``n`` (S).

    Estimand: an RNA-mass fraction exported as a cell fraction under an equal-RNA-yield
    assumption. Returns (p, diagnostics); rows of p sum to 1.

    Edge rules: spots with n_s = 0 are excluded from profile learning but still refined;
    spots with zero library (or all mass on frozen types) return p0; unsupported types
    (profile row all zero, or sum_s n_s p0_st < 1) are frozen at p0 and the others are refit
    on the remaining mass; fewer than 2 supported types raises ValueError.
    The EM uses a fixed iteration budget (``max_iter``, default 150; tol=1e-5); hitting the cap is NOT a failure
    (``em_converged`` / ``em_iters`` are reported diagnostics only). The single-pass result is returned
    and ``diagnostics["fallback"] == "single_pass"`` only if a pass after the first produces non-finite
    values, learns rank-deficient (rank < supported types) or ill-conditioned (cond > 1e8) profiles, or
    the alternation is not contracting (with 3+ passes: final-pass movement > first-pass movement).
    """
    if model not in ("multinomial", "poisson"):
        raise ValueError("model must be 'multinomial' or 'poisson'")
    Y = np.asarray(Y, dtype=float)
    n = np.asarray(n, dtype=float)
    p0 = np.asarray(p0, dtype=float)
    validate_raw_counts(Y)
    if Y.shape[0] != p0.shape[0] or n.shape[0] != Y.shape[0]:
        raise ValueError("Y, n and p0 must have the same number of spots")
    if not np.isfinite(p0).all() or (p0 < 0).any() or not np.isfinite(n).all() or (n < 0).any():
        raise ValueError("p0 and n must be finite and non-negative")
    if gene_mask is not None:
        gene_mask = np.asarray(gene_mask, dtype=bool)
        if gene_mask.shape != (Y.shape[1],):
            raise ValueError("gene_mask must be bool[G]")
        Y = Y[:, gene_mask]
    if Y.shape[1] < _MIN_GENES:
        raise ValueError(f"Only {Y.shape[1]} genes available; need at least {_MIN_GENES} genes")

    S, K = p0.shape
    L = Y.sum(1)
    sup = (n[:, None] * p0).sum(0) >= 1.0
    diag = {
        "zero_library_spots": L <= 0,
        "zero_nuclei_spots": n <= 0,
        "fallback": None,
        "model": model,
        "anchor_thr": float(anchor_thr),
    }
    while True:
        if sup.sum() < 2:
            raise ValueError("Fewer than 2 supported cell types; refinement is not possible")
        try:
            out = _run(Y, n, p0, L, sup, c, n_passes, anchor_thr, non_anchor_weight, model, n_bootstrap, diag, max_iter)
            break
        except _DropType as e:
            sup = sup.copy()
            sup[e.k] = False
    diag["supported_types"] = sup
    return out, diag


def _run(Y, n, p0, L, sup, c, n_passes, anchor_thr, wnon, model, n_bootstrap, diag, max_iter=150):
    S, K = p0.shape
    sidx = np.flatnonzero(sup)
    rem = p0[:, sup].sum(1)
    active = (L > 0) & (rem > 0)
    Ya, na = Y[active], n[active]
    p0a = p0[active][:, sup]
    if not sup.all():
        p0a = p0a / rem[active][:, None]
    first = None
    P = p0a
    movement, iters, conv = [], [], []
    diag["fallback"] = None
    for k in range(n_passes):
        mu = learn_profiles(Ya, na, P, anchor_thr=anchor_thr, non_anchor_weight=wnon)
        rank, cond, cos = _profile_diagnostics(mu)
        if k == 0 and (mu.sum(1) <= 0).any():
            raise _DropType(int(sidx[np.flatnonzero(mu.sum(1) <= 0)[0]]))
        diag["profile_rank"], diag["profile_condition"] = rank, cond
        full_cos = np.zeros((K, K))
        full_cos[np.ix_(sidx, sidx)] = cos
        diag["profile_cosine"] = full_cos
        diag["high_cosine_pairs"] = [
            (int(sidx[a]), int(sidx[b]))
            for a in range(len(sidx))
            for b in range(a + 1, len(sidx))
            if cos[a, b] > _HIGH_COSINE
        ]
        if k == 0 and n_bootstrap > 0:
            cv = np.full(K, np.nan)
            cv[sidx] = _bootstrap_cv(Ya, na, P, n_bootstrap, anchor_thr)
            diag["profile_bootstrap_cv"] = cv
        bad = k > 0 and (rank < len(sidx) or cond > _MAX_COND)
        phi = _floor_profiles(mu)
        if model == "multinomial":
            phi = phi / phi.sum(1, keepdims=True)
        q, info = em_refine(Ya, phi, p0a, c=c, model=model, max_iter=max_iter)
        iters.append(info["iters"])
        conv.append(bool(info["converged"]))
        if k > 0 and (bad or not np.isfinite(q).all()):
            diag["fallback"] = "single_pass"
            q = first
            break
        if not np.isfinite(q).all():
            raise RuntimeError("Non-finite proportions in the first refinement pass")
        movement.append(float(np.abs(q - P).max()) if q.size else 0.0)
        if k == 0:
            first = q
        P = q
    if diag["fallback"] is None and len(movement) >= 3 and movement[-1] > movement[0]:
        diag["fallback"] = "single_pass"
        q = first
    diag["pass_movement"], diag["em_iters"], diag["em_converged"] = movement, iters, conv
    out = p0.copy()
    sub = q if sup.all() else q * rem[active][:, None]
    cols = out[active]
    cols[:, sup] = sub
    out[active] = cols
    return out


def select_prior_weight(
    Y: np.ndarray,
    n: np.ndarray,
    p0: np.ndarray,
    *,
    grid: tuple = (0.0, 0.03, 0.1, 0.3, 1.0),
    n_folds: int = 5,
    seed: int = 0,
) -> tuple:
    """Choose the prior weight ``c`` by gene-wise cross-fitting (uses only Y, n, p0; no ground truth).

    Genes are split into ``n_folds`` folds (``default_rng(seed).permutation``). For each fold, profiles
    are learned for all genes from ``p0`` (protein only); the refinement is fitted for each ``c`` with
    ``gene_mask`` = genes outside the fold; held-out counts are predicted as
    ``L_heldout * normalize(p_s * mu[:, fold])`` and scored by multinomial deviance over spots with
    held-out library > 0. Deviance is summed over folds; the minimum wins, ties go to the smaller c.

    Returns (c, {"deviance": {c: float}, "boundary": bool, "fold_masks": [bool[G] held-out masks]}).
    ``boundary`` is True when the pick is the smallest or largest grid value.
    """
    Y = np.asarray(Y, dtype=float)
    n = np.asarray(n, dtype=float)
    p0 = np.asarray(p0, dtype=float)
    grid = tuple(grid)
    G = Y.shape[1]
    folds = np.array_split(np.random.default_rng(seed).permutation(G), n_folds)
    mu = _floor_profiles(learn_profiles(Y, n, p0))
    dev = {c: 0.0 for c in grid}
    masks = []
    for f in folds:
        held = np.zeros(G, dtype=bool)
        held[f] = True
        masks.append(held)
        Yh = Y[:, held]
        Lh = Yh.sum(1)
        ok = Lh > 0
        for c in grid:
            p, _ = refine_proportions(Y, n, p0, c=c, gene_mask=~held)
            pi = p @ mu[:, held]
            pi = np.clip(pi / pi.sum(1, keepdims=True), 1e-12, None)
            y, pp, L = Yh[ok], pi[ok], Lh[ok][:, None]
            with np.errstate(divide="ignore", invalid="ignore"):
                term = np.where(y > 0, y * np.log(y / (L * pp)), 0.0)
            dev[c] += float(2.0 * term.sum())
    best = min(grid, key=lambda c: (dev[c], c))
    boundary = best in (min(grid), max(grid))
    return best, {"deviance": dev, "boundary": bool(boundary), "fold_masks": masks}
