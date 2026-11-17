"""HISTORICAL (G6 record only): permutation tests, BH and surrogate-label nulls of protocol revisions 1-4.

NOT part of the primary (descriptive) API and NOT used by any decision: G6 failed on 2026-10-15 and the hurdle
reports no p-value or q-value (protocol rev 5/6). The code is kept so the committed synthetic G6 record
(synthetic/g6, synthetic/g6_fix) stays reproducible and mutation-checked. Names carry ``historical_``.

Magnitude: Freedman-Lane within-block permutation of reduced-model residuals. Participation: within-block label
permutation with the score statistic. ``scheme="plain"`` is the G6 diagnostic (free permutation; must be inflated).
"""
from __future__ import annotations

from math import lgamma, log
from typing import Dict, Optional

import numpy as np

from hurdle import MIN_ACTIVE_PER_GROUP, MIN_PER_GROUP, _irls, _orth_basis, standardise

MIN_MIXED_BLOCKS = 10
MIN_DISTINCT = 999
DEFAULT_B = 4999
DEFAULT_SEED = 20261014


def safe_log(x: np.ndarray) -> np.ndarray:
    """log with non-positive entries replaced by half the smallest positive value (the G6-era floor; NOT used by the primary API)."""
    x = np.asarray(x, float)
    pos = x[x > 0]
    floor = pos.min() / 2 if pos.size else 1.0
    return np.log(np.where(x > 0, x, floor))


def historical_covariate_matrix(obs, drop: tuple = ()) -> np.ndarray:
    """G6-era covariate builder (4 columns incl. log(n_st_eff), ``drop`` names omissions, safe_log floor). Reproduces the committed synthetic G6 record."""
    cols = {
        "log_spot_umi": safe_log(obs["spot_umi"].to_numpy()),
        "log_rna_per_cell_equiv": safe_log(obs["rna_per_cell_equiv"].to_numpy()),
        "qp_prop": obs["qp_prop"].to_numpy(float),
        "log_n_st_eff": safe_log(obs["n_st_eff"].to_numpy()),
    }
    keep = [k for k in cols if k not in drop]
    return standardise(np.column_stack([cols[k] for k in keep])) if keep else np.zeros((len(obs), 0))


# --------------------------------------------------------------------------- permutations
def _groups(block: np.ndarray):
    order = np.argsort(block, kind="stable")
    b = block[order]
    cuts = np.flatnonzero(np.diff(b)) + 1
    return np.split(order, cuts)


def within_block_perms(block: np.ndarray, B: int, rng: np.random.Generator) -> np.ndarray:
    """(B, n) index matrix; row b maps position i -> source index, permuting only inside blocks."""
    n = len(block)
    P = np.empty((B, n), dtype=np.int32)
    for idx in _groups(np.asarray(block)):
        if len(idx) == 1:
            P[:, idx] = idx[0]
        else:
            P[:, idx] = rng.permuted(np.tile(idx.astype(np.int32), (B, 1)), axis=1)
    return P


def plain_perms(n: int, B: int, rng: np.random.Generator) -> np.ndarray:
    return rng.permuted(np.tile(np.arange(n, dtype=np.int32), (B, 1)), axis=1)


def _unique_rows(P: np.ndarray, identity: np.ndarray, keep: int) -> np.ndarray:
    seen = {identity.tobytes()}
    out = []
    for row in P:
        k = row.tobytes()
        if k not in seen:
            seen.add(k)
            out.append(row)
            if len(out) == keep:
                break
    return np.array(out)


def _perm_matrix(block, B, rng, scheme, log_distinct):
    n = len(block)
    ident = np.arange(n, dtype=np.int32)
    gen = (lambda m: within_block_perms(block, m, rng)) if scheme == "block" else (lambda m: plain_perms(n, m, rng))
    if scheme == "block" and log_distinct < log(5.0 * B):
        # few distinct permutations: sample WITHOUT replacement (unique rows), at most `distinct - 1` of them
        target = int(min(B, np.floor(np.exp(log_distinct)) - 1))
        for _ in range(50):
            P = _unique_rows(gen(max(20 * target, 2000)), ident, target)
            if len(P) >= target:
                return P
        return P
    return gen(B)


# --------------------------------------------------------------------------- estimability
def _mixed_and_distinct(block, group, kind):
    mixed = 0
    logd = 0.0
    for idx in _groups(np.asarray(block)):
        k = int(group[idx].sum())
        m = len(idx)
        if 0 < k < m:
            mixed += 1
        if kind == "participation":
            logd += lgamma(m + 1) - lgamma(k + 1) - lgamma(m - k + 1)
        else:
            logd += lgamma(m + 1)
    return mixed, logd


def historical_estimability(group, block, kind, active=None) -> Dict:
    g = np.asarray(group).astype(int)
    if kind == "magnitude":
        sel = np.asarray(active, bool)
        g, block = g[sel], np.asarray(block)[sel]
        n1, n0 = int(g.sum()), int((1 - g).sum())
        if min(n1, n0) < MIN_ACTIVE_PER_GROUP:
            return {"ok": False, "reason": f"active senders per group {n1}/{n0} < {MIN_ACTIVE_PER_GROUP}"}
    else:
        n1, n0 = int(g.sum()), int((1 - g).sum())
        if min(n1, n0) < MIN_PER_GROUP:
            return {"ok": False, "reason": f"compartments per group {n1}/{n0} < {MIN_PER_GROUP}"}
    mixed, logd = _mixed_and_distinct(block, g, kind)
    if mixed < MIN_MIXED_BLOCKS:
        return {"ok": False, "reason": f"{mixed} mixed-label blocks < {MIN_MIXED_BLOCKS}", "n_mixed_blocks": mixed}
    if logd < log(MIN_DISTINCT):
        return {"ok": False, "reason": f"distinct permutations exp({logd:.2f}) < {MIN_DISTINCT}", "n_mixed_blocks": mixed}
    return {"ok": True, "n_mixed_blocks": mixed, "log_distinct": logd}



# --------------------------------------------------------------------------- magnitude
def historical_magnitude_permutation_test(score, group, cov, block, B: int = DEFAULT_B, seed: int = DEFAULT_SEED, scheme: str = "block", active: Optional[np.ndarray] = None, label_matrix: Optional[np.ndarray] = None) -> Dict:
    """Freedman-Lane block-permutation test of the group term among active senders.

    ``active`` defaults to ``score > 0`` (passing ``~active`` is the swap_nonzero_filter mutation).
    ``label_matrix`` (B, n_all) replaces the permutation null by externally generated surrogate labels
    (CANDIDATE scheme, not in the protocol; see ``surrogate_labels``).
    """
    score = np.asarray(score, float)
    group = np.asarray(group).astype(int)
    block = np.asarray(block)
    act = (score > 0) if active is None else np.asarray(active, bool)
    est = historical_estimability(group, block, "magnitude", act)
    if not est["ok"]:
        return {"status": "not_estimable", "reason": est["reason"]}
    if (score[act] <= 0).any():
        raise ValueError("magnitude set contains non-positive scores (filter must select score > 0)")
    y = np.log(score[act])
    g = group[act].astype(float)
    blk = block[act]
    _, blk = np.unique(blk, return_inverse=True)
    n = len(y)
    cols = [np.ones(n)]
    if cov is not None and np.asarray(cov).shape[1]:
        cols += list(np.asarray(cov, float)[act].T)
    nb = blk.max() + 1
    D = np.zeros((n, nb))
    D[np.arange(n), blk] = 1.0
    X0 = np.column_stack(cols + [D[:, 1:]]) if nb > 1 else np.column_stack(cols)
    U, s_ = _orth_basis(X0)
    rank = U.shape[1]
    proj = lambda v: U @ (U.T @ v)  # noqa: E731
    gt = g - proj(g)
    ss = float(gt @ gt)
    if ss <= 1e-12:
        return {"status": "not_estimable", "reason": "group term collinear with covariates/blocks"}
    r0 = y - proj(y)
    rss0 = float(r0 @ r0)
    df = n - rank - 1
    b_obs = float(gt @ y) / ss
    rss_obs = rss0 - b_obs**2 * ss
    se = np.sqrt(max(rss_obs, 1e-300) / df / ss)
    t_obs = b_obs / se
    if label_matrix is not None:
        Gs = np.asarray(label_matrix)[:, act].astype(float).T          # (n, B) surrogate labels
        Gt = Gs - U @ (U.T @ Gs)
        ssb = (Gt**2).sum(axis=0)
        okb = ssb > 1e-12
        safe = np.where(okb, ssb, 1.0)
        bstar = np.where(okb, (Gt.T @ r0) / safe, 0.0)
        rss_star = rss0 - bstar**2 * ssb
        tstar = np.where(okb, bstar / np.sqrt(np.maximum(rss_star, 1e-300) / df / safe), 0.0)
    else:
        rng = np.random.default_rng(seed)
        P = _perm_matrix(blk, B, rng, scheme, est["log_distinct"])
        rstar = r0[P]
        bstar = (rstar @ gt) / ss
        rss_star = rss0 - bstar**2 * ss
        tstar = bstar / np.sqrt(np.maximum(rss_star, 1e-300) / df / ss)
    p = (1 + np.sum(np.abs(tstar) >= abs(t_obs) - 1e-12)) / (len(tstar) + 1)
    return {
        "status": "ok", "reason": "", "coef": b_obs, "gmr": float(np.exp(b_obs)), "se": float(se), "t": float(t_obs),
        "p": float(p), "t_null_sd": float(np.std(tstar)), "sign_positive": bool(b_obs > 0), "n_active": int(n), "n_active_high": int(g.sum()),
        "n_mixed_blocks": est["n_mixed_blocks"], "n_perm": int(len(tstar)),
    }


# --------------------------------------------------------------------------- participation
def historical_participation_permutation_test(active, group, cov, block, B: int = DEFAULT_B, seed: int = DEFAULT_SEED, scheme: str = "block", label_matrix: Optional[np.ndarray] = None) -> Dict:
    """Within-block label permutation, covariates fixed, score statistic of the group term."""
    a = np.asarray(active, float)
    g = np.asarray(group).astype(int)
    block = np.asarray(block)
    est = historical_estimability(g, block, "participation")
    if not est["ok"]:
        return {"status": "not_estimable", "reason": est["reason"]}
    n = len(a)
    X0 = np.column_stack([np.ones(n)] + ([np.asarray(cov, float)] if cov is not None and np.asarray(cov).shape[1] else []))
    beta0, mu0, conv0 = _irls(X0, a)
    if not conv0:
        return {"status": "not_estimable", "reason": "null logistic model did not converge"}
    w = mu0 * (1 - mu0)
    r = a - mu0
    Qinv = np.linalg.inv(X0.T @ (w[:, None] * X0))
    WX = w[:, None] * X0

    def score(G):  # G: (m, n) 0/1 labels
        U = G @ r
        cross = G @ WX
        Iinfo = G @ w - np.einsum("ij,jk,ik->i", cross, Qinv, cross)
        return U**2 / np.maximum(Iinfo, 1e-300)

    s_obs = float(score(g[None, :].astype(float))[0])
    if label_matrix is not None:
        s_star = score(np.asarray(label_matrix).astype(float))
    else:
        rng = np.random.default_rng(seed)
        P = _perm_matrix(block, B, rng, scheme, est["log_distinct"])
        s_star = score(g[P].astype(float))
    p = (1 + np.sum(s_star >= s_obs - 1e-12)) / (len(s_star) + 1)
    # effect: adjusted odds ratio from the full model (NaN + note under separation)
    beta, _, conv = _irls(np.column_stack([X0, g]), a)
    sep = (not conv) or abs(beta[-1]) > 15
    return {
        "status": "ok", "reason": "", "score": s_obs, "p": float(p), "odds_ratio": float("nan") if sep else float(np.exp(beta[-1])),
        "sign_positive": bool(beta[-1] > 0) if not sep else bool((g @ r) > 0), "separation": bool(sep),
        "n_mixed_blocks": est["n_mixed_blocks"], "n_perm": int(len(s_star)), "n_active": int(a.sum()),
    }


def historical_h1_permutation_test(score, group, cov, block, B=DEFAULT_B, seed=DEFAULT_SEED, scheme="block", label_matrix=None) -> Dict:
    score = np.asarray(score, float)
    return {
        "magnitude": historical_magnitude_permutation_test(score, group, cov, block, B, seed, scheme, label_matrix=label_matrix),
        "participation": historical_participation_permutation_test(score > 0, group, cov, block, B, seed, scheme, label_matrix=label_matrix),
    }


def historical_surrogate_labels(kfield, spectral_null, frac_high: float, B: int, rng: np.random.Generator) -> np.ndarray:
    """CANDIDATE null (not in the protocol): spectral sign-flip surrogates of the continuous keratin field, re-thresholded.

    ``kfield`` is the continuous keratin score the label is cut from (rank-normal-scored so skew does not
    drive the surrogate); ``spectral_null`` is a ``morans.SpectralNull`` on the eligible compartments' KNN
    graph. Each surrogate keeps the field's spectrum (hence label clustering) and the observed fraction of
    high compartments. Returns (B, n) int8 labels.
    """
    from scipy.stats import norm, rankdata

    n = len(kfield)
    z = norm.ppf((rankdata(kfield) - 0.5) / n)
    zc = z - z.mean()
    c = spectral_null.V.T @ zc
    S = rng.integers(0, 2, size=(n, B), dtype=np.int8) * 2 - 1
    Z = spectral_null.V @ (S * c[:, None])
    k = max(1, int(round(frac_high * n)))
    thr = np.partition(Z, n - k, axis=0)[n - k]
    return (Z >= thr[None, :]).T.astype(np.int8)


def historical_benjamini_hochberg(p) -> np.ndarray:
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = np.empty(m)
    ranked = p[o] * m / (np.arange(m) + 1)
    q[o] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    return q


