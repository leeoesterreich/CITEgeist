"""HISTORICAL, synthetic-calibration-only (rev 7). NOT the production interval: the production estimator is ``descriptive.gmr_interval``, and
no production module imports this one (tested). This module produced the 0.874-0.900 coverage record (synthetic/g6_fix/boot_*_s6.json) and the
SE-underestimate diagnostic; it ridge-regularises singular fits, accepts an interval with > 50% finite resamples and also intervals the odds ratio,
all of which the production contract forbids.

Block-bootstrap intervals for the magnitude GMR and the participation odds ratio (no p-values).

Resample spatial blocks with replacement (2000, seed 20261014 in production); each resample is a weighted refit of the
SAME model (magnitude: log score among actives ~ group + covariates + block fixed effects; participation: logistic,
group + covariates). Percentile 95% interval.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

import hurdle as H


def _mag_design(score, group, cov, block):
    act = score > 0
    y = np.log(score[act])
    g = group[act].astype(float)
    b = block[act]
    _, b = np.unique(b, return_inverse=True)
    n = len(y)
    cols = [np.ones(n), g]
    if cov is not None and np.asarray(cov).shape[1]:
        cols += list(np.asarray(cov, float)[act].T)
    D = np.zeros((n, b.max() + 1))
    D[np.arange(n), b] = 1.0
    return y, np.column_stack(cols + [D[:, 1:]]), np.flatnonzero(act), block[act]


def _wls_group(y, X, w):
    """Group coefficient of the weighted LS fit (normal equations; tiny ridge for empty/duplicated block columns)."""
    keep = w > 0
    Xk, wk = X[keep], w[keep]
    A = Xk.T @ (wk[:, None] * Xk) + 1e-9 * np.eye(X.shape[1])
    return np.linalg.solve(A, Xk.T @ (wk * y[keep]))[1]


def block_bootstrap(score, group, cov, block, B: int = 2000, seed: int = 20261014) -> Dict:
    score = np.asarray(score, float)
    group = np.asarray(group).astype(int)
    block = np.asarray(block)
    ub = np.unique(block)
    rng = np.random.default_rng(seed)
    y, X, act_idx, blk_act = _mag_design(score, group, cov, block)
    a = (score > 0).astype(float)
    Xp = np.column_stack([np.ones(len(a)), group] + ([np.asarray(cov, float)] if cov is not None and np.asarray(cov).shape[1] else []))
    pos_of = {u: i for i, u in enumerate(ub)}
    blk_idx_all = np.array([pos_of[v] for v in block])
    blk_idx_act = np.array([pos_of[v] for v in blk_act])
    gm, orr = np.full(B, np.nan), np.full(B, np.nan)
    for k in range(B):
        mult = np.bincount(rng.integers(0, len(ub), len(ub)), minlength=len(ub)).astype(float)
        try:
            gm[k] = np.exp(_wls_group(y, X, mult[blk_idx_act]))
        except np.linalg.LinAlgError:
            pass
        w = mult[blk_idx_all]
        beta, _, conv = _irls_w(Xp, a, w)
        if conv and abs(beta[1]) < 15:
            orr[k] = np.exp(beta[1])
    out = {}
    for name, v in (("gmr", gm), ("odds_ratio", orr)):
        ok = v[np.isfinite(v)]
        out[name] = {"lo": float(np.quantile(ok, 0.025)), "hi": float(np.quantile(ok, 0.975)), "n_ok": int(len(ok))} if len(ok) > B // 2 else {"lo": np.nan, "hi": np.nan, "n_ok": int(len(ok))}
    return out


def _irls_w(X, a, w, iters=40, ridge=1e-8):
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        mu = 1 / (1 + np.exp(-np.clip(X @ beta, -30, 30)))
        ww = w * mu * (1 - mu)
        H_ = X.T @ (ww[:, None] * X) + ridge * np.eye(X.shape[1])
        step = np.linalg.solve(H_, X.T @ (w * (a - mu)))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-8:
            return beta, mu, True
    return beta, mu, False
