"""Negative-binomial mixture for protein count channels (CPU only).

Counts are modelled as background/foreground NB components, as totalVI does for
protein counts (Gayoso et al. 2021, Nat Methods). Each spot has a size factor
sf_s, and count_s ~ sum_k w_k NB(mean = sf_s * mu_k, dispersion theta_k) with the
NB2 variance m + m^2 / theta. One or two components are compared by BIC.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize
from scipy.special import digamma, gammaln

_THETA_MIN, _THETA_MAX, _MU_MIN = 1e-3, 1e6, 1e-8


@dataclass
class NBMixtureFit:
    """Fit of a 1- or 2-component NB mixture.

    means are per unit size factor, ascending. positive_mask marks spots whose
    posterior for the highest-mean component is > 0.5 (all False if 1 component).
    """

    n_components: int
    means: np.ndarray
    dispersions: np.ndarray
    weights: np.ndarray
    positive_mask: np.ndarray
    bic: dict = field(default_factory=dict)


def size_factors_from_counts(X: np.ndarray) -> np.ndarray:
    """Per-spot total counts / median total; zero-total spots get 1.0."""
    X = np.asarray(X, dtype=float)
    tot = X.sum(axis=1)
    med = np.median(tot)
    if med <= 0:
        return np.ones_like(tot)
    sf = tot / med
    sf[tot <= 0] = 1.0
    return sf


def _nb_logpmf(c, m, theta):
    """NB2 log pmf: mean m, dispersion theta (var = m + m^2/theta)."""
    return (
        gammaln(c + theta)
        - gammaln(theta)
        - gammaln(c + 1.0)
        + theta * (np.log(theta) - np.log(theta + m))
        + c * (np.log(m) - np.log(theta + m))
    )


def _fit_component(c, sf, w, x0):
    """Maximize the weighted NB log-likelihood over (log mu, log theta)."""

    def nll(p):
        mu = max(np.exp(p[0]), _MU_MIN)
        th = float(np.clip(np.exp(p[1]), _THETA_MIN, _THETA_MAX))
        m = sf * mu
        ll = _nb_logpmf(c, m, th)
        d_lm = c - (c + th) * m / (th + m)
        d_lt = th * (digamma(c + th) - digamma(th) + np.log(th) - np.log(th + m) + (m - c) / (th + m))
        return -np.sum(w * ll), -np.array([np.sum(w * d_lm), np.sum(w * d_lt)])

    lo, hi = np.log(_MU_MIN), np.log(1e9)
    res = minimize(
        nll,
        x0,
        jac=True,
        method="L-BFGS-B",
        bounds=[(lo, hi), (np.log(_THETA_MIN), np.log(_THETA_MAX))],
        options={"maxiter": 30},
    )
    return float(np.exp(res.x[0])), float(np.clip(np.exp(res.x[1]), _THETA_MIN, _THETA_MAX))


def _init_params(c, sf, r):
    mu = max(float(np.sum(r * c) / max(np.sum(r * sf), 1e-12)), 1e-3)
    return np.array([np.log(mu), np.log(1.0)])


def _em(c, sf, resp, max_iter, tol):
    n, k = resp.shape
    params = [_init_params(c, sf, resp[:, j]) for j in range(k)]
    mus = np.zeros(k)
    ths = np.ones(k)
    prev = -np.inf
    ll = -np.inf
    w = resp.mean(axis=0)
    for _ in range(max_iter):
        for j in range(k):
            mus[j], ths[j] = _fit_component(c, sf, resp[:, j], params[j])
            params[j] = np.array([np.log(mus[j]), np.log(ths[j])])
        w = np.clip(resp.mean(axis=0), 1e-6, None)
        w /= w.sum()
        logp = np.column_stack([np.log(w[j]) + _nb_logpmf(c, sf * mus[j], ths[j]) for j in range(k)])
        mx = logp.max(axis=1, keepdims=True)
        lse = mx[:, 0] + np.log(np.exp(logp - mx).sum(axis=1))
        ll = float(lse.sum())
        resp = np.exp(logp - lse[:, None])
        if abs(ll - prev) < tol * max(1.0, abs(ll)):
            break
        prev = ll
    return ll, mus, ths, w, resp


def fit_nb_mixture(
    counts: np.ndarray,
    size_factors: np.ndarray | None = None,
    max_components: int = 2,
    seed: int = 1234,
    max_iter: int = 200,
    tol: float = 1e-6,
) -> NBMixtureFit:
    """Fit 1 and 2 component NB mixtures by EM and keep the lowest BIC.

    max_components > 2 is not supported (only 1 or 2 are compared).
    Constant or all-zero input returns one component without fitting.
    """
    c = np.asarray(counts, dtype=float).ravel()
    n = c.shape[0]
    sf = np.ones(n) if size_factors is None else np.asarray(size_factors, dtype=float).ravel()
    sf = np.where(sf > 0, sf, 1.0)
    none = np.zeros(n, dtype=bool)
    if n == 0 or np.std(c) < 1e-10:
        mu = float(np.mean(c / sf)) if n else 0.0
        return NBMixtureFit(1, np.array([mu]), np.array([np.nan]), np.array([1.0]), none, {1: float("nan")})

    # k = 1
    mu0 = max(float(c.sum() / sf.sum()), 1e-3)
    v = float(np.var(c / sf))
    th0 = np.clip(mu0**2 / max(v - mu0, 1e-3), 0.1, 1e3)
    mu1, th1 = _fit_component(c, sf, np.ones(n), np.array([np.log(mu0), np.log(th0)]))
    ll1 = float(np.sum(_nb_logpmf(c, sf * mu1, th1)))
    bic = {1: 2 * np.log(n) - 2 * ll1}
    best = (1, np.array([mu1]), np.array([th1]), np.array([1.0]), none)
    if max_components < 2:
        return NBMixtureFit(*best, bic=bic)

    # k = 2
    rng = np.random.default_rng(seed)
    norm = c / sf
    zero_frac = float(np.mean(c == 0))
    thr0 = 0.0 if zero_frac > 0.5 else float(np.median(norm))
    best_ll, best2 = -np.inf, None
    for i in range(3):
        thr = thr0 if i == 0 else thr0 * rng.uniform(0.6, 1.6) + (rng.uniform(0.2, 1.0) if thr0 == 0 else 0.0)
        hi = norm > thr
        if hi.all() or (~hi).all():
            continue
        resp = np.column_stack([np.where(hi, 0.1, 0.9), np.where(hi, 0.9, 0.1)])
        try:
            ll, mus, ths, w, resp = _em(c, sf, resp, max_iter, tol)
        except (FloatingPointError, ValueError):
            continue
        if np.isfinite(ll) and ll > best_ll:
            best_ll, best2 = ll, (mus.copy(), ths.copy(), w.copy(), resp.copy())
    if best2 is not None:
        mus, ths, w, resp = best2
        bic[2] = 5 * np.log(n) - 2 * best_ll
        if bic[2] < bic[1]:
            order = np.argsort(mus)
            mus, ths, w, resp = mus[order], ths[order], w[order], resp[:, order]
            best = (2, mus, ths, w, resp[:, 1] > 0.5)
    return NBMixtureFit(*best, bic=bic)
