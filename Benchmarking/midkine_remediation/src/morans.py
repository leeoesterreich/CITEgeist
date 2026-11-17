"""Bivariate Moran's I and Moran spectral randomization (MSR) (protocol ``morans``; spec 9.5, G7).

Statistic: I(x, y) = (1/n) * sum_i zx_i * (W zy)_i with z = (v - mean)/sd (ddof=0) and W the
row-standardised KNN(k=6) weights. This is the convention of the legacy manual fallback in
the original ``spatial_morans.py``; the ``esda`` value differs by a known factor (tested).

Null (primary): MSR with sign flips (Wagner & Dray 2015, singleton variant).
  * spectrum from W_sym = (W + W^T)/2 doubly centred, ``canonical_eigh`` -> V, lambda (``eigh`` followed by a deterministic choice of
    the basis inside tied eigenspaces and of every vector's sign, so the surrogates do not depend on the BLAS thread count);
  * the continuous field y is replaced by  mean(y) + V diag(s) V^T (y - mean(y)),  s_k = +-1 iid;
    the zero-inflated x is held fixed;
  * because I is linear in y*, I(x, y*) = sum_k s_k c_k b_k / (n sd_y) with c = V^T yc and
    b = V^T (W^T zx): each surrogate costs O(n) after two matrix-vector products.

Definition A (cells inside spots): the spot-mean field is surrogated on the spot graph and
each cell keeps its within-spot deviation (``SpotNull``).
"""
from __future__ import annotations

import resource
import time
from typing import Dict, Optional

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

DENSE_EIGH_MAX_N = 6000          # validated ceiling for the dense eigendecomposition (see run_g7 timing)


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    sd = v.std()
    if sd == 0:
        raise ValueError("constant vector")
    return (v - v.mean()) / sd


def morans_bv(x, y, W: sp.csr_matrix, rank: bool = False) -> float:
    """I_xy = (1/n) zx . (W zy)."""
    if rank:
        x, y = rankdata(x), rankdata(y)
    zx, zy = zscore(x), zscore(y)
    return float(zx @ (W @ zy) / len(zx))


def morans_bv_manual(x, y, nbrs: np.ndarray) -> float:
    """Independent loop implementation from a neighbour index matrix (n, k) (row-standardised = mean over neighbours)."""
    n = len(x)
    xm, ym = sum(x) / n, sum(y) / n
    sx = (sum((v - xm) ** 2 for v in x) / n) ** 0.5
    sy = (sum((v - ym) ** 2 for v in y) / n) ** 0.5
    tot = 0.0
    for i in range(n):
        lag = sum((y[j] - ym) / sy for j in nbrs[i]) / len(nbrs[i])
        tot += (x[i] - xm) / sx * lag
    return tot / n


EIG_RESID_TOL = 1e-7             # self-check: max|A V - V lam|; genuine noise ~1e-12, a chained cluster may legitimately carry its width (a few x 1e-8), a lost-rank QR gives O(1)
EIG_ORTH_TOL = 1e-10             # self-check: max|V^T V - I|; genuine ~1e-13
EIG_TIE_TOL = 1e-8               # eigenvalues closer than this (chained, ascending order) form one tied eigenspace; the spectrum lies in [-1.1, 1.1]


def _probe(n: int, m: int, first: int = 0) -> np.ndarray:
    """Columns ``first``..``first+m-1`` of a fixed hash matrix, entries in [-0.5, 0.5): a splitmix-style integer hash, exact and identical on every machine (no RNG, no libm)."""
    i = np.arange(n, dtype=np.uint64)[:, None]
    j = np.arange(m, dtype=np.uint64)[None, :]
    h = i * np.uint64(0x9E3779B97F4A7C15) + (j + np.uint64(first + 1)) * np.uint64(0xC2B2AE3D27D4EB4F)
    h ^= h >> np.uint64(30)
    h *= np.uint64(0xBF58476D1CE4E5B9)
    h ^= h >> np.uint64(27)
    h *= np.uint64(0x94D049BB133111EB)
    h ^= h >> np.uint64(31)
    return (h >> np.uint64(11)).astype(np.float64) / 2.0**53 - 0.5


def canonical_eigh(A: np.ndarray, tol: float = EIG_TIE_TOL, _eigh=np.linalg.eigh):
    """``eigh`` with a basis that is a function of the matrix only (LAPACK returns an arbitrary basis in a tied eigenspace and an arbitrary sign,
    and which one depends on the BLAS thread count).

    1. Eigenvalues are clustered in ascending order (a gap > ``tol`` starts a new cluster). Eigenvalues and non-tied vectors are kept.
    2. A tied cluster with columns Vc spans an eigenspace whose projector P = Vc Vc^T does not depend on the basis. It is replaced by the
       Gram-Schmidt (Householder QR) orthonormalisation of P applied to the columns of a fixed hash matrix (``_probe`` columns 1..m, in column order; column 0 is the sign probe and is not reused, else it
       would be exactly orthogonal to every Gram-Schmidt column after the first).
    3. Every vector gets the sign for which its dot product with a fixed probe vector is positive.
    4. Self-check (raises ValueError): max|A V - V lam| <= EIG_RESID_TOL and max|V^T V - I| <= EIG_ORTH_TOL, so a rank-deficient QR step (vectors
       outside the eigenspace) or a corrupted ``eigh`` answer can never pass silently.
    ``_eigh`` is a test hook (inject a different but equally valid eigh answer).
    """
    lam, V = _eigh(A)
    lam, V = np.array(lam), np.array(V)
    n = len(lam)
    cuts = np.r_[0, np.flatnonzero(np.diff(lam) > tol) + 1, n]
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b - a > 1:
            Vc = V[:, a:b]
            V[:, a:b] = np.linalg.qr(Vc @ (Vc.T @ _probe(n, b - a, first=1)))[0]
    d = _probe(n, 1)[:, 0] @ V
    if (np.abs(d) < 1e-10).any():
        raise ValueError("cannot fix the sign of an eigenvector (probe is orthogonal to it)")
    V = V * np.sign(d)
    resid = float(np.abs(A @ V - V * lam).max())
    orth = float(np.abs(V.T @ V - np.eye(n)).max())
    if resid > EIG_RESID_TOL or orth > EIG_ORTH_TOL:
        raise ValueError(f"canonical_eigh self-check failed: max|A V - V lam| = {resid:.3g} (tol {EIG_RESID_TOL}), max|V^T V - I| = {orth:.3g} (tol {EIG_ORTH_TOL})")
    return lam, V


class SpectralNull:
    """Eigen-decomposition of the doubly-centred symmetrised weights (``canonical_eigh``: thread-count independent)."""

    def __init__(self, W: sp.csr_matrix):
        n = W.shape[0]
        if n > DENSE_EIGH_MAX_N:
            raise ValueError(f"n={n} exceeds the validated dense-eigh ceiling {DENSE_EIGH_MAX_N}; use the spot-level construction")
        A = ((W + W.T) / 2).toarray()
        A -= A.mean(axis=0, keepdims=True)
        A -= A.mean(axis=1, keepdims=True)
        self.lam, self.V = canonical_eigh(A)
        self.n = n
        self.W = W

    def pvalue(self, x, y, B: int = 999, seed: int = 42, rank: bool = False) -> Dict[str, float]:
        """MSR p-values for I(x, y); x fixed, y surrogated. Two-sided and one-sided (positive)."""
        if rank:
            x, y = rankdata(x), rankdata(y)
        n = self.n
        zx = zscore(x)
        yc = np.asarray(y, float) - np.mean(y)
        sd_y = yc.std()
        a = self.W.T @ zx
        c = self.V.T @ yc
        b = self.V.T @ a
        w = c * b
        I_obs = float(w.sum() / (n * sd_y))
        rng = np.random.default_rng(seed)
        S = rng.integers(0, 2, size=(B, n), dtype=np.int8) * 2 - 1
        Istar = (S @ w) / (n * sd_y)
        p2 = (1 + np.sum(np.abs(Istar) >= abs(I_obs) - 1e-12)) / (B + 1)
        p1 = (1 + np.sum(Istar >= I_obs - 1e-12)) / (B + 1)
        return {"I": I_obs, "p_two_sided": float(p2), "p_upper": float(p1), "n_surrogates": B}

    def surrogates(self, y, B: int, rng: np.random.Generator) -> np.ndarray:
        """Explicit surrogate fields, shape (n, B) (used by the definition-A construction and tests)."""
        yc = np.asarray(y, float) - np.mean(y)
        c = self.V.T @ yc
        S = rng.integers(0, 2, size=(self.n, B), dtype=np.int8) * 2 - 1
        return np.mean(y) + self.V @ (S * c[:, None])


def naive_permutation_pvalue(x, y, W: sp.csr_matrix, B: int = 999, seed: int = 42, rank: bool = False) -> Dict[str, float]:
    """Free permutation of y (the null the protocol rejects); diagnostic for the G7 mutation check."""
    if rank:
        x, y = rankdata(x), rankdata(y)
    n = len(x)
    zx, zy = zscore(x), zscore(y)
    a = W.T @ zx
    I_obs = float(a @ zy / n)
    rng = np.random.default_rng(seed)
    P = rng.permuted(np.tile(np.arange(n), (B, 1)), axis=1)
    Istar = zy[P] @ a / n
    return {"I": I_obs, "p_two_sided": float((1 + np.sum(np.abs(Istar) >= abs(I_obs) - 1e-12)) / (B + 1)),
            "p_upper": float((1 + np.sum(Istar >= I_obs - 1e-12)) / (B + 1))}


class SpotNull:
    """Definition-A null: spot-level spectral surrogates + each cell's within-spot deviation.

    ``spot_of_cell`` maps every cell to a spot row of the spot graph used by ``spot_null``.
    The cell-level statistic uses ``W_cell`` (row-standardised KNN among cells) and cell-level z-scores.
    """

    def __init__(self, spot_null: SpectralNull, spot_of_cell: np.ndarray, W_cell: sp.csr_matrix):
        self.sn = spot_null
        self.spot = np.asarray(spot_of_cell)
        self.Wc = W_cell
        self.S = spot_null.n
        self.m = np.bincount(self.spot, minlength=self.S).astype(float)
        if (self.m == 0).any():
            raise ValueError("every spot of the spot graph must contain at least one cell")

    def pvalue(self, x_cell, y_cell, B: int = 999, seed: int = 42, chunk: int = 200) -> Dict[str, float]:
        n = len(x_cell)
        zx = zscore(x_cell)
        y = np.asarray(y_cell, float)
        mu = np.bincount(self.spot, weights=y, minlength=self.S) / self.m
        dev = y - mu[self.spot]                                  # sums to zero inside every spot
        a = self.Wc.T @ zx
        A = np.bincount(self.spot, weights=a, minlength=self.S)  # per-spot sum of a
        a_dev = float(a @ dev)
        sd_obs = y.std()
        I_obs = float((a @ (y - y.mean())) / (n * sd_obs))
        mu_c = mu - mu.mean()
        c = self.sn.V.T @ mu_c
        bA = self.sn.V.T @ A
        dev2 = float(dev @ dev)
        rng = np.random.default_rng(seed)
        Istar = np.empty(B)
        for s0 in range(0, B, chunk):
            m_ = min(chunk, B - s0)
            Sg = rng.integers(0, 2, size=(self.S, m_), dtype=np.int8) * 2 - 1
            Z = self.sn.V @ (Sg * c[:, None])                     # centred spot-mean surrogates (S, m_)
            mu_star = mu.mean() + Z
            wm = self.m @ mu_star / n                             # cell-level mean of y*
            var = (dev2 + self.m @ (mu_star**2)) / n - wm**2
            num = a_dev + (Sg * c[:, None] * bA[:, None]).sum(axis=0)
            Istar[s0:s0 + m_] = num / (n * np.sqrt(var))
        p2 = (1 + np.sum(np.abs(Istar) >= abs(I_obs) - 1e-12)) / (B + 1)
        p1 = (1 + np.sum(Istar >= I_obs - 1e-12)) / (B + 1)
        return {"I": I_obs, "p_two_sided": float(p2), "p_upper": float(p1), "n_surrogates": B}


def naive_cell_permutation_pvalue(x_cell, y_cell, W_cell: sp.csr_matrix, B: int = 999, seed: int = 42) -> Dict[str, float]:
    return naive_permutation_pvalue(x_cell, y_cell, W_cell, B, seed)


def benchmark_eigh(W: sp.csr_matrix) -> Dict[str, float]:
    """Wall time and process peak RSS (MB) of building the SpectralNull for ``W``."""
    t = time.time()
    sn = SpectralNull(W)
    return {"n": int(W.shape[0]), "eigh_seconds": time.time() - t, "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0,
            "spectrum_min": float(sn.lam.min()), "spectrum_max": float(sn.lam.max())}
