"""CD3E-anchored CD4/CD8 split (Amendment B): pure numpy, no AnnData.

A spot gives one number per antibody: cells x signal per cell. CD3E sits on every T cell, so regressing CD3E on the
two subtype markers gives slopes proportional to 1 / (signal per cell); their ratio is the relative brightness. After
removing background (median over T-low spots), ``z_m = beta_m * y_m`` is in CD3E count units and the split is
``pi = (z_a + 0.5) / (z_a + z_b + 1)``. A sample that cannot be fitted is NOT estimable and raises; there is no
fallback to raw counts or to 50/50.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

SMOOTHING = 0.5  # fixed, in CD3E count units (not a Jeffreys prior)


class AnchorNotEstimable(ValueError):
    """The CD3E anchor cannot be estimated for this sample (fail closed; no fallback)."""


@dataclass
class AnchorFit:
    beta_a: float
    beta_b: float
    intercept: float
    se_a: float
    se_b: float
    a_cd3e: float
    a_a: float
    a_b: float
    n_low: int
    rank: int
    cond: float
    max_leverage: float
    zero_frac_a: float
    zero_frac_b: float
    estimable: bool
    reason: str | None

    @property
    def ratio(self) -> float:
        """beta_a / beta_b, an estimate of b_b / b_a (relative brightness of marker b over marker a)."""
        return float(self.beta_a / self.beta_b) if self.beta_b else float("nan")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["ratio"] = self.ratio
        return d


def _as_vec(x, name):
    v = np.asarray(x, dtype=float)
    if v.ndim != 1:
        raise ValueError(f"{name} must be 1-D")
    return v


def fit_cd3e_anchor(x_cd3e, x_a, x_b, *, low_pct: float = 10, min_low_spots: int = 20) -> AnchorFit:
    """OLS of CD3E on markers a and b (with intercept) plus T-low background medians; flags non-estimable samples."""
    y, xa, xb = _as_vec(x_cd3e, "x_cd3e"), _as_vec(x_a, "x_a"), _as_vec(x_b, "x_b")
    if not (len(y) == len(xa) == len(xb)):
        raise ValueError("x_cd3e, x_a and x_b must have the same length")
    nan = float("nan")
    n = len(y)

    def _fail(reason, **kw):
        base = dict(
            beta_a=nan, beta_b=nan, intercept=nan, se_a=nan, se_b=nan, a_cd3e=nan, a_a=nan, a_b=nan, n_low=0,
            rank=0, cond=nan, max_leverage=nan,
            zero_frac_a=float(np.mean(xa == 0)) if n else nan,
            zero_frac_b=float(np.mean(xb == 0)) if n else nan,
            estimable=False, reason=reason,
        )  # fmt: skip
        base.update(kw)
        return AnchorFit(**base)

    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(xa)) and np.all(np.isfinite(xb))):
        return _fail("non-finite input")
    if n < 4:
        return _fail(f"fewer than 4 spots (n={n})")
    design = np.column_stack([np.ones(n), xa, xb])
    rank = int(np.linalg.matrix_rank(design))
    cond = float(np.linalg.cond(design))
    if rank < 3:
        return _fail("rank-deficient design", rank=rank, cond=cond)
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    resid = y - design @ coef
    sigma2 = float(resid @ resid) / (n - 3)
    xtx_inv = np.linalg.inv(design.T @ design)
    se = np.sqrt(sigma2 * np.diag(xtx_inv))
    max_lev = float(np.max(np.einsum("ij,jk,ik->i", design, xtx_inv, design)))
    thr = np.percentile(y, low_pct)
    low = y <= thr  # ties at the percentile included
    n_low = int(low.sum())
    a_cd3e, a_a, a_b = (float(np.median(v[low])) for v in (y, xa, xb))
    fields = dict(
        beta_a=float(coef[1]), beta_b=float(coef[2]), intercept=float(coef[0]), se_a=float(se[1]), se_b=float(se[2]),
        a_cd3e=a_cd3e, a_a=a_a, a_b=a_b, n_low=n_low, rank=rank, cond=cond, max_leverage=max_lev,
        zero_frac_a=float(np.mean(xa == 0)), zero_frac_b=float(np.mean(xb == 0)),
    )  # fmt: skip
    if not all(np.isfinite(v) for v in (*coef, *se, a_cd3e, a_a, a_b)):
        return AnchorFit(**fields, estimable=False, reason="non-finite coefficient")
    if coef[1] <= 0 or coef[2] <= 0:
        return AnchorFit(**fields, estimable=False, reason="non-positive slope")
    if n_low < min_low_spots:
        return AnchorFit(**fields, estimable=False, reason=f"fewer than {min_low_spots} T-low spots (n_low={n_low})")
    return AnchorFit(**fields, estimable=True, reason=None)


def anchored_z(x_a, x_b, fit: AnchorFit) -> tuple[np.ndarray, np.ndarray]:
    """Background-subtracted, slope-scaled evidence ``z = beta * max(x - a, 0)`` in CD3E count units."""
    if not fit.estimable:
        raise AnchorNotEstimable(f"CD3E anchor not estimable: {fit.reason}")
    xa, xb = _as_vec(x_a, "x_a"), _as_vec(x_b, "x_b")
    if not (np.all(np.isfinite(xa)) and np.all(np.isfinite(xb))):
        raise ValueError("x_a and x_b must be finite")
    return fit.beta_a * np.maximum(xa - fit.a_a, 0.0), fit.beta_b * np.maximum(xb - fit.a_b, 0.0)


def anchored_pi(x_a, x_b, fit: AnchorFit) -> np.ndarray:
    """``pi = (z_a + 0.5) / (z_a + z_b + 1)``; both z zero gives exactly 0.5. Raises if the fit is not estimable."""
    za, zb = anchored_z(x_a, x_b, fit)
    return (za + SMOOTHING) / (za + zb + 2 * SMOOTHING)
