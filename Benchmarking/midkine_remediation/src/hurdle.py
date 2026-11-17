"""Hurdle building blocks for the DESCRIPTIVE H1 analysis (protocol rev 6; spec 9.4.7).

Spatial blocks, the primary covariate set, the estimability gates and the moving-block label-shift pool used by the
negative control. The estimator itself is ``descriptive.py``. The permutation tests, BH and surrogate nulls of
revisions 1-4 live in ``hurdle_g6_historical.py`` (G6 record only; no p-value exists in the primary API).
"""
from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Dict

import numpy as np

MIN_MIXED_BLOCKS = 10
MIN_PER_GROUP = 30
MIN_ACTIVE_PER_GROUP = 10
BLOCK_SIDE_PITCHES = 6.0
BLOCK_MIN_N = 8
PRIMARY_COVARIATES = ("log_spot_umi", "log_rna_per_cell_equiv", "qp_prop")      # protocol hurdle.models.covariates
PRIMARY_COVARIATES_TOTAL = ("log_spot_umi", "log_rna_total", "qp_prop")           # protocol deviation 2026-11-06 (rna="total"): column 1 = log of the raw allocated-count total
OPTIONAL_COVARIATE = "log_n_st_eff"                                               # explicit opt-in only (collider; A12.3 fix round)
RNA_MODES = ("per_cell", "total")
DEVIATION_FILE = "PROTOCOL_DEVIATION_2026-11-06_rna_total.md"                     # beside protocol.json (Benchmarking/midkine_remediation/)
DEVIATION_ID = "rna_total_2026-11-06"


class CovariateDomainError(ValueError):
    """A log-transformed covariate is not strictly positive / finite (typed not_estimable upstream, never floored)."""


# --------------------------------------------------------------------------- blocks and covariates
def make_blocks(coords: np.ndarray, pitch: float, side: float = 6.0, min_n: int = 8) -> np.ndarray:
    """Square tiles of ``side`` pitches from coordinates only; tiles with < ``min_n`` members merge into the nearest tile."""
    xy = np.asarray(coords, float)
    cell = side * pitch
    ij = np.floor((xy - xy.min(axis=0)) / cell).astype(int)
    _, lab = np.unique(ij[:, 0] * (ij[:, 1].max() + 1) + ij[:, 1], return_inverse=True)
    lab = lab.copy()
    while True:
        ids, counts = np.unique(lab, return_counts=True)
        small = ids[counts < min_n]
        if len(small) == 0 or len(ids) == 1:
            break
        s = small[np.argmin(counts[counts < min_n])]
        cents = {i: xy[lab == i].mean(axis=0) for i in ids}
        others = [i for i in ids if i != s]
        d = [np.linalg.norm(cents[i] - cents[s]) for i in others]
        lab[lab == s] = others[int(np.argmin(d))]
    _, lab = np.unique(lab, return_inverse=True)
    return lab


def median_pitch(coords: np.ndarray) -> float:
    from scipy.spatial import cKDTree

    d, _ = cKDTree(coords).query(coords, k=2)
    return float(np.median(d[:, 1]))


def standardise(Z: np.ndarray) -> np.ndarray:
    Z = np.asarray(Z, float)
    sd = Z.std(axis=0)
    sd[sd == 0] = 1.0
    return (Z - Z.mean(axis=0)) / sd


def primary_covariates(rna: str = "per_cell") -> tuple:
    """The three primary covariate names for ``rna`` ('per_cell': protocol; 'total': protocol deviation 2026-11-06)."""
    if rna not in RNA_MODES:
        raise ValueError(f"rna must be one of {RNA_MODES}, got {rna!r}")
    return PRIMARY_COVARIATES if rna == "per_cell" else PRIMARY_COVARIATES_TOTAL


def rna_total(obs) -> np.ndarray:
    """rna_total = rna_per_cell_equiv * n_st_eff * depth_ratio per compartment: the raw allocated-count total over the genes not in keratins + MDK.

    Exact for the stored primary columns (it equals the float64 row sum of E over those genes to about 1e-16, the float32 layer E within 1e-6 relative). It has no
    nuclei term, so the nuclei-free run (R11) uses the same value. The one definition for every stage (``d5_fixed_moran.rna_total_gate`` checks it against layer E)."""
    return obs["rna_per_cell_equiv"].to_numpy(float) * obs["n_st_eff"].to_numpy(float) * obs["depth_ratio"].to_numpy(float)


def deviation_record() -> dict:
    """{id, path (repo-relative), sha256} of the committed protocol deviation file that licenses ``rna="total"``. Raises FileNotFoundError if it is missing."""
    f = Path(__file__).resolve().parents[1] / DEVIATION_FILE
    return {"id": DEVIATION_ID, "path": f"Benchmarking/midkine_remediation/{DEVIATION_FILE}", "sha256": hashlib.sha256(f.read_bytes()).hexdigest()}


def covariate_matrix(obs, include_log_n_st_eff: bool = False, drop: tuple = (), nuclei_free: bool = False, rna: str = "per_cell") -> np.ndarray:
    """Primary covariates, standardised: log(spot_umi), the RNA covariate, qp_prop.

    * ``rna="per_cell"`` (default, the frozen protocol; old runs stay byte-reproducible): column 1 = log(rna_per_cell_equiv). ``rna="total"`` (protocol deviation
      2026-11-06, ``PROTOCOL_DEVIATION_2026-11-06_rna_total.md``): column 1 = log(``rna_total(obs)``) and nothing else changes. Every production caller passes ``rna=`` explicitly.
    * ``log(n_st_eff)`` is added ONLY by ``include_log_n_st_eff=True`` (explicit opt-in; it is not a primary covariate).
    * ``drop`` names columns to omit (mutation / sensitivity use); the names of ``primary_covariates(rna)``.
    * ``nuclei_free=True`` (R11) swaps ``rna_per_cell_equiv`` for ``rna_per_cell_equiv_nf`` and nothing else when ``rna="per_cell"``; with ``rna="total"`` it changes nothing (the raw total has no nuclei term).
    * Frozen domain rule (replaces any floor): every log-transformed covariate must be finite and strictly > 0, otherwise
      ``CovariateDomainError`` (the caller records a typed not_estimable). qp_prop must be finite.
    * Frozen standardisation population: the rows passed in (the specimen's eligible cancer compartments), population SD
      (ddof=0), computed ONCE on the full set, never per bootstrap resample or per subset.
    """
    names = primary_covariates(rna)
    unknown = [d for d in drop if d not in names + ((OPTIONAL_COVARIATE,) if include_log_n_st_eff else ())]
    if unknown:                                                    # a name of the other mode (e.g. log_rna_per_cell_equiv with rna="total") would silently drop nothing
        raise ValueError(f"drop names {unknown} are not columns of the rna={rna!r} matrix {names}")
    rna_col = "rna_per_cell_equiv_nf" if nuclei_free else "rna_per_cell_equiv"
    raw = {"log_spot_umi": "spot_umi", "log_rna_per_cell_equiv": rna_col, "log_n_st_eff": "n_st_eff"}     # column names; read lazily
    cols = {}
    for k in names + ((OPTIONAL_COVARIATE,) if include_log_n_st_eff else ()):
        if k == "qp_prop":
            v = obs["qp_prop"].to_numpy(float)
            if not np.isfinite(v).all():
                raise CovariateDomainError("qp_prop is not finite")
            cols[k] = v
        else:
            v = rna_total(obs) if k == "log_rna_total" else obs[raw[k]].to_numpy(float)
            if not (np.isfinite(v).all() and (v > 0).all()):
                raise CovariateDomainError(f"{k}: values must be finite and strictly positive before the log")
            cols[k] = np.log(v)
    keep = [k for k in cols if k not in drop]
    return standardise(np.column_stack([cols[k] for k in keep])) if keep else np.zeros((len(obs), 0))


# --------------------------------------------------------------------------- groups, mixed blocks, estimability
def _groups(block: np.ndarray):
    order = np.argsort(block, kind="stable")
    b = block[order]
    cuts = np.flatnonzero(np.diff(b)) + 1
    return np.split(order, cuts)


def n_mixed_blocks(block, group) -> int:
    """Blocks containing at least one compartment of each group."""
    g = np.asarray(group).astype(int)
    return int(sum(1 for idx in _groups(np.asarray(block)) if 0 < int(g[idx].sum()) < len(idx)))


def estimability(group, block, kind, active=None) -> Dict:
    """Typed gate. Returns {ok, reason_code, reason, n_mixed_blocks}. No permutation-count requirement (rev 5)."""
    g = np.asarray(group).astype(int)
    blk = np.asarray(block)
    if kind == "magnitude":
        sel = np.asarray(active, bool)
        g, blk = g[sel], blk[sel]
        n1, n0 = int(g.sum()), int((1 - g).sum())
        if min(n1, n0) < MIN_ACTIVE_PER_GROUP:
            return {"ok": False, "reason_code": "active_senders", "reason": f"active senders per group {n1}/{n0} < {MIN_ACTIVE_PER_GROUP}"}
    else:
        n1, n0 = int(g.sum()), int((1 - g).sum())
        if min(n1, n0) < MIN_PER_GROUP:
            return {"ok": False, "reason_code": "group_size", "reason": f"compartments per group {n1}/{n0} < {MIN_PER_GROUP}"}
    mixed = n_mixed_blocks(blk, g)
    if mixed < MIN_MIXED_BLOCKS:
        return {"ok": False, "reason_code": "mixed_blocks", "reason": f"{mixed} mixed-label blocks < {MIN_MIXED_BLOCKS}", "n_mixed_blocks": mixed}
    return {"ok": True, "reason_code": "", "reason": "", "n_mixed_blocks": mixed}


def _orth_basis(X):
    """Orthonormal basis of the column space (rank-revealing SVD; falls back to the slower gesvd driver)."""
    try:
        U, s, _ = np.linalg.svd(X, full_matrices=False)
    except np.linalg.LinAlgError:
        from scipy.linalg import svd as _svd

        U, s, _ = _svd(X, full_matrices=False, lapack_driver="gesvd")
    keep = s > 1e-9 * s.max()
    return U[:, keep], s




def _irls(X, a, offset_cols=None, iters=60, ridge=1e-8):
    """Logistic regression by Newton/IRLS. Returns beta, mu, converged."""
    beta = np.zeros(X.shape[1])
    for it in range(iters):
        eta = X @ beta
        mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
        w = mu * (1 - mu)
        H = X.T @ (w[:, None] * X) + ridge * np.eye(X.shape[1])
        step = np.linalg.solve(H, X.T @ (a - mu))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-9:
            eta = X @ beta
            mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
            return beta, mu, True
    eta = X @ beta
    return beta, 1 / (1 + np.exp(-np.clip(eta, -30, 30))), False




# --------------------------------------------------------------------------- between-block null: moving-block label shifts
class ShiftPoolUnfillable(RuntimeError):
    """The acceptance rule yielded fewer than the requested number of shifts (a data limitation, e.g. few compartments)."""

    def __init__(self, n_found: int, n_requested: int, acceptance: dict = None):
        super().__init__(f"only {n_found} acceptable shifts found")
        self.n_found, self.n_requested, self.acceptance = int(n_found), int(n_requested), dict(acceptance or {})


def shift_pool(xy, pitch: float, n_shifts: int = 4000, seed: int = 20261015, min_shift_pitch: float = 6.0, min_cover: float = 0.5, cover_radius_pitch: float = 1.5, max_tries: int = 200000) -> np.ndarray:
    """(n_shifts, n) int32 source indices: row s maps compartment i to the compartment whose label it receives under shift s.

    Coordinates are translated by a random vector, wrapped on a torus over the bounding box, and matched to the nearest
    observed compartment. A shift is accepted if its length (torus) is >= ``min_shift_pitch`` pitches and >= ``min_cover``
    of compartments land within ``cover_radius_pitch`` pitches of an observed compartment (uncovered ones take the nearest).
    """
    from scipy.spatial import cKDTree

    xy = np.asarray(xy, float)
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = hi - lo + pitch
    tree = cKDTree(xy)
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(max_tries):
        d = rng.uniform(0, 1, 2) * span
        dt = np.minimum(d, span - d)                      # torus displacement
        if np.hypot(*dt) < min_shift_pitch * pitch:
            continue
        q = lo + np.mod(xy - lo + d, span)
        dist, idx = tree.query(q)
        if (dist <= cover_radius_pitch * pitch).mean() >= min_cover:
            out.append(idx.astype(np.int32))
            if len(out) == n_shifts:
                break
    if len(out) < n_shifts:
        raise ShiftPoolUnfillable(len(out), n_shifts, {"min_shift_pitch": min_shift_pitch, "min_cover": min_cover,
                                                       "cover_radius_pitch": cover_radius_pitch, "max_tries": max_tries})
    return np.array(out)


def shifted_labels(group, pool: np.ndarray, B: int, rng: np.random.Generator) -> np.ndarray:
    """(B, n) int8 label matrix: observed labels moved by B shifts sampled from ``pool``."""
    sel = rng.choice(len(pool), size=B, replace=False)
    return np.asarray(group)[pool[sel]].astype(np.int8)
