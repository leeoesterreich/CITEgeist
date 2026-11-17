"""G7 synthetic Moran gates: MSR type-I (definition C), the naive-permutation failure, runtime/memory at n~5000,
and the definition-A spot-level null on synthetic cell-level data.

Independent SAR fields x, y (rho in 0.5/0.8/0.95) on real P4-S2 / P1-S2 lattices and QP-shaped masks.
x variants: continuous, zero-inflated 50%, zero-inflated 80% (point mass at 0, spatially structured
intensity); y is continuous. The test surrogates y only and keeps x fixed.
"""
from __future__ import annotations

import time
from multiprocessing import Pool
from typing import Dict, List, Optional

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import morans as M
import spatial_utils as S
from synthetic_g6 import wilson

ALPHA = 0.05
MAX_TYPE1 = 0.07
LATTICES = {            # name -> (specimen, min cancer proportion)
    "P4S2_full": ("HCC22-088-P4-S2_1i_rep", 0.0),
    "P4S2_ge0.10": ("HCC22-088-P4-S2_1i_rep", 0.10),
    "P4S2_ge0.20": ("HCC22-088-P4-S2_1i_rep", 0.20),
    "P1S2_full": ("HCC22-088-P1-S2", 0.0),
    "P1S2_ge0.10": ("HCC22-088-P1-S2", 0.10),
    "P1S2_ge0.20": ("HCC22-088-P1-S2", 0.20),
}
RHOS = (0.5, 0.8, 0.95)
XTYPES = ("continuous", "zi50", "zi80")
_G: Dict[str, dict] = {}


class Setup:
    def __init__(self, name):
        spec, thr = LATTICES[name]
        self.name = name
        self.xy, _ = S.load_lattice(spec, thr)
        self.pitch = S.median_pitch(self.xy)
        self.W = S.knn_weights(self.xy, 6)
        self.zero_smoother = S.gaussian_smoother(self.xy, 4 * self.pitch)
        self.null = None
        self.build_seconds = None

    def build_null(self):
        t = time.time()
        self.null = M.SpectralNull(self.W)
        self.build_seconds = time.time() - t
        import resource

        self.peak_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0
        self.lu = {r: spla.splu((sp.identity(len(self.xy), format="csc") - r * self.W.tocsc())) for r in RHOS}


def sar_field(st: Setup, rho: float, rng) -> np.ndarray:
    return st.lu[rho].solve(rng.standard_normal(len(st.xy)))


def make_x(st: Setup, rho: float, xtype: str, rng) -> np.ndarray:
    base = np.exp(0.5 * sar_field(st, rho, rng) / 1.0)
    if xtype == "continuous":
        return base
    q = {"zi50": 0.5, "zi80": 0.8}[xtype]
    f = S.smooth_field(st.zero_smoother, rng)
    return np.where(f > np.quantile(f, q), base, 0.0)


def _pair(args):
    name, rho, xtype, seed = args
    st = _G[name]
    rng = np.random.default_rng(seed)
    x = make_x(st, rho, xtype, rng)
    y = sar_field(st, rho, rng)
    out = {}
    for rank in (False, True):
        a = st.null.pvalue(x, y, B=999, seed=42, rank=rank)
        b = M.naive_permutation_pvalue(x, y, st.W, B=999, seed=42, rank=rank)
        tag = "rank" if rank else "raw"
        out[f"msr_{tag}"] = (a["p_two_sided"], a["p_upper"])
        out[f"naive_{tag}"] = (b["p_two_sided"], b["p_upper"])
    return out


def _setting(args):
    name, rho, xtype, n_pairs, seed0 = args
    from threadpoolctl import threadpool_limits

    with threadpool_limits(1):
        res = [_pair((name, rho, xtype, seed0 + i)) for i in range(n_pairs)]
    row = {"lattice": name, "n": len(_G[name].xy), "rho": rho, "x": xtype, "n_pairs": n_pairs}
    for k in res[0]:
        p2 = np.array([r[k][0] for r in res])
        p1 = np.array([r[k][1] for r in res])
        k2, k1 = int((p2 < ALPHA).sum()), int((p1 < ALPHA).sum())
        row[k] = {"type1_two_sided": k2 / n_pairs, "ci95": list(wilson(k2, n_pairs)), "type1_upper": k1 / n_pairs}
    return row


def run_definition_c(n_pairs: int = 1000, procs: int = 8, lattices=None, rhos=RHOS, xtypes=XTYPES, seed0: int = 7000) -> dict:
    lattices = list(lattices or LATTICES)
    timing = {}
    for name in lattices:
        st = Setup(name)
        st.build_null()
        _G[name] = st
        timing[name] = {"n": len(st.xy), "eigh_and_lu_seconds": st.build_seconds, "peak_rss_mb_so_far": st.peak_rss_mb}
    jobs = [(nm, r, xt, n_pairs, seed0 + 100000 * i) for i, (nm, r, xt) in enumerate((n, r, x) for n in lattices for r in rhos for x in xtypes)]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_setting, jobs)
    return {"settings": rows, "timing": timing, "seconds": time.time() - t0, "n_pairs_per_setting": n_pairs, "B": 999, "alpha": ALPHA, "max_type1_allowed": MAX_TYPE1}


def summarize(res: dict) -> dict:
    """Gate verdict. Empty input fails (never passes vacuously)."""
    rows = res["settings"]
    if not rows:
        raise ValueError("no settings: gate cannot pass vacuously")
    out = {"n_settings": len(rows)}
    for k in ("msr_raw", "msr_rank", "naive_raw", "naive_rank"):
        rates = np.array([r[k]["type1_two_sided"] for r in rows])
        out[k] = {"max": float(rates.max()), "mean": float(rates.mean()), "n_over_limit": int((rates > MAX_TYPE1).sum())}
    out["msr_gate_pass"] = bool(out["msr_raw"]["n_over_limit"] == 0)
    hi = [r for r in rows if r["rho"] == max(RHOS)]
    out["naive_inflated_at_high_rho"] = bool(np.mean([r["naive_raw"]["type1_two_sided"] for r in hi]) > MAX_TYPE1)
    out["naive_max_at_high_rho"] = float(max(r["naive_raw"]["type1_two_sided"] for r in hi))
    return out


# --------------------------------------------------------------------------- definition C, non-SAR fields
GAUSS_BW = (3.0, 6.0, 10.0)             # Gaussian-kernel bandwidth in lattice pitches (fields NOT generated by the statistic's W)
_GS: Dict[str, dict] = {}


def _pair_g(args):
    name, bw, xtype, seed = args
    st, sm = _G[name], _GS[(name, bw)]
    rng = np.random.default_rng(seed)
    base = np.exp(0.5 * S.smooth_field(sm, rng))
    if xtype == "continuous":
        x = base
    else:
        q = {"zi50": 0.5, "zi80": 0.8}[xtype]
        f = S.smooth_field(st.zero_smoother, rng)
        x = np.where(f > np.quantile(f, q), base, 0.0)
    y = S.smooth_field(sm, rng)
    a = st.null.pvalue(x, y, B=999, seed=42)
    b = M.naive_permutation_pvalue(x, y, st.W, B=999, seed=42)
    return a["p_two_sided"], b["p_two_sided"]


def _setting_g(args):
    name, bw, xtype, n_pairs, seed0 = args
    from threadpoolctl import threadpool_limits

    with threadpool_limits(1):
        res = np.array([_pair_g((name, bw, xtype, seed0 + i)) for i in range(n_pairs)])
    k = [int((res[:, j] < ALPHA).sum()) for j in (0, 1)]
    return {"lattice": name, "n": len(_G[name].xy), "bw_pitches": bw, "x": xtype, "n_pairs": n_pairs,
            "msr": {"type1_two_sided": k[0] / n_pairs, "ci95": list(wilson(k[0], n_pairs))},
            "naive": {"type1_two_sided": k[1] / n_pairs, "ci95": list(wilson(k[1], n_pairs))}}


def run_definition_c_gauss(n_pairs: int = 1000, procs: int = 8, lattices=("P4S2_ge0.10", "P1S2_ge0.10"), seed0: int = 31000) -> dict:
    for name in lattices:
        st = Setup(name)
        st.build_null()
        _G[name] = st
        for bw in GAUSS_BW:
            _GS[(name, bw)] = S.gaussian_smoother(st.xy, bw * st.pitch)
    jobs = [(nm, bw, xt, n_pairs, seed0 + 100000 * i) for i, (nm, bw, xt) in enumerate((n, b, x) for n in lattices for b in GAUSS_BW for x in XTYPES)]
    t0 = time.time()
    with Pool(procs) as pool:
        rows = pool.map(_setting_g, jobs)
    mx = max(r["msr"]["type1_two_sided"] for r in rows)
    return {"settings": rows, "seconds": time.time() - t0, "n_pairs_per_setting": n_pairs, "B": 999, "alpha": ALPHA,
            "summary": {"msr_max": mx, "msr_pass": bool(mx <= MAX_TYPE1), "naive_min": min(r["naive"]["type1_two_sided"] for r in rows),
                        "n_settings": len(rows)}}


# --------------------------------------------------------------------------- definition A
def make_cells(specimen: str, rng, per_spot_mean: float = 16.8, jitter_frac: float = 0.4):
    xy, _ = S.load_lattice(specimen, 0.0)
    pitch = S.median_pitch(xy)
    m = 1 + rng.poisson(per_spot_mean - 1, size=len(xy))
    spot = np.repeat(np.arange(len(xy)), m)
    ang, rad = rng.uniform(0, 2 * np.pi, len(spot)), jitter_frac * pitch * np.sqrt(rng.uniform(0, 1, len(spot)))
    pos = xy[spot] + np.column_stack([rad * np.cos(ang), rad * np.sin(ang)])
    return xy, spot, pos, pitch


def run_definition_a(specimen: str = "HCC22-088-P1-S2", n_pairs: int = 300, B: int = 999, rhos=RHOS, xtypes=("continuous", "zi80"), structures=("additive", "multiplicative"), seed0: int = 90000) -> dict:
    rng0 = np.random.default_rng(seed0)
    xy_s, spot, pos, pitch = make_cells(specimen, rng0)
    n_cells = len(spot)
    Ws = S.knn_weights(xy_s, 6)
    t = time.time()
    Wc = S.knn_weights(pos, 6)
    knn_seconds = time.time() - t
    t = time.time()
    sn = M.SpectralNull(Ws)
    eigh_seconds = time.time() - t
    sp_null = M.SpotNull(sn, spot, Wc)
    lus = {r: spla.splu(sp.identity(len(xy_s), format="csc") - r * Ws.tocsc()) for r in rhos}
    zero_sm = S.gaussian_smoother(xy_s, 4 * pitch)
    rows = []
    t0 = time.time()
    for rho in rhos:
        for xt in xtypes:
            for struct in structures:
                rej = {"spot_msr": 0, "naive_cell": 0}
                for i in range(n_pairs):
                    rng = np.random.default_rng(seed0 + 17 * i + 1)
                    base = np.exp(0.5 * lus[rho].solve(rng.standard_normal(len(xy_s))))
                    if xt == "continuous":
                        xs = base
                    else:
                        f = S.smooth_field(zero_sm, rng)
                        xs = np.where(f > np.quantile(f, 0.8), base, 0.0)
                    mu = lus[rho].solve(rng.standard_normal(len(xy_s)))
                    x_cell = xs[spot]
                    gamma = rng.standard_normal(3)
                    typ = rng.integers(0, 3, n_cells)
                    if struct == "additive":
                        y_cell = mu[spot] + 0.5 * gamma[typ] + 0.7 * rng.standard_normal(n_cells)
                    else:  # deviation structure scales with the spot field (type-specific slopes)
                        y_cell = mu[spot] * (1.0 + 0.8 * gamma[typ]) + 0.7 * rng.standard_normal(n_cells)
                    a = sp_null.pvalue(x_cell, y_cell, B=B, seed=42)
                    b = M.naive_permutation_pvalue(x_cell, y_cell, Wc, B=199, seed=42)
                    rej["spot_msr"] += a["p_two_sided"] < ALPHA
                    rej["naive_cell"] += b["p_two_sided"] < ALPHA
                rows.append({"rho": rho, "x": xt, "dev_structure": struct, "n_pairs": n_pairs,
                             "spot_msr": {"type1": rej["spot_msr"] / n_pairs, "ci95": list(wilson(int(rej["spot_msr"]), n_pairs))},
                             "naive_cell_permutation": {"type1": rej["naive_cell"] / n_pairs, "ci95": list(wilson(int(rej["naive_cell"]), n_pairs))}})
                print(rows[-1], flush=True)
    return {"specimen_lattice": specimen, "n_spots": len(xy_s), "n_cells": int(n_cells), "knn_cell_seconds": knn_seconds,
            "spot_eigh_seconds": eigh_seconds, "settings": rows, "seconds": time.time() - t0, "B": B, "alpha": ALPHA}
