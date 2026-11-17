"""G6 fix round (rule: g6_fix_round_selection_rule.md): moving-block label-shift null, k_share label, 4 variants."""
from __future__ import annotations

import time
from multiprocessing import Pool

import numpy as np

import hurdle as H
import hurdle_g6_historical as HG
import synthetic_g6 as G

K0_SHARE_FROZEN = 0.0143791
SEED_POOL = 20261015
VARIANTS = {"v_b3_nocov": (3.0, False), "v_b3_cov": (3.0, True), "v_b6_nocov": (6.0, False), "v_b6_cov": (6.0, True)}
CASES = ("null_autocorrelated", "null_autocorrelated_stress", "null_noisy_nst_kshare", "null_depth_confounded")
_CTX = {}


def _init(lat, pool):
    _CTX["lat"], _CTX["pool"] = lat, pool
    G.K0_SHARE = K0_SHARE_FROZEN


def _blocks(lat, side):
    return lat.blocks if side == 6.0 else lat.alt_blocks[side]


def _cov(obs, use_nst):
    return HG.historical_covariate_matrix(obs, drop=() if use_nst else ("log_n_st_eff",))


def _t1(args):
    case, variant, seed, B = args
    lat, pool = _CTX["lat"], _CTX["pool"]
    side, use_nst = VARIANTS[variant]
    rng = np.random.default_rng(seed)
    sim = G.simulate(case, lat, rng)
    Z = _cov(sim["obs"], use_nst)
    L = H.shifted_labels(sim["group"], pool, B, np.random.default_rng(seed + 7))
    r = HG.historical_h1_permutation_test(sim["scores"][0], sim["group"], Z, _blocks(lat, side), B=B, seed=seed, label_matrix=L)
    return {k: (v["p"] if v["status"] == "ok" else None) for k, v in r.items()} | {"frac_high": float(sim["group"].mean())}


def type1(case, variant, n=1000, B=999, seed0=1000, procs=1, lat=None, pool=None):
    G.K0_SHARE = K0_SHARE_FROZEN
    lat = lat or G.make_lattice()
    pool = pool if pool is not None else H.shift_pool(lat.xy, lat.pitch, 4000, SEED_POOL)
    jobs = [(case, variant, seed0 + i, B) for i in range(n)]
    t0 = time.time()
    if procs > 1:
        with Pool(procs, initializer=_init, initargs=(lat, pool)) as p:
            res = p.map(_t1, jobs, chunksize=4)
    else:
        _init(lat, pool)
        res = [_t1(j) for j in jobs]
    out = {"case": case, "variant": variant, "n_datasets": n, "B": B, "seconds": time.time() - t0, "mean_frac_high": float(np.mean([r["frac_high"] for r in res]))}
    for part in ("magnitude", "participation"):
        ps = np.array([r[part] for r in res if r[part] is not None])
        k = int((ps < 0.05).sum())
        lo, hi = G.wilson(k, len(ps))
        out[part] = {"n_ok": int(len(ps)), "n_not_estimable": n - int(len(ps)), "rejections": k, "type1": k / max(len(ps), 1), "ci95": [lo, hi], "upper95": hi}
    return out


def _inj(args):
    seed, B, fold, variant = args
    lat, pool = _CTX["lat"], _CTX["pool"]
    side, use_nst = VARIANTS[variant]
    sim = G.simulate("null_autocorrelated", lat, np.random.default_rng(seed), inject_fold=fold, n_pairs=5)
    Z = _cov(sim["obs"], use_nst)
    L = H.shifted_labels(sim["group"], pool, B, np.random.default_rng(seed + 7))
    ps, gmr0 = [], np.nan
    for j, sc in enumerate(sim["scores"]):
        r = HG.historical_h1_permutation_test(sc, sim["group"], Z, _blocks(lat, side), B=B, seed=seed, label_matrix=L)
        for part in ("magnitude", "participation"):
            ps.append(r[part]["p"] if r[part]["status"] == "ok" else 1.0)
        if j == 0:
            gmr0 = r["magnitude"].get("gmr", np.nan)
    q = HG.historical_benjamini_hochberg(np.array(ps))
    return {"q_ncl_mag": float(q[0]), "gmr": float(gmr0), "q_others_min": float(q[2:].min())}   # q[1] = pair-1 participation (not injected)


def injection(variant, n=50, B=999, fold=3.0, seed0=5000, procs=1, lat=None, pool=None):
    G.K0_SHARE = K0_SHARE_FROZEN
    lat = lat or G.make_lattice()
    pool = pool if pool is not None else H.shift_pool(lat.xy, lat.pitch, 4000, SEED_POOL)
    jobs = [(seed0 + i, B, fold, variant) for i in range(n)]
    if procs > 1:
        with Pool(procs, initializer=_init, initargs=(lat, pool)) as p:
            res = p.map(_inj, jobs)
    else:
        _init(lat, pool)
        res = [_inj(j) for j in jobs]
    gmr = np.array([r["gmr"] for r in res])
    sig = np.array([r["q_ncl_mag"] < 0.05 for r in res])
    acc = np.abs(gmr / fold - 1) <= 0.20
    rec = sig & acc
    fp = np.array([r["q_others_min"] < 0.05 for r in res])
    return {"variant": variant, "n_draws": n, "B": B, "fold": fold, "recovery_rate": float(rec.mean()), "recovery_ci95": list(G.wilson(int(rec.sum()), n)),
            "frac_q_lt_0.05": float(sig.mean()), "frac_within_20pct": float(acc.mean()), "gmr_mean": float(np.nanmean(gmr)), "gmr_sd": float(np.nanstd(gmr)),
            "bh_false_positive_rate_other_pairs": float(fp.mean()), "bh_fp_ci95": list(G.wilson(int(fp.sum()), n))}
