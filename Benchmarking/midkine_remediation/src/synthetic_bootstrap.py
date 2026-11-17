"""Coverage of the block-bootstrap intervals on the G6 null and injected cases (descriptive-path check)."""
import time
from multiprocessing import Pool

import numpy as np

import block_bootstrap as BB
import hurdle_g6_historical as HG
import synthetic_g6 as G
import synthetic_g6fix as F

_C = {}


def _init(lat):
    _C["lat"] = lat
    G.K0_SHARE = F.K0_SHARE_FROZEN


def _one(args):
    case, seed, side, B, fold = args
    lat = _C["lat"]
    sim = G.simulate(case, lat, np.random.default_rng(seed), inject_fold=fold)
    blocks = lat.blocks if side == 6.0 else lat.alt_blocks[side]
    Z = HG.historical_covariate_matrix(sim["obs"], drop=("log_n_st_eff",))
    r = BB.block_bootstrap(sim["scores"][0], sim["group"], Z, blocks, B=B, seed=seed)
    est = HG.historical_h1_permutation_test(sim["scores"][0], sim["group"], Z, blocks, B=19, seed=seed)
    return {"gmr": est["magnitude"].get("gmr", np.nan), "gmr_ci": (r["gmr"]["lo"], r["gmr"]["hi"]),
            "or": est["participation"].get("odds_ratio", np.nan), "or_ci": (r["odds_ratio"]["lo"], r["odds_ratio"]["hi"])}


def coverage(case, side=6.0, n=300, B=500, fold=1.0, seed0=1000, procs=1, lat=None):
    lat = lat or G.make_lattice()
    G.K0_SHARE = F.K0_SHARE_FROZEN
    jobs = [(case, seed0 + i, side, B, fold) for i in range(n)]
    t0 = time.time()
    if procs > 1:
        with Pool(procs, initializer=_init, initargs=(lat,)) as p:
            res = p.map(_one, jobs)
    else:
        _init(lat)
        res = [_one(j) for j in jobs]
    out = {"case": case, "tile_side_pitches": side, "n_datasets": n, "n_boot": B, "fold": fold, "seconds": time.time() - t0}
    truth = {"gmr": fold, "or": 1.0}
    for k, ck in (("gmr", "gmr_ci"), ("or", "or_ci")):
        ok = [(r[ck][0], r[ck][1]) for r in res if np.isfinite(r[ck][0]) and np.isfinite(r[ck][1])]
        cov = np.mean([lo <= truth[k] <= hi for lo, hi in ok]) if ok else float("nan")
        out[k] = {"n_ok": len(ok), "coverage": float(cov), "ci95_wilson": list(G.wilson(int(round(cov * len(ok))), len(ok))) if ok else None,
                  "median_width_ratio": float(np.median([hi / lo for lo, hi in ok if lo > 0])) if ok else None}
    return out
