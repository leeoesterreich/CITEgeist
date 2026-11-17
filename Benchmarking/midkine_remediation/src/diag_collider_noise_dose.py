"""Diagnostic (not a gate): noise-dose dependence of the group-effect bias when log(n_st_eff) is a covariate.

Synthetic noisy-n_st null with the k_share label (label independent of the outcome by construction). For
sigma_n in {0.01, 0.4, 0.8}, over N datasets: mean plain-OLS t of the group term among active senders
(log score ~ group + covariates; no block effects), with and without log(n_st_eff) in the covariates.
Reports the dose-response only. It does NOT demonstrate the causal path.

Usage: diag_collider_noise_dose.py [--n 300] [--procs 8] [--out PATH]
"""
import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
import hurdle_g6_historical as HG  # noqa: E402
import provenance  # noqa: E402
import synthetic_g6 as G  # noqa: E402
import synthetic_g6fix as F  # noqa: E402

SIGMAS = (0.01, 0.4, 0.8)
SEED0 = 1000
_LAT = {}


def _t_group(y, g, Z):
    X = np.column_stack([np.ones(len(y)), g, Z])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    res = y - X @ b
    s2 = res @ res / (len(y) - X.shape[1])
    se = np.sqrt(s2 * np.linalg.inv(X.T @ X)[1, 1])
    return float(b[1] / se)


def _one(args):
    sigma, seed = args
    G.K0_SHARE = F.K0_SHARE_FROZEN
    sim = G.simulate(f"diag_sigma_{sigma}", _LAT["lat"], np.random.default_rng(seed))
    sc, g = sim["scores"][0], sim["group"]
    act = sc > 0
    y = np.log(sc[act])
    out = {}
    for name, drop in (("without_log_n_st_eff", ("log_n_st_eff",)), ("with_log_n_st_eff", ())):
        Z = HG.historical_covariate_matrix(sim["obs"], drop=drop)
        out[name] = _t_group(y, g[act], Z[act])
    return out


def run(n=300, procs=8):
    _LAT["lat"] = G.make_lattice()
    for s in SIGMAS:
        G.CASES[f"diag_sigma_{s}"] = dict(G.CASES["null_noisy_nst_kshare"], sigma_n=s)
    res = {}
    with Pool(procs) as p:
        for s in SIGMAS:
            r = p.map(_one, [(s, SEED0 + i) for i in range(n)])
            res[str(s)] = {k: {"mean_t": float(np.mean([x[k] for x in r])), "sd_t": float(np.std([x[k] for x in r]))} for k in r[0]}
    return {"n_datasets": n, "seed0": SEED0, "case": "null_noisy_nst_kshare (sigma_n varied)", "statistic": "plain OLS t of group among actives", "by_sigma_n": res}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--out", type=Path, default=HERE.parent / "synthetic" / "g6_fix" / "diag_collider_noise_dose.json")
    a = ap.parse_args()
    res = run(a.n, a.procs)
    provenance.write_stamped_json(a.out, {"run": "diag_collider_noise_dose", "result": res}, "Benchmarking/midkine_remediation/src/diag_collider_noise_dose.py")
    print(res, flush=True)


if __name__ == "__main__":
    main()
