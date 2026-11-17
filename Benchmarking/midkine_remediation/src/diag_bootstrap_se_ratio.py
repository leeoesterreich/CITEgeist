"""Diagnostic (not a gate): empirical SD of the point estimate vs the block-bootstrap SE, on synthetic nulls.

Per dataset: point estimate (h1_test, as in synthetic_bootstrap) and block-bootstrap 95% interval
(6-pitch tiles, primary covariates, no log(n_st_eff)). Bootstrap SE on the log scale = (log hi - log lo) / (2 * 1.96).
Reports, per case and per outcome (GMR, OR): SD of the log point estimate over datasets, median and mean bootstrap SE,
and ratio = SD / median SE (>1 means the bootstrap SE is too small by that factor).
Seeds: dataset i uses seed SEED0 + i (SEED0 = 1000), also as its block-bootstrap seed (same as run_bootstrap.py).

Usage: diag_bootstrap_se_ratio.py [--n 300] [--B 500] [--procs 8] [--out PATH]
"""
import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
import synthetic_bootstrap as SB  # noqa: E402
import provenance  # noqa: E402

CASES = ("null_autocorrelated", "null_autocorrelated_stress")
SEED0 = 1000
SIDE = 6.0
Z95 = 1.959964


def ratio_table(res):
    """res: list of _one() dicts -> per-outcome SD / bootstrap-SE summary (log scale)."""
    out = {}
    for k, ck in (("gmr", "gmr_ci"), ("or", "or_ci")):
        pt, se = [], []
        for r in res:
            lo, hi = r[ck]
            if np.isfinite(r[k]) and np.isfinite(lo) and np.isfinite(hi) and lo > 0 and hi > 0 and r[k] > 0:
                pt.append(np.log(r[k]))
                se.append((np.log(hi) - np.log(lo)) / (2 * Z95))
        pt, se = np.array(pt), np.array(se)
        sd = float(np.std(pt, ddof=1))
        out[k] = {"n_ok": int(len(pt)), "sd_log_point_estimate": sd, "mean_log_point_estimate": float(np.mean(pt)),
                  "median_bootstrap_se": float(np.median(se)), "mean_bootstrap_se": float(np.mean(se)),
                  "ratio_sd_over_median_se": sd / float(np.median(se))}
    return out


def run(n, B, procs):
    lat = SB.G.make_lattice()
    out = {}
    for case in CASES:
        SB.G.K0_SHARE = SB.F.K0_SHARE_FROZEN
        jobs = [(case, SEED0 + i, SIDE, B, 1.0) for i in range(n)]
        if procs > 1:
            with Pool(procs, initializer=SB._init, initargs=(lat,)) as p:
                res = p.map(SB._one, jobs)
        else:
            SB._init(lat)
            res = [SB._one(j) for j in jobs]
        out[case] = ratio_table(res)
        print(case, out[case], flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--B", type=int, default=500)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--out", type=Path, default=HERE.parent / "synthetic" / "g6_fix" / "diag_bootstrap_se_ratio.json")
    a = ap.parse_args()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    res = run(a.n, a.B, a.procs)
    payload = {"settings": {"n_datasets": a.n, "n_boot": a.B, "seed0": SEED0, "tile_side_pitches": SIDE, "cases": list(CASES),
                            "se_definition": "(log hi - log lo) / (2 * 1.959964) per dataset"}, "result": res}
    provenance.write_stamped_json(a.out, payload, "Benchmarking/midkine_remediation/src/diag_bootstrap_se_ratio.py")


if __name__ == "__main__":
    main()
