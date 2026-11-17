"""Run one G6 synthetic gate run and write a stamped JSON (scripts/provenance.py).

Usage: python -u run_g6.py --run <name> [--n 1000] [--procs 8] [--B 4999] [--out DIR]
Runs: t1_<case> (block, B=4999), plain_<case> (free permutation diagnostic), mut_<name>, inj, inj_shuffle.
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))

import provenance  # noqa: E402
import synthetic_g6 as G  # noqa: E402

RUNS = {}
for _c in G.CASES:
    RUNS[f"t1_{_c}"] = ("type1", dict(case=_c, scheme="block"))
    RUNS[f"plain_{_c}"] = ("type1", dict(case=_c, scheme="plain"))
    RUNS[f"sur_{_c}"] = ("type1", dict(case=_c, scheme="surrogate"))
    RUNS[f"b12_{_c}"] = ("type1", dict(case=_c, scheme="blockS12"))
    for _sd in (2, 3, 4):
        RUNS[f"b{_sd}_{_c}"] = ("type1", dict(case=_c, scheme=f"blockS{_sd}"))
RUNS["mut_drop_log_n_st_eff_covariate"] = ("type1", dict(case="null_noisy_nst", scheme="block", mutation="drop_log_n_st_eff_covariate"))
RUNS["mut_drop_depth_covariates"] = ("type1", dict(case="null_depth_confounded", scheme="block", mutation="drop_depth_covariates"))
RUNS["mut_break_block_assignment"] = ("type1", dict(case="null_autocorrelated", scheme="block", mutation="break_block_assignment"))
RUNS["b3m_drop_log_n_st_eff"] = ("type1", dict(case="null_noisy_nst", scheme="blockS3", mutation="drop_log_n_st_eff_covariate"))
RUNS["b3m_drop_depth"] = ("type1", dict(case="null_depth_confounded", scheme="blockS3", mutation="drop_depth_covariates"))
RUNS["b3m_break_blocks"] = ("type1", dict(case="null_autocorrelated", scheme="blockS3", mutation="break_block_assignment"))
RUNS["inj"] = ("injection", dict(mutation="none"))
RUNS["inj_sur"] = ("injection", dict(mutation="none", scheme="surrogate"))
RUNS["inj_shuffle"] = ("injection", dict(mutation="shuffle_group_labels_between_spatial_halves"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, choices=sorted(RUNS))
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--B", type=int, default=None)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--out", type=Path, default=HERE.parent / "synthetic")
    a = ap.parse_args()
    kind, kw = RUNS[a.run]
    kw = dict(kw)
    if kind == "type1":
        kw["n_datasets"] = a.n or 1000
        kw["B"] = a.B or (4999 if kw["scheme"] == "block" and kw.get("mutation", "none") == "none" else 999)
        kw["seed0"] = 1000
        res = G.type1(procs=a.procs, **kw)
    else:
        kw["n_draws"] = a.n or 50
        kw["B"] = a.B or 4999
        res = G.injection(procs=a.procs, **kw)
    a.out.mkdir(parents=True, exist_ok=True)
    path = a.out / f"g6_{a.run}.json"
    provenance.write_stamped_json(path, {"run": a.run, "result": res, "protocol_revision": 4}, "Benchmarking/midkine_remediation/src/run_g6.py")
    print(path, res, flush=True)


if __name__ == "__main__":
    main()
