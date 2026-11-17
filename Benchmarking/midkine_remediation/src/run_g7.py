"""Run the G7 synthetic Moran gates and write stamped JSON.

Usage: python -u run_g7.py --run {c|a|bench} [--n N] [--procs P] [--out DIR] [--lattices a,b]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))

import provenance  # noqa: E402
import synthetic_g7 as G  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, choices=["c", "a", "bench", "cg"])
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--out", type=Path, default=HERE.parent / "synthetic")
    ap.add_argument("--lattices", default=None)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.run == "c":
        lat = a.lattices.split(",") if a.lattices else None
        res = G.run_definition_c(n_pairs=a.n or 1000, procs=a.procs, lattices=lat)
        res["summary"] = G.summarize(res)
        name = "g7_definition_c.json"
    elif a.run == "cg":
        res = G.run_definition_c_gauss(n_pairs=a.n or 1000, procs=a.procs)
        name = "g7_definition_c_gauss.json"
    elif a.run == "a":
        res = G.run_definition_a(n_pairs=a.n or 300)
        name = "g7_definition_a.json"
    else:
        import numpy as np
        import morans as M
        import spatial_utils as S

        res = {}
        for spec, thr in (("HCC22-088-P1-S2", 0.0), ("HCC22-088-P4-S2_1i_rep", 0.0)):
            xy, _ = S.load_lattice(spec, thr)
            res[spec] = M.benchmark_eigh(S.knn_weights(xy, 6))
        name = "g7_runtime_n5000.json"
    path = a.out / name
    provenance.write_stamped_json(path, {"result": res, "protocol_revision": 4}, "Benchmarking/midkine_remediation/src/run_g7.py")
    print(path, res.get("summary", ""), flush=True)


if __name__ == "__main__":
    main()
