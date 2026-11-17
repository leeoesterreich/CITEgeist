"""Reproduce the RNA-refinement spike headline from committed code.

Run:  python -u reproduce_spike.py OUT_DIR   (OUT_DIR must be new)
Writes OUT_DIR/metrics.json and OUT_DIR/loro_selections.json (+ prediction CSVs).
Model inputs: protein, raw spot counts, h5ad n_cells. GT is used only by score.py.
"""
import hashlib
import itertools
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / "rna_refinement_spike"))
sys.path.insert(0, str(REPO))

import lib  # noqa: E402  (spike reference code; its `load` is replaced below)
from Benchmarking.xenium_benchmarking.rna_refinement import data_io, score  # noqa: E402

lib.load = data_io.load_xenium
learn_mu, em_fit, run_variant, ms_share, PROF = lib.learn_mu, lib.em_fit, lib.run_variant, lib.ms_share, lib.PROF

THRS = (0.0, 0.4, 0.6, 0.8)
CS = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0)
ROUNDS = (1, 2, 3)
GRID = list(itertools.product(THRS, CS, ROUNDS))  # (thr, c, rounds); frozen default = (0.0, 0.1, 3)
FROZEN = {"V1": (0.0, 0.1, 1), "V3": (0.0, 0.1, 3)}

# Reference values reproduced with h5ad obs["n_cells"] nuclei (official scorer), tolerance +-0.002.
EXPECTED = {
    "6type": {"V0": 0.6496, "V1": 0.7649, "V3": 0.7924, "LORO_V3": 0.7906},
    "7type": {"V0": 0.6037, "V1": 0.7184, "V3": 0.7225, "LORO_V3": 0.7373},
}
TOL = 0.002
NOTES = {
    "nuclei_source": "h5ad obs['n_cells'] (GT n_cells is never used)",
    "diagnostic_gt_ncells_numbers_only": {
        "note": "Spike originally used GT-CSV n_cells; that differs from h5ad n_cells in ~90% of spots (max abs diff 22-35). "
        "Diagnostic rerun with it (not committed): V3 frozen 6t/7t 0.7931/0.7273; LORO V3 6t/7t 0.7931/0.7327.",
    },
}


def _job(r, d, p0, s):
    return r, s, run_variant(d, p0, s[0], s[1], s[2], False)


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _frame(P, data, r, types):
    return pd.DataFrame(P, index=data[r]["obs"], columns=types)


def _nanmean_r(G, P):
    return float(np.nanmean(lib.score_region(G, P)))


def main(out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    metrics = dict(
        params=dict(grid=dict(thr=THRS, c=CS, rounds=ROUNDS), frozen=FROZEN, use_yield=False, n_jobs=8),
        git_head=subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip(),
        versions={m.__name__: m.__version__ for m in (np, pd, scipy, sc)},
        input_sha256={Path(data_io.h5ad_path(r, k)).name: _sha(data_io.h5ad_path(r, k)) for r in range(5) for k in ("CITE", "GEX")},
        results={},
    )
    selections = {}
    for n in (6, 7):
        types, data = data_io.load_xenium(n)
        p0s = {r: ms_share(data[r]["X"], data[r]["names"], PROF[n]) for r in range(5)}
        res = Parallel(n_jobs=8)(delayed(_job)(r, data[r], p0s[r], s) for r in range(5) for s in GRID)
        preds = {(r, s): P for r, s, P in res}

        rows = {"V0": {r: _frame(p0s[r], data, r, types) for r in range(5)}}
        for v, s in FROZEN.items():
            rows[v] = {r: _frame(preds[(r, s)], data, r, types) for r in range(5)}
        # LORO audit: setting for held-out k chosen by mean per-type r on the other 4 regions.
        # GT is read here ONLY to rank settings on non-held-out regions (tuning audit), never as a model input.
        gts = _gt_frames(n, types, data)
        fold = {r: {s: _nanmean_r(gts[r][0], preds[(r, s)][gts[r][1]]) for s in GRID} for r in range(5)}
        sel = score.loro_select(fold)
        selections[str(n)] = {str(k): list(v) for k, v in sel.items()}
        rows["LORO_V3"] = {r: _frame(preds[(r, sel[r])], data, r, types) for r in range(5)}

        metrics["results"][f"{n}type"] = {
            name: score.score_official(n, pr, out_dir / f"preds_{n}type_{name}") for name, pr in rows.items()
        }
        for name, m in metrics["results"][f"{n}type"].items():
            print(n, name, round(m["headline"], 4), flush=True)
    metrics["notes"] = NOTES
    bad = {
        f"{k}/{v}": (metrics["results"][k][v]["headline"], e)
        for k, d in EXPECTED.items()
        for v, e in d.items()
        if abs(metrics["results"][k][v]["headline"] - e) > TOL
    }
    metrics["expected_check"] = dict(tol=TOL, expected=EXPECTED, out_of_tolerance=bad)
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2, default=str))
    (out_dir / "loro_selections.json").write_text(json.dumps(selections, indent=2))
    if bad:
        raise SystemExit(f"reproduction out of tolerance: {bad}")
    return metrics


def _gt_frames(n, types, data):
    """Per region: (GT proportions array, row positions of GT spots in the prediction)."""
    gtd = Path(score._context(n)[0])
    out = {}
    for r in range(5):
        g = pd.read_csv(gtd / "ground_truth" / f"Xenium_region_{r}_prop.csv", index_col=0)
        out[r] = (g[types].values, data[r]["obs"].get_indexer(g.index))
    return out


if __name__ == "__main__":
    main(sys.argv[1])
