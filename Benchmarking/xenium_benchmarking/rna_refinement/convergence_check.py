"""Exploratory: is the frozen v3 refinement result set by the EM iteration cap? Rule: convergence_DECISION_RULE.md.

One array task per (arm, max_iter); inputs, scorer and gate are those of c_sweep.py.
  python -u convergence_check.py task [TASK_ID] [OUT_DIR]
  python -u convergence_check.py aggregate OUT_DIR
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO))
from Benchmarking.xenium_benchmarking.rna_refinement.c_sweep import GATE, GATE_TOL, mean_abs_diff  # noqa: E402

ARMS = (("c0.1_p3", 0.1, 3), ("c0.03_p3", 0.03, 3), ("c0.3_p3", 0.3, 3), ("c0_p1", 0.0, 1))
CAPS = (150, 1000, 10000)
TOL = 0.01


def task_table():
    return [(a, c, p, m) for a, c, p in ARMS for m in CAPS]


def run_task(tid, out_dir):
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io, score
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from Benchmarking.xenium_benchmarking.rna_refinement.v2_official import FROZEN, N_REGIONS
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    f = out_dir / f"task_{tid}.json"
    if f.exists():
        raise FileExistsError(f)
    arm, c, passes, cap = task_table()[tid]
    scratch = Path(tempfile.mkdtemp(prefix=f"conv_{tid}_"))
    res = {"task_id": tid, "arm": arm, "c": c, "n_passes": passes, "max_iter": cap, "gate": {}, "types": {}}
    for nt in (6, 7):
        types, data = data_io.load_xenium(nt)
        regs = range(N_REGIONS)
        p0s = {r: rs.ms_share(data[r]["X"], data[r]["names"], rs.PROF[nt]) for r in regs}

        def fit(r, **kw):
            return refine_proportions(data[r]["Y"], data[r]["nc"], p0s[r], n_bootstrap=0, **(dict(FROZEN) | kw))

        def sc(P, tag):
            return score.score_official(nt, {r: pd.DataFrame(P[r], index=data[r]["obs"], columns=types) for r in regs}, scratch / tag)

        gate = sc({r: fit(r)[0] for r in regs}, "gate")["headline"]
        ok = abs(gate - GATE[nt]) <= GATE_TOL
        res["gate"][str(nt)] = {"got": gate, "want": GATE[nt], "ok": bool(ok)}
        if not ok:
            f.write_text(json.dumps(res, indent=2))
            raise SystemExit(f"STOP: gate failed {nt}")
        ref = {r: fit(r, c=c, n_passes=passes)[0] for r in regs}  # same arm at max_iter = 150
        out, dg = {}, {}
        for r in regs:
            out[r], g = fit(r, c=c, n_passes=passes, max_iter=cap)
            dg[f"region_{r}"] = {"em_iters": g["em_iters"], "em_converged": g["em_converged"], "fallback": g["fallback"],
                                 "mean_abs_diff_from_cap150": mean_abs_diff(out[r], ref[r])}
        s = sc(out, "arm")
        res["types"][str(nt)] = {"headline": s["headline"], "per_type": s["per_type"], "diagnostics": dg,
                                 "headline_cap150": sc(ref, "ref")["headline"]}
        print(nt, arm, cap, s["headline"], flush=True)
    f.write_text(json.dumps(res, indent=2))


def aggregate(out_dir):
    out_dir = Path(out_dir)
    T = [json.loads((out_dir / f"task_{i}.json").read_text()) for i in range(len(task_table()))]
    lines = ["| arm | max_iter | 6t | 7t | all passes converged | max EM iters | mean abs p - p(cap 150) |", "|---|---|---|---|---|---|---|"]
    for t in T:
        dgs = [d for nt in ("6", "7") for d in t["types"][nt]["diagnostics"].values()]
        conv = all(all(d["em_converged"]) for d in dgs)
        lines.append(f"| {t['arm']} | {t['max_iter']} | {t['types']['6']['headline']:.4f} | {t['types']['7']['headline']:.4f} | {conv} | "
                     f"{max(max(d['em_iters']) for d in dgs)} | {np.mean([d['mean_abs_diff_from_cap150'] for d in dgs]):.4f} |")
    fz = {t["max_iter"]: t for t in T if t["arm"] == "c0.1_p3"}
    big = fz[10000]
    moved = max(abs(big["types"][nt]["headline"] - fz[150]["types"][nt]["headline"]) for nt in ("6", "7"))
    conv = all(all(d["em_converged"]) for nt in ("6", "7") for d in big["types"][nt]["diagnostics"].values())
    verdict = "CAP_LOAD_BEARING" if moved > TOL else ("NOT_CONVERGED" if not conv else "CAP_NOT_LOAD_BEARING")
    lines += ["", f"Verdict: **{verdict}** (frozen c = 0.1 / 3 passes: max headline change at max_iter 10000 = {moved:.4f}; all converged = {conv})"]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "task"
    if cmd == "aggregate":
        aggregate(sys.argv[2])
    else:
        tid = int(sys.argv[2]) if len(sys.argv) > 2 else int(os.environ["SLURM_ARRAY_TASK_ID"])
        od = sys.argv[3] if len(sys.argv) > 3 else HERE / "runs" / f"convergence_{os.environ['SLURM_ARRAY_JOB_ID']}"
        run_task(tid, od)
