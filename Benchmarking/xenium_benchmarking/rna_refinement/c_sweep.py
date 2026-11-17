"""Exploratory: can the prior weight c be dropped from the frozen-v2 RNA refinement?

c = 0 is allowed by rna_refinement.py: the prior pseudo-counts are alpha = c * L * p0 = 0, the EM update is
(R + alpha) / (denom + alpha.sum) = R / denom (plain multiplicative EM started from p0), and nothing divides by c or validates c > 0.

Only one setting changes per arm relative to the frozen run (c = 0.1, n_passes = 3, anchor_thr = 0, multinomial, no gene_mask):
  real_c{c}   real n, c in C_GRID
  lib_c{c}    n = library depth (n_s = L_s * median(n_real) / median(L)), c in C_GRID
  c0_passes{k} real n, c = 0, n_passes in PASSES
Every task re-runs the real-n gate (both type sets) and refuses to report if it fails.
  python -u c_sweep.py task [TASK_ID] [OUT_DIR]   (defaults: SLURM_ARRAY_TASK_ID, runs/c_sweep_<SLURM_ARRAY_JOB_ID>)
  python -u c_sweep.py aggregate OUT_DIR
  python c_sweep.py ntasks
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
GATE = {6: 0.7924411383912634, 7: 0.7224664231030145}
GATE_TOL = 1e-9
MS_QUOTED = {6: 0.6496, 7: 0.6037}
C_GRID = (0.0, 0.001, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
PASSES = (1, 2, 3, 5, 10)


def n_library(L, n_real):
    L = np.asarray(L, float)
    return L * (np.median(np.asarray(n_real, float)) / np.median(L))


def task_table():
    """[(family, value)]; family in real_c, lib_c, c0_passes."""
    return [("real_c", c) for c in C_GRID] + [("lib_c", c) for c in C_GRID] + [("c0_passes", k) for k in PASSES]


def arm_name(fam, v):
    return f"{fam}{v:g}"


def mean_abs_diff(a, b):
    return float(np.abs(np.asarray(a, float) - np.asarray(b, float)).mean())


def run_task(task_id, out_dir):
    sys.path.insert(0, str(REPO))
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io, score
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from Benchmarking.xenium_benchmarking.rna_refinement.v2_official import FROZEN, N_REGIONS
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"task_{task_id}.json"
    if out_file.exists():
        raise FileExistsError(out_file)
    fam, val = task_table()[task_id]
    scratch = Path(tempfile.mkdtemp(prefix=f"csweep_{task_id}_"))
    res = {"task_id": task_id, "arm": arm_name(fam, val), "family": fam, "value": val, "frozen_settings": FROZEN, "gate": {}, "types": {}}

    def refine(d, p0, n, types, **over):
        kw = dict(FROZEN) | over
        try:
            p, g = refine_proportions(d["Y"], n, p0, n_bootstrap=0, **kw)
        except ValueError as e:
            return None, {"error": str(e)}
        return p, {
            "fallback": g["fallback"],
            "supported": [bool(x) for x in g["supported_types"]],
            "frozen_types": [t for t, s in zip(types, g["supported_types"]) if not s],
            "em_iters": g["em_iters"],
            "em_converged": g["em_converged"],
            "pass_movement": g["pass_movement"],
            "profile_condition": float(g["profile_condition"]),
        }

    def frame(P, data, types):
        return {r: pd.DataFrame(P[r], index=data[r]["obs"], columns=types) for r in P}

    for n_types in (6, 7):
        types, data = data_io.load_xenium(n_types)
        regs = range(N_REGIONS)
        p0s = {r: rs.ms_share(data[r]["X"], data[r]["names"], rs.PROF[n_types]) for r in regs}
        lib = {r: n_library(data[r]["Y"].sum(1), data[r]["nc"]) for r in regs}
        ms = score.score_official(n_types, frame(p0s, data, types), scratch / "ms")
        base = {r: refine(data[r], p0s[r], data[r]["nc"], types) for r in regs}  # real n, frozen settings = gate
        real = score.score_official(n_types, frame({r: base[r][0] for r in regs}, data, types), scratch / "real")
        gate = abs(real["headline"] - GATE[n_types])
        res["gate"][str(n_types)] = {"got": real["headline"], "want": GATE[n_types], "abs_diff": gate, "ok": bool(gate <= GATE_TOL)}
        print("GATE", n_types, res["gate"][str(n_types)], flush=True)
        if gate > GATE_TOL:
            res["GATE_FAILED"] = True
            out_file.write_text(json.dumps(res, indent=2))
            raise SystemExit(f"STOP: gate failed for {n_types} types")
        r_ms, r_real = ms["headline"], real["headline"]
        lib_base = {r: refine(data[r], p0s[r], lib[r], types) for r in regs}  # library n, c = 0.1
        over = {"c": val} if fam in ("real_c", "lib_c") else {"c": 0.0, "n_passes": val}
        nn = (lambda r: data[r]["nc"]) if fam != "lib_c" else (lambda r: lib[r])
        out, dg, failed = {}, {}, None
        for r in regs:
            out[r], dg[f"region_{r}"] = refine(data[r], p0s[r], nn(r), types, **over)
            if out[r] is None:
                failed = f"region {r}: {dg[f'region_{r}']['error']}"
        run = {"diagnostics": dg}
        if failed:
            run["failed"] = failed
        else:
            s = score.score_official(n_types, frame(out, data, types), scratch / "arm")
            lib_ref = lib_base if fam == "lib_c" else None
            run.update(
                headline=s["headline"], per_type=s["per_type"], per_region=s["per_region"],
                fraction_of_gain_kept=(s["headline"] - r_ms) / (r_real - r_ms),
                mean_abs_diff_from_p0={f"region_{r}": mean_abs_diff(out[r], p0s[r]) for r in regs},
                mean_abs_diff_from_real_n_c0p1={f"region_{r}": mean_abs_diff(out[r], base[r][0]) for r in regs},
                mean_abs_diff_from_same_n_c0p1=({f"region_{r}": mean_abs_diff(out[r], lib_ref[r][0]) for r in regs} if lib_ref else
                                                {f"region_{r}": mean_abs_diff(out[r], base[r][0]) for r in regs}),
            )
        res["types"][str(n_types)] = {"types": types, "reference": {"marker_share": ms, "marker_share_quoted": MS_QUOTED[n_types], "real_n_c0p1": real,
                                                                  "real_n_c0p1_mean_abs_diff_from_p0": {f"region_{r}": mean_abs_diff(base[r][0], p0s[r]) for r in regs},
                                                                  "library_c0p1_mean_abs_diff_from_p0": {f"region_{r}": mean_abs_diff(lib_base[r][0], p0s[r]) for r in regs}},
                                     "run": run}
        print(n_types, res["arm"], run.get("headline"), run.get("failed"), flush=True)
    out_file.write_text(json.dumps(res, indent=2))
    print("wrote", out_file, flush=True)


def _avg(d):
    return float(np.mean(list(d.values())))


def aggregate(out_dir):
    out_dir = Path(out_dir)
    tasks = []
    for i in range(len(task_table())):
        f = out_dir / f"task_{i}.json"
        if not f.exists():
            raise SystemExit(f"STOP: {f} is missing")
        t = json.loads(f.read_text())
        if t.get("GATE_FAILED") or set(t["gate"]) != {"6", "7"} or not all(g["ok"] for g in t["gate"].values()):
            raise SystemExit(f"STOP: task {i} gate failed or incomplete")
        tasks.append(t)
    res = {"gate": tasks[0]["gate"], "frozen_settings": tasks[0]["frozen_settings"], "types": {}}
    lines = []
    for nt in ("6", "7"):
        types = tasks[0]["types"][nt]["types"]
        ref = tasks[0]["types"][nt]["reference"]
        arms = {t["arm"]: t["types"][nt]["run"] | {"family": t["family"], "value": t["value"]} for t in tasks}
        res["types"][nt] = {"types": types, "reference": ref, "arms": arms}
        lines += [f"## {nt}-type", "",
                  f"marker share {ref['marker_share']['headline']:.4f}; real n c=0.1 {ref['real_n_c0p1']['headline']:.4f}; "
                  f"mean|p-p0| at c=0.1: real n {_avg(ref['real_n_c0p1_mean_abs_diff_from_p0']):.4f}, library n {_avg(ref['library_c0p1_mean_abs_diff_from_p0']):.4f}", "",
                  "| arm | headline | gain kept | mean abs p-p0 | mean abs vs same-n c=0.1 | mean abs vs real-n c=0.1 | fallback | frozen | EM not conv | " + " | ".join(types) + " |",
                  "|---|---|---|---|---|---|---|---|---|" + "---|" * len(types)]
        for name, a in arms.items():
            if "headline" not in a:
                lines.append(f"| {name} | FAILED: {a['failed']} |")
                continue
            dg = a["diagnostics"].values()
            lines.append(
                f"| {name} | {a['headline']:.4f} | {a['fraction_of_gain_kept']:.2f} | {_avg(a['mean_abs_diff_from_p0']):.4f} | "
                f"{_avg(a['mean_abs_diff_from_same_n_c0p1']):.4f} | {_avg(a['mean_abs_diff_from_real_n_c0p1']):.4f} | "
                f"{sum(1 for d in dg if d['fallback'])} | {sum(len(d['frozen_types']) for d in dg)} | "
                f"{sum(1 for d in dg for c in d['em_converged'] if not c)} | " + " | ".join(f"{a['per_type'][k]:.3f}" for k in types) + " |")
        lines.append("")
    out = out_dir / "results.json"
    if out.exists():
        raise FileExistsError(out)
    out.write_text(json.dumps(res, indent=2))
    (out_dir / "summary.md").write_text("\n".join(lines))
    print("\n".join(lines))
    print("wrote", out, flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "task"
    if cmd == "ntasks":
        print(len(task_table()))
    elif cmd == "aggregate":
        aggregate(sys.argv[2])
    elif cmd == "task":
        tid = int(sys.argv[2]) if len(sys.argv) > 2 else int(os.environ["SLURM_ARRAY_TASK_ID"])
        od = sys.argv[3] if len(sys.argv) > 3 else HERE / "runs" / f"c_sweep_{os.environ['SLURM_ARRAY_JOB_ID']}"
        run_task(tid, od)
    else:
        raise SystemExit(f"unknown command {cmd}")
