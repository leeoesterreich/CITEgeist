"""v42 stage X12: EM convergence check of the default (frozen v2) refinement, stamped.

Question (convergence_DECISION_RULE.md, committed before any number): is the reported result set by the EM iteration cap (max_iter = 150,
tol = 1e-5)? The computation is ``convergence_check.run_task``, called unchanged for all 12 (arm, max_iter) tasks (same inputs, settings,
regions, types, gate). This stage only runs them in one job, checks them, and writes ONE stamped ``convergence.json`` that holds the per-task
results and the summary the manuscript builder reads (tol, cap, max_iters_to_converge, r_cap_{6,7}, r_converged_{6,7}, max_abs_change).

Gate (fail closed; a miss writes gate_failure.json and exits 1, nothing is stamped, _COMPLETE.json is not written):
  * every task passed run_task's own gate (the cap-150 frozen run reproduces the c_sweep constants within GATE_TOL);
  * the frozen arm c0.1_p3 has exactly the caps 150 / 1000 / 10000 and every pass converged at 10000;
  * the frozen cap-150 headline equals the stamped X2 per_type_mean_r within TOL (1e-6), expected value read from the X2 file;
  * the frozen cap-150 predictions, refit here, equal the stamped X1 predictions within PRED_TOL (1e-6), expected values read from the X1 files;
  * the package defaults are tol 1e-5 and max_iter 150 (read from the function signature, then written into the result).

Code runs from a clean detached worktree (.commands/submit-x12-convergence.sh); data and outputs stay in the main checkout ($V2_NAMESPACE_ROOT).
CLI (CITEgeist_env python, inside an sbatch job):
  run [--out DIR] [--stamp]   all 12 tasks into DIR/tasks, gate, write DIR/convergence.json + DIR/_COMPLETE.json; --stamp stamps last.
                              Default DIR = <namespace root>/.../runs/<X12 id>.
  paths                       shell assignments (CONV_OUT, RID_X12) for the submit script
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (str(REPO), str(XB), str(REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from Benchmarking.xenium_benchmarking.rna_refinement import convergence_check as cc  # noqa: E402
from Benchmarking.xenium_benchmarking.rna_refinement import v2_official as vo  # noqa: E402

FROZEN_ARM = "c0.1_p3"
EXPECTED_CAPS = {150, 1000, 10000}
TOL = vo.HEADLINE_TOL  # 1e-6: frozen cap-150 headline vs the stamped X2
PRED_TOL = 1e-6
EM_TOL = 1e-5
EM_CAP = 150
OUT_NAME = "convergence.json"
REL_RUNS = "Benchmarking/xenium_benchmarking/rna_refinement/runs"


def namespace() -> Path:
    import v42_stamp as W

    return W.namespace_root(REPO)


def out_dir() -> Path:
    """Output dir inside the namespace root (main checkout); its name IS the X12 run id."""
    return namespace() / REL_RUNS / vo.run_id("X12")


def package_defaults() -> dict:
    from CITEgeist.model.deconvolution.rna_refinement import em_refine, refine_proportions

    # refine_proportions passes no tol to em_refine, so em_refine's default is the tolerance in force
    return {"tol": inspect.signature(em_refine).parameters["tol"].default, "cap": inspect.signature(refine_proportions).parameters["max_iter"].default}


def summarize(tasks: list, x2: dict, defaults: dict, x1_diff: dict, tol: float = TOL) -> dict:
    """Gate the 12 task results and build the summary. Raises ValueError (problems joined) on any miss; never returns a partial result."""
    problems = []
    if [t.get("task_id") for t in tasks] != list(range(len(cc.task_table()))):
        problems.append(f"task ids {[t.get('task_id') for t in tasks]} != 0..{len(cc.task_table()) - 1}")
    for t in tasks:
        if set(t.get("gate", {})) != {"6", "7"} or not all(g["ok"] for g in t["gate"].values()):
            problems.append(f"task {t.get('task_id')} ({t.get('arm')}, {t.get('max_iter')}): run_task gate not ok")
    if defaults != {"tol": EM_TOL, "cap": EM_CAP}:
        problems.append(f"package defaults {defaults} != tol {EM_TOL}, cap {EM_CAP}")
    runs = {t["max_iter"]: t for t in tasks if t.get("arm") == FROZEN_ARM}
    if set(runs) != EXPECTED_CAPS or len([t for t in tasks if t.get("arm") == FROZEN_ARM]) != len(EXPECTED_CAPS):
        problems.append(f"{FROZEN_ARM}: expected caps {sorted(EXPECTED_CAPS)}, got {sorted(runs)}")
    if problems:
        raise ValueError("; ".join(problems))
    big, base = runs[10000], runs[EM_CAP]
    diags = [dg for nt in ("6", "7") for dg in big["types"][nt]["diagnostics"].values()]
    if not diags or not all(all(dg["em_converged"]) for dg in diags):
        problems.append(f"{FROZEN_ARM} did not converge at max_iter 10000")
    for nt in ("6", "7"):
        want = float(x2[nt]["per_type_mean_r"])
        got = float(base["types"][nt]["headline"])
        if abs(got - want) > tol:
            problems.append(f"{nt}t frozen cap-150 headline {got!r} vs stamped X2 {want!r}: |diff| {abs(got - want):.3g} > {tol:g}")
    for nt in ("6", "7"):
        if not x1_diff[nt] <= PRED_TOL:  # written so that NaN fails
            problems.append(f"{nt}t frozen predictions vs stamped X1: max abs diff {x1_diff[nt]!r} > {PRED_TOL:g}")
    if problems:
        raise ValueError("; ".join(problems))
    out = {"tol": defaults["tol"], "cap": defaults["cap"], "max_iters_to_converge": max(max(dg["em_iters"]) for dg in diags)}
    for nt in ("6", "7"):
        out[f"r_cap_{nt}"] = base["types"][nt]["headline"]
        out[f"r_converged_{nt}"] = big["types"][nt]["headline"]
    out["max_abs_change"] = max(abs(out[f"r_converged_{nt}"] - out[f"r_cap_{nt}"]) for nt in ("6", "7"))
    out["verdict"] = "CAP_LOAD_BEARING" if out["max_abs_change"] > cc.TOL else "CAP_NOT_LOAD_BEARING"  # all converged was checked above
    return out


def read_x2(ns: Path, n: int) -> dict:
    f = ns / vo._rel(vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json")
    return json.loads(f.read_text())["summary"]


def x1_max_abs(ns: Path) -> dict:
    """Max abs difference of the frozen default fit (vo.compute_frozen: cap 150) vs the stamped X1 predictions, per type set."""
    import pandas as pd

    out = {}
    for n in (6, 7):
        fit = vo.compute_frozen(n)
        out[str(n)] = max(
            vo.max_abs_diff_frames(df, pd.read_csv(ns / vo._rel(vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"), index_col=0))
            for r, df in fit.items())
    return out


def stamp_inputs(ns: Path) -> list:
    """Frozen-v2 model inputs (h5ad), the stamped X2 results it is gated against, and the stamped X1 predictions (the upstream the frozen arm reproduces)."""
    rel = lambda p: ns / vo._rel(p)  # noqa: E731
    ins = [rel(vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json") for n in (6, 7)]
    ins += [rel(vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz") for n in (6, 7) for r in range(vo.N_REGIONS)]
    ins += [rel(f) for r in range(vo.N_REGIONS) for f in vo.x1_inputs(r)]
    return ins


def run(out: Path, stamp: bool) -> dict:
    import v42_stamp as W

    out = Path(out)
    ns = namespace()
    vo.refuse_existing(out)
    art = out / OUT_NAME
    started = W.utc_now()
    if stamp:
        W.preflight(vo._rel(__file__), out_paths=[art], repo_root=ns, git_root=REPO)  # clean tree, tracked producer, SLURM job: before the compute
    t0 = time.time()
    out.mkdir(parents=True)
    tdir = out / "tasks"
    for tid in range(len(cc.task_table())):
        cc.run_task(tid, tdir)
    files = [tdir / f"task_{i}.json" for i in range(len(cc.task_table()))]
    tasks = [json.loads(f.read_text()) for f in files]
    x1 = x1_max_abs(ns)
    sha = os.environ.get("SHA") or W._git(REPO, "rev-parse", "HEAD").strip()
    try:
        summary = summarize(tasks, {str(n): read_x2(ns, n) for n in (6, 7)}, package_defaults(), x1)
    except ValueError as e:
        (out / "gate_failure.json").write_text(json.dumps({"status": "GATE_FAILED", "problem": str(e), "sha": sha}, indent=2))
        raise SystemExit(f"STOP: {e}")
    result = {
        "stage": "X12", "sha": sha, "frozen_arm": FROZEN_ARM,
        "note": "EM convergence check of the default refinement (convergence_DECISION_RULE.md): same fits at max_iter 150 / 1000 / 10000. Exploratory; "
                "no official output or manuscript number changes.",
        **summary,
        "tasks": tasks,
        "task_files_sha256": {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
        "gate": {"tol": TOL, "pred_tol": PRED_TOL, "x1_frozen_predictions_max_abs": x1, "x2_per_type_mean_r": {n: float(read_x2(ns, int(n))["per_type_mean_r"]) for n in ("6", "7")}},
    }
    art.write_text(json.dumps(result, indent=2, sort_keys=True))
    if stamp:
        W.write_stamp(art, run_id=vo.run_id("X12"), protocol_sha12="none", inputs=stamp_inputs(ns), producer=vo._rel(__file__), started_at=started,
                      repo_root=ns, git_root=REPO, stage_role="x12_convergence")
    (out / "_COMPLETE.json").write_text(json.dumps(
        {"status": "COMPLETE", "stage": "X12", "stamped": bool(stamp), "sha": sha, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
         "seconds": round(time.time() - t0, 1), "verdict": summary["verdict"]}, indent=2))
    print("X12 done", out, summary, flush=True)
    return result


def cmd_run(args) -> None:
    run(Path(args.out) if args.out else out_dir(), args.stamp)


def cmd_paths(_args) -> None:
    print(f"CONV_OUT={out_dir()}")
    print(f"RID_X12={vo.run_id('X12')}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    pr = sub.add_parser("run")
    pr.add_argument("--out", default=None)
    pr.add_argument("--stamp", action="store_true")
    pr.set_defaults(fn=cmd_run)
    sub.add_parser("paths").set_defaults(fn=cmd_paths)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
