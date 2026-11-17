"""v42 stage X11: the leave-one-region-out (LORO) tuning audit of the frozen v2 Xenium arm, recomputed from frozen-v2 inputs.

What it is: frozen v2 is ONE setting (thr 0, c 0.1, 3 passes) fixed before any scoring. The audit asks what a tuner would have
gained: for each held-out region k, take the setting (of the 72-setting grid) with the best per-type mean r on the OTHER four regions
(GT is read only to rank settings, never as a model input), score region k with it, pool the five regions with the official strict
scorer. Stored result (variants.json, candidate a_V3_multinomial, git_head 65d435c3): 6-type 0.7906, 7-type 0.7373, beside the frozen
0.7924 / 0.7225. This stage recomputes it and stamps it, so the number cited beside the frozen v2 headline has stamped lineage.
Settings were chosen on this benchmark (sensitivity_DECISION_RULE.md); this audit changes none of them.

Inputs are the frozen-v2 inputs of X1: protein (CITE h5ad), raw spot counts (GEX h5ad) and h5ad obs["n_cells"]. The fits and the scoring are the
experiment's own code, called not copied: ``run_variants`` (candidate a, ``_job``, ``fold_scores``), ``reproduce_spike`` (marker share, GT
frames), ``score`` (``loro_select``, ``score_official``).

Gate (fail closed; a miss writes gate_failure.json and exits 1, nothing is stamped, _COMPLETE.json is not written):
  * the recomputed frozen and LORO rows (headline, every per-type r, every per-region value) equal the stored r8 row within TOL (1e-6, the tolerance
    X1/X2 use against the same r8 file) and the LORO setting chosen for each held-out region is IDENTICAL (exact list equality);
  * the number of EM fallbacks of the grid equals the stored count;
  * the frozen row equals the official X2 per_type_mean_r (headline and per type) within TOL, and the frozen predictions equal the stamped
    X1 predictions within PRED_TOL (1e-9).
Expected values come from the stored r8 file and the stamped X1/X2 files, never from the recompute.

CLI (CITEgeist_env python, repo root, inside an sbatch job):
  run [--out DIR] [--stamp]   compute, gate, write DIR/loro_audit.json + DIR/_COMPLETE.json; --stamp (the official run: SHA env = the commit that
                              holds this code) proves the code tree first and stamps loro_audit.json as stage X11 last. Default DIR = runs/<X11 id>.
  paths                       shell assignments (LORO_OUT, RID_X11) for the submit script
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (str(REPO), str(XB), str(REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from Benchmarking.xenium_benchmarking.rna_refinement import v2_official as vo  # noqa: E402

N_REGIONS = 5
N_SETTINGS = 72  # reproduce_spike.GRID: 4 thr x 6 c x 3 rounds (fixed constant, not derived from the grid being checked)
CANDIDATE = "a_V3_multinomial"
DEFAULT_ROW = "package default (frozen)"
LORO_ROW = "LORO tuning audit"
FROZEN_SETTING = (0.0, 0.1, 3)
TOL = vo.HEADLINE_TOL  # 1e-6: recomputed rows vs the stored r8 row and vs X2
PRED_TOL = 1e-9  # recomputed frozen predictions vs the stamped X1 predictions (the prediction tolerance of v2_stress / v2_heldout)
REFERENCE = vo.VARIANTS_JSON  # runs/variants_20261031_1404_r8/variants.json
OUT_NAME = "loro_audit.json"


def out_dir() -> Path:
    """The official output dir: its name IS the X11 run id (the lineage path rules need a v2 run-name component)."""
    return HERE / "runs" / vo.run_id("X11")


# ----------------------------------------------------------------------------- compute
def _imports():
    from joblib import Parallel, delayed

    from Benchmarking.xenium_benchmarking.rna_refinement import data_io, score
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from Benchmarking.xenium_benchmarking.rna_refinement import run_variants as rv

    return Parallel, delayed, data_io, score, rs, rv


def fold_margins(fold: dict) -> dict:
    """Per held-out region: best minus second-best mean over the other regions (a near tie is where thread count could flip an argmax)."""
    regions = sorted(fold)
    out = {}
    for k in regions:
        means = sorted((float(np.mean([fold[r][s] for r in regions if r != k])) for s in fold[regions[0]]), reverse=True)
        out[str(k)] = means[0] - means[1]
    return out


def compute_type(n: int, scored_root: Path, n_jobs: int) -> dict:
    """Frozen row, LORO row, the LORO selections and the frozen predictions for n types. Candidate a_V3_multinomial of run_variants."""
    Parallel, delayed, data_io, score, rs, rv = _imports()
    model, grid, frozen = rv.CANDIDATES[CANDIDATE]
    if tuple(frozen) != FROZEN_SETTING or len(grid) != N_SETTINGS:
        raise SystemExit(f"STOP: candidate {CANDIDATE} is not the frozen setting {FROZEN_SETTING} over {N_SETTINGS} settings")
    types, data = data_io.load_xenium(n)
    p0s = {r: rs.ms_share(data[r]["X"], data[r]["names"], rs.PROF[n]) for r in range(N_REGIONS)}
    gts = rs._gt_frames(n, types, data)
    settings = list(grid)
    res = Parallel(n_jobs=n_jobs)(delayed(rv._job)(r, data[r], p0s[r], model, s, None) for r in range(N_REGIONS) for s in settings)
    preds = {(r, s): P for r, s, P, *_ in res}
    n_fallback = sum(1 for _r, _s, _P, fb, _c in res if fb)
    frozen_fb = [r for r, s, _P, fb, _c in res if fb and tuple(s) == FROZEN_SETTING]
    if frozen_fb:
        raise SystemExit(f"STOP: the frozen setting fell back in regions {frozen_fb} ({n} types): not a frozen-v2 result")
    frozen_preds = {r: rv._frame(preds[(r, FROZEN_SETTING)], data, r, types) for r in range(N_REGIONS)}
    m_def = score.score_official(n, frozen_preds, scored_root / f"default_{n}type")
    fold = rv.fold_scores(gts, preds, settings)
    sel = score.loro_select(fold)
    loro_preds = {r: rv._frame(preds[(r, sel[r])], data, r, types) for r in range(N_REGIONS)}
    m_loro = score.score_official(n, loro_preds, scored_root / f"loro_{n}type")
    return {
        "frozen": m_def, "loro": m_loro, "selections": {str(k): list(v) for k, v in sel.items()}, "fallbacks": n_fallback,
        "fold_margin_best_minus_second": fold_margins(fold), "n_settings": len(settings),
        "fold_scores": {str(r): {",".join(str(x) for x in s): v for s, v in d.items()} for r, d in fold.items()},
        "_frozen_preds": frozen_preds,
    }


# ----------------------------------------------------------------------------- gate
def _cmp(label, got, want, tol, out):
    d = abs(float(got) - float(want))
    out.append({"what": label, "got": float(got), "want": float(want), "abs_diff": d, "tol": tol, "ok": bool(d <= tol)})


def _compare_row(prefix, got: dict, want: dict, tol, out, problems):
    """headline, every per-type r and every per-region value of two scored rows; a differing key set or length is a problem, not a skip."""
    _cmp(f"{prefix} headline", got["headline"], want["headline"], tol, out)
    if set(got["per_type"]) != set(want["per_type"]):
        problems.append(f"{prefix}: per-type keys differ: {sorted(got['per_type'])} vs {sorted(want['per_type'])}")
    for t in sorted(set(got["per_type"]) & set(want["per_type"])):
        _cmp(f"{prefix} per_type {t}", got["per_type"][t], want["per_type"][t], tol, out)
    if len(got["per_region"]) != len(want["per_region"]) or len(want["per_region"]) != N_REGIONS:
        problems.append(f"{prefix}: per-region length {len(got['per_region'])} vs stored {len(want['per_region'])} (expected {N_REGIONS})")
    for i, (g, w) in enumerate(zip(got["per_region"], want["per_region"])):
        _cmp(f"{prefix} per_region {i}", g, w, tol, out)


def gate(computed: dict, stored: dict, x2: dict, x1_max_abs: dict, tol: float = TOL, pred_tol: float = PRED_TOL) -> dict:
    """computed: {n: compute_type result}; stored: the r8 variants.json dict; x2: {n: official X2 summary dict}; x1_max_abs: {n: max abs diff of the
    frozen predictions vs X1}. Returns {"problems": [...], "comparisons": [...], "max_abs_diff": float}; ok iff problems is empty and every comparison ok."""
    problems, comps = [], []
    for n in (6, 7):
        c = computed[n]
        for row, key in ((DEFAULT_ROW, "frozen"), (LORO_ROW, "loro")):
            hit = [r for r in stored["rows"] if r["variant"] == CANDIDATE and r["row"] == row and r["n_types"] == n]
            if len(hit) != 1:
                problems.append(f"stored r8 file has {len(hit)} rows for {CANDIDATE} / {row} / {n}t (expected exactly 1)")
                continue
            _compare_row(f"{n}t {row} vs r8", c[key], hit[0], tol, comps, problems)
        want_sel = stored["loro_selections"].get(f"{CANDIDATE}/{n}")
        if want_sel is None:
            problems.append(f"stored r8 file has no loro_selections for {CANDIDATE}/{n}")
        elif c["selections"] != {str(k): list(v) for k, v in want_sel.items()}:
            problems.append(f"{n}t LORO selections differ from r8: got {c['selections']}, stored {want_sel}")
        want_fb = stored["fallbacks"].get(f"{CANDIDATE}/{n}")
        if want_fb is None or c["fallbacks"] != want_fb:
            problems.append(f"{n}t EM fallbacks: got {c['fallbacks']}, stored {want_fb}")
        # frozen row vs the official X2 (per_type_mean_r and every {type}_mean_r)
        s2 = x2[n]
        _cmp(f"{n}t frozen vs X2 per_type_mean_r", c["frozen"]["headline"], s2["per_type_mean_r"], tol, comps)
        for t, v in c["frozen"]["per_type"].items():
            if f"{t}_mean_r" not in s2:
                problems.append(f"{n}t X2 summary has no {t}_mean_r")
            else:
                _cmp(f"{n}t frozen vs X2 {t}_mean_r", v, s2[f"{t}_mean_r"], tol, comps)
        _cmp(f"{n}t frozen predictions vs stamped X1 (max abs)", x1_max_abs[n], 0.0, pred_tol, comps)
    problems += [f"{x['what']}: got {x['got']!r}, stored {x['want']!r}, |diff| {x['abs_diff']:.3g} > {x['tol']:g}" for x in comps if not x["ok"]]
    return {"problems": problems, "comparisons": comps, "max_abs_diff": max([x["abs_diff"] for x in comps] or [float("nan")])}


def read_x2(n: int) -> dict:
    f = vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
    return json.loads(f.read_text())["summary"]


def x1_diff(n: int, frozen_preds: dict) -> float:
    worst = 0.0
    for r, df in frozen_preds.items():
        f = vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
        worst = max(worst, vo.max_abs_diff_frames(df, pd.read_csv(f, index_col=0)))
    return worst


# ----------------------------------------------------------------------------- run
def _sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stamp_inputs() -> list:
    """The inputs of the audit: the frozen-v2 model inputs (h5ad), the stamped X1 predictions and X2 results it is gated against, and the stored r8 row.
    The GT files are not listed (as in the X2 stamp: the scorer reads them, the stamp lists the predictions)."""
    ins = [REFERENCE] + [vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json" for n in (6, 7)]
    ins += [vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz" for n in (6, 7) for r in range(N_REGIONS)]
    ins += [f for r in range(N_REGIONS) for f in vo.x1_inputs(r)]
    return ins


def run(out: Path, stamp: bool, n_jobs: int) -> dict:
    out = Path(out)
    vo.refuse_existing(out)
    started = None
    if stamp:
        import v42_stamp as W

        started = W.utc_now()
        vo.assert_code_tree()  # before the compute: a stamped run on uncommitted code is refused up front (the stamp step proves it again)
    t0 = time.time()
    stored = json.loads(REFERENCE.read_text())
    out.mkdir(parents=True)
    computed = {n: compute_type(n, out / "scored", n_jobs) for n in (6, 7)}
    x1_max_abs = {n: x1_diff(n, computed[n]["_frozen_preds"]) for n in (6, 7)}
    g = gate(computed, stored, {n: read_x2(n) for n in (6, 7)}, x1_max_abs)
    head = {n: (computed[n]["frozen"]["headline"], computed[n]["loro"]["headline"]) for n in (6, 7)}
    print("frozen/LORO headline", head, "max |diff| vs references", g["max_abs_diff"], "problems", len(g["problems"]), flush=True)
    sha = os.environ.get("SHA") or subprocess.check_output(["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True).strip()
    if g["problems"]:
        (out / "gate_failure.json").write_text(json.dumps(
            {"status": "GATE_FAILED", "problems": g["problems"], "comparisons": g["comparisons"], "tol": TOL, "pred_tol": PRED_TOL,
             "recomputed_headlines": {str(n): {"frozen": head[n][0], "loro": head[n][1]} for n in (6, 7)},
             "recomputed_selections": {str(n): computed[n]["selections"] for n in (6, 7)}, "reference": str(REFERENCE.relative_to(REPO)), "sha": sha},
            indent=2))
        raise SystemExit("STOP: the recompute does not reproduce the stored r8 row / X2 / X1: " + "; ".join(g["problems"][:6]))
    ref = json.loads(REFERENCE.read_text())
    result = {
        "stage": "X11", "candidate": CANDIDATE, "sha": sha, "n_settings": N_SETTINGS, "frozen_setting": {"thr": 0.0, "c": 0.1, "rounds": 3},
        "note": "LORO tuning audit of the frozen v2 arm: per held-out region the setting is chosen on the other four regions (GT ranks settings only). "
                "Settings were chosen on this benchmark; the audit changes none. The frozen v2 headline (X2) is unchanged.",
        "headline": {str(n): {"frozen_default": head[n][0], "loro_audit": head[n][1]} for n in (6, 7)},
        "frozen_default": {str(n): {k: computed[n]["frozen"][k] for k in ("headline", "per_type", "per_region")} for n in (6, 7)},
        "loro_audit": {str(n): {k: computed[n]["loro"][k] for k in ("headline", "per_type", "per_region")} for n in (6, 7)},
        "loro_selections": {str(n): computed[n]["selections"] for n in (6, 7)},
        "fold_margin_best_minus_second": {str(n): computed[n]["fold_margin_best_minus_second"] for n in (6, 7)},
        "fallbacks": {str(n): computed[n]["fallbacks"] for n in (6, 7)},
        "fold_scores": {str(n): computed[n]["fold_scores"] for n in (6, 7)},
        "gate": {"tol": TOL, "pred_tol": PRED_TOL, "max_abs_diff": g["max_abs_diff"], "n_comparisons": len(g["comparisons"]),
                 "reference": str(REFERENCE.relative_to(REPO)), "reference_sha256": _sha256(REFERENCE), "reference_git_head": ref["git_head"],
                 "x1_frozen_predictions_max_abs": {str(n): x1_max_abs[n] for n in (6, 7)}},
    }
    art = out / OUT_NAME
    art.write_text(json.dumps(result, indent=2, sort_keys=True))
    if stamp:
        vo.stamp_artifact(art, run_id=vo.run_id("X11"), inputs=stamp_inputs(), producer=vo._rel(__file__), started_at=started, role="x11_loro_audit")
    (out / "_COMPLETE.json").write_text(json.dumps(
        {"status": "COMPLETE", "stage": "X11", "stamped": bool(stamp), "sha": sha, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
         "seconds": round(time.time() - t0, 1), "max_abs_diff_vs_references": g["max_abs_diff"]}, indent=2))
    print("X11 done", out, flush=True)
    return result


def cmd_run(args) -> None:
    n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
    run(Path(args.out) if args.out else out_dir(), args.stamp, n_jobs)


def cmd_paths(_args) -> None:
    print(f"LORO_OUT={out_dir()}")
    print(f"RID_X11={vo.run_id('X11')}")


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
