"""v42 plan Task 9: X8 (settings table) and the withheld-gene (RNA-reuse) test of frozen v2 under the fair GEX standard.

Commands (CITEgeist_env python, repo root, inside an sbatch job; SHA proves the code tree as in v2_official):
  x8                     settings table from runs/sensitivity_20261101_1510 (a re-aggregation, no refit), v42 stamp as stage X8
  wh-region OUT R        fold fits + SACE no_round for region R into the NEW run dir OUT (rule: heldout_genes_DECISION_RULE.md)
  wh-score OUT           fair scoring on the fold genes and the mechanical verdict -> OUT/heldout_genes.json (NOT v42-stamped: v42_run_ids.yaml has no
                         stage id for this analysis; the file records the commit, the rule hash and git blob, the roots it read and the sha256 of
                         every input). Exit 0 only for HELD_OUT_HOLDS; HELD_OUT_WORSE exits 3 and NO_VERDICT exits 4, after the files are written.
  Env: SHA = the FIRST stage's commit (region 0's record); WH_RERUN=1 moves a started region aside and reruns it.
X7 varies the PROTEIN input (v2_stress.py); X8 varies the SETTINGS (thr, c, passes) and nothing else.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
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
K_FOLDS = 5
N_SETTINGS = 72  # reproduce_spike.GRID: 4 thr x 6 c x 3 rounds (fixed constant)
SENS_JSON = HERE / "runs" / "sensitivity_20261101_1510" / "sensitivity.json"
CELL_COUNT_SACE = HERE / "runs" / "cell_count_20261101_1424" / "sace"
OLD_FOLD_PROPS = HERE / "runs" / "gexgate_20261031_1832" / "props"
X8_OUT_REL = "Benchmarking/xenium_benchmarking/evaluation/rnaprot_v2_settings"
FOLD_TOL = 1e-8
FULL_TOL = 1e-9
REF_TOL = 1e-6
COMPETITORS = ("Cell2Location", "Tangram")
DECISION_RULE = HERE / "heldout_genes_DECISION_RULE.md"
# Fixed constants (rule file: "fewer than 5 folds x 5 regions x 7 types blocks"): never derived from the outputs being checked.
N_CELL_TYPES = 7
N_GENES = 405
BLOCKS_PER_TYPE = K_FOLDS * N_REGIONS  # 25 per (method, type)
EXPECTED_BLOCKS = K_FOLDS * N_REGIONS * N_CELL_TYPES  # 175 per method
METHODS = ("crossfit", "v2_full", "v41", "Cell2Location", "Tangram")
CONTROL_TOL = 1e-9  # fold SACE path given the full-fit proportions vs the X3 layers: max abs difference of any layer value (measured 0.0 exactly, regions 0 and 3)
EXIT_WORSE, EXIT_NO_VERDICT = 3, 4


def x8_dir() -> Path:
    return REPO / X8_OUT_REL / vo.run_id("X8")


# ----------------------------------------------------------------------------- X8
def build_settings_table(sens: dict, x2_refs: dict, n_settings: int = N_SETTINGS) -> dict:
    """The compact settings table of the v2 manuscript from the sensitivity run. Raises unless the run is complete (72 settings), its
    frozen row reproduces the official X2 per_type_mean_r within 1e-6 at 6 and 7 types, and its own parity record is ok.

    The v41 comparator is dropped on purpose (a v41 number does not belong in a v42 table); RCTD and C2L stay."""
    table = sens["table"]
    if sens.get("n_settings") != n_settings or len(table) != n_settings:
        raise ValueError(f"sensitivity run has {len(table)} settings (n_settings {sens.get('n_settings')}), expected {n_settings}")
    if not sens["frozen_parity"]["ok"]:
        raise ValueError("the sensitivity run's own frozen-parity record is not ok")
    frozen = sens["frozen_reference"]
    bad = {n: (frozen[str(n)], x2_refs[n]) for n in (6, 7) if abs(frozen[str(n)] - x2_refs[n]) > REF_TOL}
    if bad:
        raise ValueError(f"frozen row differs from the official X2 per_type_mean_r beyond {REF_TOL}: {bad}")
    row = next((r for r in table if (r["thr"], r["c"], r["rounds"]) == (0.0, 0.1, 3)), None)
    if row is None:
        raise ValueError("no frozen (thr 0, c 0.1, rounds 3) row in the table")
    for n in (6, 7):
        if abs(row[f"r{n}"] - x2_refs[n]) > REF_TOL:
            raise ValueError(f"table frozen row r{n} differs from X2")
    comps = {k: v for k, v in sens["comparators"].items() if k in ("RCTD", "C2L")}
    if set(comps) != {"RCTD", "C2L"}:
        raise ValueError("comparators RCTD and C2L are required")
    return {
        "stage": "X8", "source_run": "runs/sensitivity_20261101_1510/sensitivity.json", "rule": sens["rule"], "n_settings": n_settings,
        "frozen_setting": {"thr": 0.0, "c": 0.1, "rounds": 3}, "frozen_reference_x2": {str(n): x2_refs[n] for n in (6, 7)},
        "comparators": comps, "table": table, "one_at_a_time": sens["one_at_a_time"], "marginal": sens["marginal"], "summary": sens["summary"],
        "frac_above_RCTD_and_C2L_both_types": sens["frac_above_RCTD_and_C2L_both_types"], "bray_curtis": sens["bray_curtis"],
        "a_priori_point": sens["a_priori_point"], "verdict": sens["verdict"], "gex": sens["gex"],
        "note": "re-aggregation of the existing sensitivity run (no refit); X7 (protein-input stress test) is a different analysis",
    }


def cmd_x8(_args) -> None:
    import v42_stamp as W

    started = W.utc_now()
    vo.assert_code_tree()
    out = x8_dir()
    vo.refuse_existing(out)
    refs = {}
    for n in (6, 7):
        f = vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
        refs[n] = float(json.loads(f.read_text())["summary"]["per_type_mean_r"])
    res = build_settings_table(json.loads(SENS_JSON.read_text()), refs)
    out.mkdir(parents=True)
    js = out / "settings_table.json"
    js.write_text(json.dumps(res, indent=2, sort_keys=True))
    cs = out / "settings_table.csv"
    with open(cs, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["thr", "c", "rounds", "r6", "r7", "bray_curtis_7t"])
        for r in res["table"]:
            w.writerow([r["thr"], r["c"], r["rounds"], r["r6"], r["r7"], r["bray_curtis_7t"]])
    ins = [SENS_JSON] + [vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json" for n in (6, 7)]
    ins += [vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz" for n in (6, 7) for r in range(N_REGIONS)]
    for art, role in ((js, "x8_settings_table"), (cs, "x8_settings_table_csv")):
        vo.stamp_artifact(art, run_id=vo.run_id("X8"), inputs=ins, producer=vo._rel(__file__), started_at=started, role=role)
    print("X8 done", out, "verdict", res["verdict"]["verdict"], flush=True)


# ----------------------------------------------------------------------------- withheld genes: pure logic
def fold_genes_upper(genes, folds) -> list:
    """Uppercase gene-name set per fold (the scorer works on uppercase names)."""
    return [{str(genes[i]).upper() for i in f} for f in folds]


def restrict(df: pd.DataFrame, genes: set) -> pd.DataFrame:
    return df[[c for c in df.columns if c in genes]]


def verdict(per_set: dict, block_counts: dict, regressions_ok: bool, expected_blocks: int = EXPECTED_BLOCKS, methods=METHODS,
            block_problems=()) -> dict:
    """per_set: {"S2": {method: {type: rmse}}}; block_counts: {method: number of scored (finite-rmse) S2 blocks}.
    The rule of heldout_genes_DECISION_RULE.md, mechanically (set S2). NO_VERDICT unless EVERY method has exactly the expected number
    of scored blocks (a method averaged over fewer regions or folds would be compared with one averaged over all of them)."""
    if not regressions_ok:
        return {"verdict": "NO_VERDICT", "reason": "regression check failed"}
    off = {m: block_counts.get(m, 0) for m in methods if block_counts.get(m, 0) != expected_blocks}
    if off:
        return {"verdict": "NO_VERDICT", "reason": f"scored blocks != {expected_blocks} for {off}"}
    if block_problems:
        return {"verdict": "NO_VERDICT", "reason": f"block completeness: {list(block_problems)[:5]}"}

    def n_below_both(m):
        mine = per_set["S2"][m]
        n = 0
        for t, v in mine.items():
            others = [per_set["S2"][k].get(t, np.nan) for k in COMPETITORS]
            n += int(bool(np.isfinite([v] + others).all() and all(v < o for o in others)))
        return n

    nb = {m: n_below_both(m) for m in ("crossfit", "v2_full", "v41")}
    below_v41 = int(sum(per_set["S2"]["crossfit"][t] < per_set["S2"]["v41"].get(t, np.nan) for t in per_set["S2"]["crossfit"]))
    return {"verdict": "HELD_OUT_HOLDS" if nb["crossfit"] >= nb["v2_full"] else "HELD_OUT_WORSE", "n_types_below_both_C2L_Tangram": nb,
            "n_types_crossfit_below_v41": below_v41}


# ----------------------------------------------------------------------------- withheld genes: record keeping and guards
def region_leftovers(out: Path, region: int) -> list:
    """Everything a region job of this run dir writes (props, SACE layers, control, region files); non-empty = the region was started."""
    out = Path(out)
    cand = [out / "region_checks" / f"FAILED_region_{region}.json", out / "region_checks" / f"region_{region}.json"]
    for k in range(K_FOLDS):
        cand += [out / "props" / f"fold{k}" / f"Xenium_region_{region}", out / "sace" / f"fold{k}" / f"Xenium_region_{region}"]
    cand += [out / "control" / "props" / f"Xenium_region_{region}", out / "control" / "sace" / f"Xenium_region_{region}"]
    return [c for c in cand if c.exists()]


def move_region_aside(out: Path, region: int, stamp: str) -> Path:
    """The explicit rerun path (WH_RERUN=1): move a started region's files to out/_rerun_aside_<stamp>_region<R>/ (nothing deleted)."""
    out = Path(out)
    dest = out / f"_rerun_aside_{stamp}_region{region}"
    for f in region_leftovers(out, region):
        t = dest / f.relative_to(out)
        t.parent.mkdir(parents=True, exist_ok=True)
        f.rename(t)
    return dest


def regressions_problems(out, sha=None, repo=REPO, paths=vo.CODE_TREE) -> list:
    """Re-read the numbers of the five region files and check them (existence is not enough): exact names region_0..4, each with all 5
    folds, a full fit equal to X1 (<= 1e-9), every fold equal to the Gate 4 fold proportions (<= 1e-8), no fold fallback, the fold
    withheld counts summing to the 405-gene panel, the SACE control against X3 within CONTROL_TOL, one code-tree content throughout
    (and equal to ``sha``'s when given). Returns problem strings; empty = clean."""
    out, errs = Path(out), []
    d = out / "region_checks"
    names = sorted(f.name for f in d.glob("region_*.json")) if d.is_dir() else []
    want = [f"region_{r}.json" for r in range(N_REGIONS)]
    if names != want:
        return [f"region files are {names}, expected exactly {want}"]
    sigs = {}
    for r in range(N_REGIONS):
        rec = json.loads((d / f"region_{r}.json").read_text())
        tag = f"region {r}"
        if rec.get("region") != r:
            errs.append(f"{tag}: file records region {rec.get('region')!r}")
        x = rec.get("full_vs_x1_max_abs")
        if not (isinstance(x, (int, float)) and np.isfinite(x) and x <= FULL_TOL):
            errs.append(f"{tag}: full fit vs X1 = {x!r} (> {FULL_TOL})")
        folds = rec.get("folds") or {}
        if sorted(folds) != [str(k) for k in range(K_FOLDS)]:
            errs.append(f"{tag}: folds recorded are {sorted(folds)}")
        for k, f in folds.items():
            dd = f.get("vs_old_gate4_props_max_abs")
            if not (isinstance(dd, (int, float)) and np.isfinite(dd) and dd <= FOLD_TOL):
                errs.append(f"{tag} fold {k}: vs Gate 4 fold props = {dd!r} (> {FOLD_TOL})")
            if f.get("fallback") is not None:
                errs.append(f"{tag} fold {k}: fold fit fell back ({f.get('fallback')!r})")
        n_w = sum(int(f.get("n_withheld", 0)) for f in folds.values())
        if n_w != N_GENES:
            errs.append(f"{tag}: withheld genes sum to {n_w}, expected {N_GENES}")
        c = (rec.get("control_x3") or {}).get("max_abs_diff")
        if not (isinstance(c, (int, float)) and np.isfinite(c) and c <= CONTROL_TOL):
            errs.append(f"{tag}: SACE control vs X3 = {c!r} (> {CONTROL_TOL})")
        try:
            sigs[r] = vo.code_tree_signature(rec.get("sha") or "", repo, paths)
        except ValueError as e:
            errs.append(f"{tag}: no usable commit sha in the region file: {e}")
    if len(set(sigs.values())) > 1:
        errs.append("regions ran on different code (code-tree content differs between region files)")
    if sha:
        try:
            if sigs and set(sigs.values()) != {vo.code_tree_signature(sha, repo, paths)}:
                errs.append(f"region files' code differs from the scoring job's commit {sha[:12]}")
        except ValueError as e:
            errs.append(str(e))
    return errs


def block_counts(rows, methods=METHODS, set_name="S2") -> tuple:
    """(counts, problems). counts[method] = number of SCORED blocks of ``set_name``: rows with a finite rmse, one per
    (fold, region, cell type). Expected per method: 5 folds x 5 regions x 7 types = 175 (25 per type), fixed constants. Rows with a NaN
    rmse, duplicate (fold, region, type) blocks, extra or missing methods and a cell-type count other than 7 are problems."""
    df = pd.DataFrame(rows, columns=["set", "fold", "method", "region", "cell_type", "rmse"])
    df = df[df["set"] == set_name]
    counts, probs = {}, []
    extra = sorted(set(df["method"]) - set(methods))
    if extra:
        probs.append(f"unexpected methods scored: {extra}")
    for m in methods:
        g = df[(df["method"] == m) & np.isfinite(pd.to_numeric(df["rmse"], errors="coerce"))]
        counts[m] = int(len(g))
        if len(g) != EXPECTED_BLOCKS:
            probs.append(f"{m}: {len(g)} scored blocks, expected {EXPECTED_BLOCKS}")
        if len(g.drop_duplicates(["fold", "region", "cell_type"])) != len(g):
            probs.append(f"{m}: duplicate (fold, region, type) blocks")
        per = g.groupby("cell_type").size()
        if len(per) != N_CELL_TYPES or (per != BLOCKS_PER_TYPE).any():
            probs.append(f"{m}: blocks per cell type {per.to_dict()}, expected {N_CELL_TYPES} types x {BLOCKS_PER_TYPE}")
    return counts, probs


def rule_problems(sha: str, repo=REPO, rule=DECISION_RULE) -> dict:
    """Verify the decision rule against git, not only by recording its hash: the working file equals the blob committed at ``sha`` and has
    had exactly one blob across its whole history (never edited after it was committed). Returns {"blob", "problems", "sha256"}."""
    rel = str(Path(rule).resolve().relative_to(Path(repo).resolve()))

    def git(*a):
        return subprocess.run(["git", "--no-optional-locks", "-C", str(repo), *a], capture_output=True, text=True, check=False)

    errs = []
    committed = git("rev-parse", f"{sha}:{rel}")
    working = git("hash-object", str(rule))
    if committed.returncode != 0:
        errs.append(f"rule file is not committed at {sha[:12]}: {committed.stderr.strip()[:120]}")
    elif committed.stdout.strip() != working.stdout.strip():
        errs.append("rule file on disk differs from the committed blob")
    hist = git("log", "--format=%H", "--", rel)
    blobs = set()
    for c in (x for x in hist.stdout.split() if x):
        b = git("rev-parse", f"{c}:{rel}")
        if b.returncode == 0:
            blobs.add(b.stdout.strip())
    if len(blobs) != 1:
        errs.append(f"rule file has {len(blobs)} distinct committed blobs in its history (it must never change)")
    return {"blob": committed.stdout.strip(), "problems": errs, "sha256": hashlib.sha256(Path(rule).read_bytes()).hexdigest()}


def competitor_pin_problems(roots: dict, repo=REPO, pin_path=None) -> dict:
    """confirm_pinned on the reused competitor layers against v2_reused_inputs.json, and each competitor root must BE the pinned dir."""
    pin = json.loads(Path(pin_path or vo.REUSED_PIN).read_text())["groups"]
    groups = {"Cell2Location": "cell2location_layers", "Tangram": "tangram_native_layers"}
    sub = {g: pin[g] for g in groups.values()}
    res = vo.confirm_pinned(sub, repo)
    errs = []
    if not res["ok"]:
        errs.append(f"pinned competitor layers differ: changed {res['changed'][:3]} missing {res['missing'][:3]} added {res['added'][:3]}")
    for m, g in groups.items():
        if Path(roots[m]).resolve() != (Path(repo) / pin[g]["dir"]).resolve():
            errs.append(f"{m} root {roots[m]} is not the pinned dir {pin[g]['dir']}")
    return {"result": res, "problems": errs}


def x3_layer_problems(repo=REPO) -> list:
    """The X3 layers (v2_full) against the hashes of their own layers_manifest.json."""
    man = json.loads((vo.sace_root() / "layers_manifest.json").read_text())["layers"]
    if len(man) != N_REGIONS * N_CELL_TYPES:
        return [f"X3 manifest lists {len(man)} layers, expected {N_REGIONS * N_CELL_TYPES}"]
    return [f"X3 layer changed or missing: {rel}" for rel, h in man.items()
            if not (Path(repo) / rel).is_file() or hashlib.sha256((Path(repo) / rel).read_bytes()).hexdigest() != h]


def layer_control(a_dir, b_dir) -> dict:
    """Max abs difference over all 7 layer files of two flat sace dirs (Xenium_region_R/layers/*_layer.csv); a differing file set, spot
    index or gene columns is an error."""
    a_dir, b_dir = Path(a_dir), Path(b_dir)
    fa = sorted(p.name for p in (a_dir / "layers").glob("*_layer.csv"))
    fb = sorted(p.name for p in (b_dir / "layers").glob("*_layer.csv"))
    if fa != fb or len(fa) != N_CELL_TYPES:
        raise ValueError(f"layer files differ or are not {N_CELL_TYPES}: {fa} vs {fb}")
    worst = 0.0
    for n in fa:
        a, b = pd.read_csv(a_dir / "layers" / n, index_col=0), pd.read_csv(b_dir / "layers" / n, index_col=0)
        if a.shape != b.shape or list(a.columns) != list(b.columns) or list(a.index) != list(b.index):
            raise ValueError(f"layer {n}: shape/index/columns differ")
        if a.size:
            worst = max(worst, float(np.nanmax(np.abs(a.values - b.values))))
    return {"n_layers": len(fa), "max_abs_diff": worst}


# ----------------------------------------------------------------------------- withheld genes: jobs
def cmd_wh_region(args) -> None:
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    import time

    vo.assert_code_tree()
    sha = os.environ["SHA"]
    out, region = Path(args.out), int(args.region)
    # A finished region is never overwritten. A started one (props, SACE layers, control or a FAILED file) is refused unless this is the
    # explicit rerun path (WH_RERUN=1), which moves it aside (nothing deleted) and records that it was rerun.
    vo.refuse_existing(out / "region_checks" / f"region_{region}.json")
    rerun = os.environ.get("WH_RERUN") == "1"
    if region_leftovers(out, region):
        if not rerun:
            raise SystemExit(f"STOP: region {region} was already started in {out} ({[str(x.relative_to(out)) for x in region_leftovers(out, region)][:3]}); "
                             "set WH_RERUN=1 to move it aside and rerun it")
        print("rerun: moved aside to", move_region_aside(out, region, time.strftime("%Y%m%d_%H%M%S")), flush=True)
    rec0 = out / "region_checks" / "region_0.json"
    if region != 0:
        if not rec0.is_file():
            raise SystemExit("STOP: region 0 (whfirst) has not finished; later regions run at its commit")
        if vo.code_tree_signature(json.loads(rec0.read_text())["sha"]) != vo.code_tree_signature(sha):
            raise SystemExit("STOP: this job's commit holds different code than region 0's commit; later regions must run at the first stage's commit")
    types, data = data_io.load_xenium(7)
    genes = list(data[region]["genes"])
    for r in range(N_REGIONS):
        assert list(data[r]["genes"]) == genes, "gene panel differs across regions; folds would not be shared"
    d = data[region]
    p0 = rs.ms_share(d["X"], d["names"], rs.PROF[7])
    p_full, diag = refine_proportions(d["Y"], d["nc"], p0, n_bootstrap=0, **vo.FROZEN)
    if diag["fallback"]:
        raise SystemExit(f"STOP: full fit fell back for region {region}")
    x1 = pd.read_csv(vo.pred_root(7) / f"Xenium_region_{region}" / "predictions" / "baseline.csv.gz", index_col=0)
    full_diff = vo.max_abs_diff_frames(pd.DataFrame(p_full, index=d["obs"], columns=types), x1)
    print(f"region {region}: full fit vs X1 max abs diff {full_diff:.3e}", flush=True)
    if not full_diff <= FULL_TOL:
        raise SystemExit(f"STOP: full fit differs from the stamped X1 predictions by {full_diff}")
    sys.path.insert(0, str(XB / "spots_baseline"))
    import sace_props_ablation as sp

    # Control (F4): the fold SACE path, given the FULL-FIT proportions of this run, must reproduce the X3 layers. Without it the
    # crossfit-minus-v2_full difference could contain SACE run-to-run noise.
    cdir = out / "control" / "props" / f"Xenium_region_{region}" / "predictions"
    cdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(p_full, index=d["obs"], columns=types).to_csv(cdir / "baseline.csv.gz")
    sp.run_region("marker_share", region, out / "control" / "sace", str(out / "control" / "props"), "ms_v41", cell_mode="no_round")
    ctl = layer_control(out / "control" / "sace" / f"Xenium_region_{region}", vo.sace_root() / f"Xenium_region_{region}")
    ctl["tolerance"] = CONTROL_TOL
    print(f"region {region}: SACE control (full-fit props, fold path) vs X3 layers max abs diff {ctl['max_abs_diff']:.3e}", flush=True)
    if not ctl["max_abs_diff"] <= CONTROL_TOL:
        raise SystemExit(f"STOP: the fold SACE path does not reproduce the X3 layers for the full-fit proportions (max abs diff {ctl['max_abs_diff']})")
    checks = {"region": region, "sha": sha, "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "rerun": rerun, "n_genes": len(genes),
              "full_vs_x1_max_abs": full_diff, "control_x3": ctl, "folds": {}}
    for k, (fold, mask) in enumerate(gg.crossfit_masks(len(genes))):
        if mask[fold].any():
            raise SystemExit("STOP: a fold gene is still in the fit mask")
        p_k, dk = refine_proportions(d["Y"], d["nc"], p0, n_bootstrap=0, gene_mask=mask, **vo.FROZEN)
        old = pd.read_csv(OLD_FOLD_PROPS / f"fold{k}" / f"Xenium_region_{region}" / "predictions" / "baseline.csv.gz", index_col=0)
        dd = vo.max_abs_diff_frames(pd.DataFrame(p_k, index=d["obs"], columns=types), old)
        checks["folds"][str(k)] = {"n_withheld": int(len(fold)), "vs_old_gate4_props_max_abs": dd, "fallback": dk["fallback"]}
        print(f"region {region} fold {k}: vs Gate 4 fold props max abs diff {dd:.3e}", flush=True)
        if not dd <= FOLD_TOL:
            (out / "region_checks").mkdir(parents=True, exist_ok=True)
            (out / "region_checks" / f"FAILED_region_{region}.json").write_text(json.dumps(checks, indent=2))
            raise SystemExit(f"STOP: fold {k} proportions differ from the Gate 4 fold proportions by {dd}")
        pdir = out / "props" / f"fold{k}" / f"Xenium_region_{region}" / "predictions"
        pdir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(p_k, index=d["obs"], columns=types).to_csv(pdir / "baseline.csv.gz")
        sp.run_region("marker_share", region, out / "sace" / f"fold{k}", str(out / "props" / f"fold{k}"), "ms_v41", cell_mode="no_round")
    (out / "region_checks").mkdir(parents=True, exist_ok=True)
    with open(out / "region_checks" / f"region_{region}.json", "x") as fh:  # exclusive: never overwrites
        fh.write(json.dumps(checks, indent=2))
    print("withheld region done", region, flush=True)


def score_blocks(out: Path, roots: dict, folds_genes: list, cx, bp, sc_) -> list:
    """rows: one per (set, fold, method, region, type) block, fair standard (all GT spots, absent = 0), fold genes only."""
    gt_root = bp.gt_gex_dir("rnaprot_v41")
    rows = []
    for r in range(N_REGIONS):
        for stem in cx.CELL_TYPES:
            gt_df = sc_.norm_df(cx.load_csv_spots_genes(gt_root / f"Xenium_region_{r}" / f"{stem}_GT.csv"))
            full_preds = {m: sc_.load_pred(roots[m], r, stem, cx, gt_df.index) for m in roots if m != "crossfit"}
            for k, fg in enumerate(folds_genes):
                preds = {**full_preds, "crossfit": sc_.load_pred(out / "sace" / f"fold{k}", r, stem, cx, gt_df.index)}
                gt_f = restrict(gt_df, fg)
                preds_f = {m: (None if d is None else restrict(d, fg)) for m, d in preds.items()}
                for s in ("S1", "S2"):
                    info, blocks = sc_.build_blocks(gt_f, preds_f, s)
                    if info is None:
                        continue
                    for m, blk in blocks.items():
                        if blk is None:
                            continue
                        pv, gv, g = blk
                        met = cx.compute_metrics(pv, gv)
                        rows.append(dict(set=s, fold=k, method=m, region=r, cell_type=stem, rmse=met["rmse"], flattened_r=met["flattened_r"],
                                         n_genes=len(g), n_scored_spots=pv.shape[0]))
    return rows


def summarize(rows: list) -> dict:
    df = pd.DataFrame(rows)
    out = {}
    for s, g in df.groupby("set"):
        per_fold = g.groupby(["method", "cell_type", "fold"])[["rmse", "flattened_r"]].mean().reset_index()  # mean over regions
        per = per_fold.groupby(["method", "cell_type"])[["rmse", "flattened_r"]].mean()  # mean over folds
        out[s] = {"rmse": {m: per.loc[m]["rmse"].to_dict() for m in per.index.get_level_values(0).unique()},
                  "flattened_r": {m: per.loc[m]["flattened_r"].to_dict() for m in per.index.get_level_values(0).unique()}}
    return out


def heldout_inputs(out, roots: dict) -> list:
    """Every file the withheld-gene result is computed from (stamp-shaped: a flat list of existing paths, for ``vo.stamp_artifact(inputs=...)``):
    the five region files, X1 predictions, the Gate 4 fold props, this run's fold props and fold SACE layers, X3 layers, the v41 layers,
    the pinned competitor layers, the GT GEX files, the h5ad inputs and the decision rule. Raises if any is missing."""
    out = Path(out)
    ins = [out / "region_checks" / f"region_{r}.json" for r in range(N_REGIONS)]
    ins += [vo.pred_root(7) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz" for r in range(N_REGIONS)]
    for k in range(K_FOLDS):
        for r in range(N_REGIONS):
            ins += [OLD_FOLD_PROPS / f"fold{k}" / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz",
                    out / "props" / f"fold{k}" / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"]
            ins += sorted((out / "sace" / f"fold{k}" / f"Xenium_region_{r}" / "layers").glob("*_layer.csv"))
    for m in ("v2_full", "v41", *COMPETITORS):
        ins += sorted(Path(roots[m]).glob("Xenium_region_*/layers/*_layer.csv"))
    pin = json.loads(Path(vo.REUSED_PIN).read_text())["groups"]
    ins += [REPO / rel for rel in pin["gt_gex"]["files"]]
    ins += [vo.H5AD_DIR / f"Xenium_region_{r}_{kind}.h5ad" for r in range(N_REGIONS) for kind in ("CITE", "GEX")]
    ins.append(DECISION_RULE)
    missing = [str(f) for f in ins if not Path(f).is_file()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} input files missing, e.g. {missing[:3]}")
    return ins


def _finish(result_f: Path, result: dict, rows, out: Path, stamp=None) -> None:
    """Write the result (+ blocks csv), stamp it when ``stamp`` is given (a callable taking the list of files), then exit by verdict.
    A bad verdict is written AND stamped: it is the honest record. The exit code still fails the job; whether a non-HOLDS verdict may reach the
    manuscript is the lineage check's call (it refuses it at --stage complete)."""
    result_f.write_text(json.dumps(result, indent=2, sort_keys=True, default=float))
    if rows is not None:
        pd.DataFrame(rows).to_csv(out / "heldout_blocks.csv", index=False)
    if stamp is not None:
        stamp([result_f] + ([out / "heldout_blocks.csv"] if rows is not None else []))
    v = result["verdict"]["verdict"]
    print("verdict", result["verdict"], flush=True)
    # The numbers are written first; the exit code then tells a job-success glance or an afterok chain that this is NOT a finished result.
    if v == "HELD_OUT_WORSE":
        print(f"STOP: HELD_OUT_WORSE (exit {EXIT_WORSE}): report the numbers in {result_f}, do not tune, do not rephrase", file=sys.stderr, flush=True)
        raise SystemExit(EXIT_WORSE)
    if v != "HELD_OUT_HOLDS":
        print(f"STOP: NO_VERDICT (exit {EXIT_NO_VERDICT}): {result['verdict'].get('reason')}", file=sys.stderr, flush=True)
        raise SystemExit(EXIT_NO_VERDICT)


def cmd_wh_score(args) -> None:
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg
    from Benchmarking.xenium_benchmarking.rna_refinement import same_cells_rescore as sc_
    import v42_stamp as W

    started = W.utc_now()
    vo.assert_code_tree()
    sha = os.environ["SHA"]
    out = Path(args.out)
    result_f = out / "heldout_genes.json"
    vo.refuse_existing(result_f)
    vo.refuse_existing(out / "heldout_blocks.csv")
    bp, cx = gg._load_scorer()
    roots = {"v2_full": vo.sace_root(), "v41": CELL_COUNT_SACE / "v41__NO_ROUND", **{m: cx.ALL_METHODS[m] for m in COMPETITORS}}
    # Integrity of what is read (a problem here is STOP, not a verdict): the rule file against git, the pinned competitor layers, the X3 layers.
    rule = rule_problems(sha)
    pin = competitor_pin_problems(roots)
    integrity = rule["problems"] + pin["problems"] + x3_layer_problems()
    if integrity:
        raise SystemExit("STOP: input integrity: " + "; ".join(integrity[:6]))
    genes = list(data_io.load_xenium(7)[1][0]["genes"])
    folds_genes = fold_genes_upper(genes, gg.gene_folds(len(genes)))
    if len(genes) != N_GENES or sum(len(g) for g in folds_genes) != N_GENES or len(set().union(*folds_genes)) != N_GENES:
        raise SystemExit(f"STOP: gene folds are not a disjoint cover of the {N_GENES}-gene panel")
    reg_problems = regressions_problems(out, sha)
    base = {
        "rule": DECISION_RULE.name, "rule_sha256": rule["sha256"], "rule_git_blob": rule["blob"], "sha": sha,
        "roots": {m: str(p) for m, p in roots.items()}, "folds": {str(k): sorted(g) for k, g in enumerate(folds_genes)},
        "regression_problems": reg_problems, "competitor_pin": {k: pin["result"][k] for k in ("n_files", "ok")},
    }
    run_id = vo.run_id("X10")
    role = {"heldout_genes.json": "x10_heldout_genes", "heldout_blocks.csv": "x10_heldout_blocks_csv"}

    def stamper(ins):
        return lambda files: [vo.stamp_artifact(f, run_id=run_id, inputs=ins, producer=vo._rel(__file__), started_at=started, role=role[Path(f).name]) for f in files]

    if reg_problems:
        regs = [json.loads(f.read_text()) for f in sorted((out / "region_checks").glob("region_*.json"))]
        try:
            ins0 = heldout_inputs(out, roots)
        except FileNotFoundError as e:  # a failed region run may lack fold outputs; the diagnostic file is then written without a stamp
            print(f"NOTE: NO_VERDICT file not stamped, inputs incomplete: {e}", flush=True)
            ins0 = None
        _finish(result_f, {**base, "verdict": {"verdict": "NO_VERDICT", "reason": f"regression check failed: {reg_problems[:5]}"}, "regression_checks": regs},
                None, out, stamp=stamper(ins0) if ins0 else None)
    regs = [json.loads((out / "region_checks" / f"region_{r}.json").read_text()) for r in range(N_REGIONS)]
    rows = score_blocks(out, roots, folds_genes, cx, bp, sc_)
    counts, block_probs = block_counts(rows)
    summ = summarize(rows)
    v = verdict({s: summ[s]["rmse"] for s in summ}, counts, True, block_problems=block_probs)
    s2 = summ.get("S2", {}).get("rmse", {})
    delta = {t: s2["crossfit"][t] - s2["v2_full"][t] for t in s2.get("crossfit", {}) if t in s2.get("v2_full", {})}
    ins = heldout_inputs(out, roots)
    result = {
        **base, "verdict": v, "per_set": summ, "crossfit_minus_v2_full_rmse_S2": delta, "scored_blocks_S2": counts,
        "expected_blocks_per_method": EXPECTED_BLOCKS, "block_problems": block_probs, "regression_checks": regs,
        "input_files": {str(Path(f).relative_to(REPO)): hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in ins},
    }
    for t, x in delta.items():
        print(f"{t:16s} crossfit-v2_full rmse {x:+.4f}", flush=True)
    _finish(result_f, result, rows, out, stamp=stamper(ins))


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("x8").set_defaults(fn=cmd_x8)
    p = sub.add_parser("wh-region")
    p.add_argument("out")
    p.add_argument("region")
    p.set_defaults(fn=cmd_wh_region)
    p = sub.add_parser("wh-score")
    p.add_argument("out")
    p.set_defaults(fn=cmd_wh_score)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
