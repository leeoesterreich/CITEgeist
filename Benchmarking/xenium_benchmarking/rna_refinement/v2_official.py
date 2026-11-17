"""Official Xenium v2 arm (v42 plan Task 8): stages X1 (frozen v2 proportions), X3 manifest, X9 (reused inputs), and the stamp step.

The other stages reuse existing producers by option (see ``.commands/submit-v2-official.sh``):
  X2/B15  spots_baseline/score_p138_strict.py --arm v2          -> evaluation/results_rnaprot_v2/ (+ results_rnagt_v2/)
  X3      spots_baseline/sace_props_ablation.py ... no_round    -> CITEgeist/results_sace_rna_refined_v2_7type/<X3>/
  X4      evaluation/src/compare_xenium_gex_methods.py (XEN_REGISTER_V2_SACE=1, --spots all) -> evaluation/rnaprot_v2_fair/{own,matched}/
  X5      scripts/compute_gex_rmse_sd_ratio.py (XEN_V2_SACE_DIR)        -> evaluation/rnaprot_v2_fair/
  X6      CITEgeist/src/eval_single_cell_assignment.py                  -> CITEgeist/output_sc_assignment_v2/
Every output dir is NEW (``refuse_existing``); run ids come from ``Benchmarking/midkine_remediation/v42_run_ids.yaml``.

CLI (CITEgeist_env python, from the repo root, inside an sbatch job):
  v2_official.py x1                       compute frozen v2 (6t, 7t) -> output_rna_refined_v2_{6,7}type/<X1>/, stamp, check vs variants.json
  v2_official.py x3-manifest              layers manifest for the X3 SACE root + stamp
  v2_official.py x9 OUT_DIR               sha256 of every reuse-allowlist input of the lineage manifest vs disk (OUT_DIR must be new)
  v2_official.py x2-check | x4-check | x6-check   acceptance of a finished stage against the experiment record (exit 1 on a miss)
  v2_official.py paths                    shell assignments of the official dirs and run ids (eval'd by the submit script)
  v2_official.py check-code SHA           refuse (exit 1) if tracked code differs from SHA (the code tree the stages import)
  v2_official.py complete STAGE           write the stage-complete marker (a rerun of an unfinished stage is allowed, of a finished one is not)
  v2_official.py pin-reused               write v2_reused_inputs.json: sha256 of every reused input the stages read (X9 verifies it)
  v2_official.py stamp --stage X2 --producer REL --artifact PATH [--input PATH ...] [--input-glob DIR PATTERN ...] [--role R]
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
RUN_IDS = REPO / "Benchmarking" / "midkine_remediation" / "v42_run_ids.yaml"
LINEAGE_MANIFEST = HERE / "v2_reuse_allowlist.yaml"
VARIANTS_JSON = HERE / "runs" / "variants_20261031_1404_r8" / "variants.json"
N_REGIONS = 5
# Frozen v2 (plan Global Constraints): rna_refined, multinomial, c = 0.1, n_passes = 3, anchor_thr = 0. No retuning.
FROZEN = {"c": 0.1, "n_passes": 3, "anchor_thr": 0.0, "model": "multinomial"}
H5AD_DIR = REPO / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "h5ad_objects"
HEADLINE_TOL = 1e-6  # plan Step 1: X1 vs variants.json
# Everything the stage jobs import or run lives under these tracked paths; a job proves its code by requiring them to equal the commit SHA
# (a detached worktree cannot be used: the stage scripts resolve untracked data and imports from their own checkout path).
CODE_TREE = (
    "CITEgeist",
    "Benchmarking/xenium_benchmarking",
    "scripts/v42_stamp.py",
    "scripts/check_v42_lineage.py",
    "scripts/provenance.py",
    "scripts/compute_gex_rmse_sd_ratio.py",
    "scripts/compute_gex_rmse_sd_ratio_matched_c2l.py",
)
# X7 (protein-input stress test) also runs code under Benchmarking/protein_sensitivity (conditions, perturbation, strict scorer); its stages prove that tree too.
CODE_TREE_STRESS = CODE_TREE + ("Benchmarking/protein_sensitivity",)
MARKER_DIR = HERE / "runs" / "v2_official_complete"
REUSED_PIN = HERE / "v2_reused_inputs.json"
XBR = "Benchmarking/xenium_benchmarking"
# reused inputs the stages read (not produced by v2): group -> (dir relative to the repo, glob). X9 pins and verifies every file.
REUSED_GROUPS = {
    "cell2location_layers": (f"{XBR}/Cell2Location/output_singler_7type", "*_layer.csv"),
    "tangram_native_layers": (f"{XBR}/Tangram/output_singler_7type_native", "*_layer.csv"),
    "gt_gex": (f"{XBR}/ground_truth_xenium_rnaprot/ground_truth_7type/ground_truth_gex", "*_GT.csv"),
    "gt_props_7type": (f"{XBR}/ground_truth_xenium_rnaprot/ground_truth_7type/ground_truth", "*.csv"),
    "gt_props_6type": (f"{XBR}/ground_truth_xenium_rnaprot/ground_truth_6type/ground_truth", "*.csv"),
    "gt_cell_labels": (f"{XBR}/ground_truth_xenium_rnaprot/labels", "cell_labels_7type.csv"),
    "competitor_results_rnaprot": (f"{XBR}/evaluation/results_rnaprot/method_comparison_singler_6type", "*_results.json"),
    "competitor_results_rnagt": (f"{XBR}/evaluation/results_rnagt/method_comparison_singler_6type", "*_results.json"),
    "competitor_results_rnagt_v41p": (f"{XBR}/evaluation/results_rnagt_v41p/method_comparison_singler_6type", "*_results.json"),
    "sensitivity_scores_all": (f"{XBR}/reference_audit/singler_variants", "scores_all.csv"),
    "sensitivity_scores_v41p": (f"{XBR}/reference_audit/singler_variants_v41p", "scores_v41p.csv"),
    "sensitivity_variant_gt": (f"{XBR}/reference_audit/singler_variants/gt", "*.csv"),
    "sensitivity_canonical_gt": (f"{XBR}/ground_truth_singler/ground_truth_6type/ground_truth", "*.csv"),
    "h5ad_inputs": ("Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects", "*.h5ad"),
}


# ----------------------------------------------------------------------------- run ids and paths
def run_ids(path=RUN_IDS) -> dict:
    """{stage: run id} of the single protocol key of the run-id file (the same key holds B1-B12 and X1-X9)."""
    runs = (yaml.safe_load(Path(path).read_text()) or {}).get("runs") or {}
    if len(runs) != 1:
        raise ValueError(f"expected exactly one protocol key in {path}, found {sorted(runs)}")
    return dict(next(iter(runs.values())))


def run_id(stage: str, path=RUN_IDS) -> str:
    ids = run_ids(path)
    v = ids.get(stage)
    if not v or v == "NOT_RUN":
        raise ValueError(f"stage {stage} has no run id in {path} (value {v!r})")
    return str(v)


def pred_root(n_types: int) -> Path:
    """X1 output: Xenium_region_N/predictions/baseline.csv.gz below this dir (same layout the strict scorer reads for P138)."""
    return XB / "CITEgeist" / f"output_rna_refined_v2_{n_types}type" / run_id("X1")


def sace_root() -> Path:
    """X3 output: Xenium_region_N/{layers,props.csv,meta.json} below this dir (flat layout, as the marker_share arm writes)."""
    return XB / "CITEgeist" / "results_sace_rna_refined_v2_7type" / run_id("X3")


def results_dir(labelset: str = "rnaprot") -> Path:
    if labelset not in ("rnaprot", "rnagt"):
        raise ValueError(f"labelset must be rnaprot or rnagt, got {labelset!r}")
    return XB / "evaluation" / f"results_{labelset}_v2"


def fair_dir() -> Path:
    """X4/X5 output: own/ and matched/ (fair GEX scorer), plus the RMSE/sigma files."""
    return XB / "evaluation" / "rnaprot_v2_fair"


X5_DDOF1_SUFFIX = "ddof1"  # corrected X5 rerun: own dir, run id and stage-complete marker; the original x5 outputs are never touched


def x5_dir(variant=None) -> Path:
    """X5 output dir: the shared fair dir by default; variant 'ddof1' is a NEW dir beside it."""
    if variant is None:
        return fair_dir()
    if variant != X5_DDOF1_SUFFIX:
        raise ValueError(f"unknown X5 variant {variant!r}")
    return fair_dir().with_name(fair_dir().name + "_x5_" + variant)


def x5_run_id(variant=None, path=RUN_IDS) -> str:
    return run_id("X5", path) if variant is None else run_id("X51", path)  # X51 = first rerun of X5 (id <stage>1, digits only)


def sc_dir() -> Path:
    return XB / "CITEgeist" / "output_sc_assignment_v2"


def variants_scores_dir() -> Path:
    """B15: variant-GT rescore of the v2 6-type predictions."""
    return XB / "reference_audit" / "singler_variants_v2"


def refuse_existing(path) -> None:
    """Official outputs are written once; a rerun needs a new name (a dir that exists may be adopted/read-only)."""
    if Path(path).exists():
        raise FileExistsError(f"{path} exists; refuse (new dirs only, a rerun writes a NEW run id)")


def files_under(root, pattern: str) -> list:
    """Sorted files matching ``pattern`` below ``root``; a missing root or an empty match is an error, never an empty list."""
    root = Path(root)
    if not root.is_dir():
        raise FileNotFoundError(f"{root} is not a directory")
    found = sorted(p for p in root.rglob(pattern) if p.is_file())
    if not found:
        raise FileNotFoundError(f"no files match {pattern!r} under {root}")
    return found


def max_abs_diff_frames(a: pd.DataFrame, b: pd.DataFrame) -> float:
    """Max |a - b| over the same index and columns (order-free). A differing index or column set raises."""
    if set(a.index) != set(b.index) or len(a) != len(b):
        raise ValueError(f"row sets differ: {len(a)} vs {len(b)} rows")
    if set(a.columns) != set(b.columns):
        raise ValueError(f"column sets differ: {list(a.columns)} vs {list(b.columns)}")
    bb = b.loc[a.index, list(a.columns)]
    return float(np.abs(a.to_numpy(float) - bb.to_numpy(float)).max())


# ----------------------------------------------------------------------------- X1
def compute_frozen(n_types: int) -> dict:
    """Frozen v2 proportions per region from model inputs only (protein, raw spot counts, h5ad n_cells): {region: DataFrame}.

    Same calls as run_variants.py's frozen row of candidate a_V3_multinomial: p0 = marker share on the ms_v41 dictionary,
    then ``refine_proportions`` with the frozen settings."""
    sys.path.insert(0, str(REPO))
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    types, data = data_io.load_xenium(n_types)
    out = {}
    for r in range(N_REGIONS):
        d = data[r]
        p0 = rs.ms_share(d["X"], d["names"], rs.PROF[n_types])
        p, diag = refine_proportions(d["Y"], d["nc"], p0, n_bootstrap=0, **FROZEN)
        if diag["fallback"]:
            raise RuntimeError(f"refine_proportions fell back for region {r} ({n_types} types): not a frozen-v2 result")
        out[r] = pd.DataFrame(p, index=d["obs"], columns=types)
    return out


def write_predictions(preds: dict, root) -> list:
    """Write {region: DataFrame} as root/Xenium_region_N/predictions/baseline.csv.gz (never over an existing file)."""
    files = []
    for r, df in sorted(preds.items()):
        f = Path(root) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
        if f.exists():
            raise FileExistsError(f)
        f.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(f)
        files.append(f)
    return files


def x1_inputs(region: int) -> list:
    return [H5AD_DIR / f"Xenium_region_{region}_CITE.h5ad", H5AD_DIR / f"Xenium_region_{region}_GEX.h5ad"]


def headline_vs_variants(n_types: int, preds: dict, scratch) -> dict:
    """Score ``preds`` with the experiment's official scorer and compare with the variants.json frozen row (tolerance 1e-6)."""
    from Benchmarking.xenium_benchmarking.rna_refinement import score

    v = json.loads(VARIANTS_JSON.read_text())
    want = next(r["headline"] for r in v["rows"]
                if r["variant"] == "a_V3_multinomial" and r["row"] == "package default (frozen)" and r["n_types"] == n_types)
    got = score.score_official(n_types, preds, scratch)["headline"]
    return {"n_types": n_types, "got": got, "want": want, "abs_diff": abs(got - want), "ok": bool(abs(got - want) <= HEADLINE_TOL)}


def cmd_x1(_args) -> None:
    import tempfile
    import time

    sys.path.insert(0, str(REPO / "scripts"))
    import v42_stamp as W

    started = W.utc_now()
    roots = {n: pred_root(n) for n in (6, 7)}
    for p in roots.values():
        refuse_existing(p)
    assert_code_tree()
    W.preflight(__file__, [roots[n] / "Xenium_region_0" / "predictions" / "baseline.csv.gz" for n in (6, 7)], require_clean_tree=False)
    rid, checks = run_id("X1"), []
    t0 = time.time()
    all_preds = {}
    for n in (6, 7):
        preds = compute_frozen(n)
        all_preds[n] = preds
        chk = headline_vs_variants(n, preds, Path(tempfile.mkdtemp(prefix=f"x1_score_{n}_")) / "scored")
        checks.append(chk)
        print("X1", n, "types", chk, flush=True)
    if not all(c["ok"] for c in checks):
        raise SystemExit(f"STOP: X1 headline differs from variants.json beyond {HEADLINE_TOL}: {checks}")
    for n in (6, 7):
        files = write_predictions(all_preds[n], roots[n])
        for r, f in enumerate(files):
            stamp_artifact(f, run_id=rid, inputs=x1_inputs(r), producer=_rel(__file__), started_at=started, role="x1_frozen_v2_proportions")
        (roots[n] / "_COMPLETE.json").write_text(json.dumps(
            {"status": "COMPLETE", "stage": "X1", "run_id": rid, "n_types": n, "settings": FROZEN,
             "acceptance_vs_variants_json": next(c for c in checks if c["n_types"] == n), "seconds": round(time.time() - t0, 1)}, indent=2))
    print("X1 done", flush=True)


# ----------------------------------------------------------------------------- X3 manifest
def cmd_x3_manifest(_args) -> None:
    import hashlib

    sys.path.insert(0, str(REPO / "scripts"))
    import v42_stamp as W

    started = W.utc_now()
    root = sace_root()
    layers = files_under(root, "*_layer.csv")
    if len(layers) != N_REGIONS * 7:
        raise SystemExit(f"STOP: expected {N_REGIONS * 7} layer files under {root}, found {len(layers)}")
    out = root / "layers_manifest.json"
    refuse_existing(out)
    out.write_text(json.dumps(
        {"stage": "X3", "run_id": run_id("X3"), "cell_mode": "no_round", "props_dir": str(pred_root(7).relative_to(REPO)),
         "layers": {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest() for p in layers}}, indent=2, sort_keys=True))
    ins = [pred_root(7) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz" for r in range(N_REGIONS)] + layers
    stamp_artifact(out, run_id=run_id("X3"), inputs=ins, producer=_rel(__file__), started_at=started, role="x3_sace_no_round_layers")
    (root / "_COMPLETE.json").write_text(json.dumps({"status": "COMPLETE", "stage": "X3", "run_id": run_id("X3"), "n_layers": len(layers)}, indent=2))
    print("X3 manifest done", len(layers), flush=True)


# ----------------------------------------------------------------------------- code proof, stamping, completion markers
def code_tree_problems(sha: str, repo=REPO, paths=CODE_TREE) -> list:
    """Tracked non-markdown files of the code tree that differ from commit ``sha`` in the working tree or index, or between ``sha`` and HEAD."""
    import subprocess

    def names(*a):
        r = subprocess.run(["git", "--no-optional-locks", "-C", str(repo), "diff", "--name-only", *a, "--", *paths, ":(exclude)*.md"],
                           capture_output=True, text=True, check=False)
        if r.returncode != 0:
            raise RuntimeError(f"git diff failed: {r.stderr.strip()[:200]}")
        return [x for x in r.stdout.split("\n") if x]

    out = [f"differs from {sha[:12]} in the working tree/index: {n}" for n in names(sha)]
    out += [f"differs between {sha[:12]} and HEAD: {n}" for n in names(sha, "HEAD")]
    return out


def code_tree_signature(sha: str, repo=REPO, paths=CODE_TREE) -> str:
    return _code_tree_signature(sha, str(repo), tuple(paths))  # a commit's tree is immutable, so the git read is cached


@functools.lru_cache(maxsize=256)
def _code_tree_signature(sha: str, repo: str, paths: tuple) -> str:
    """sha256 over the blob ids and names of the tracked non-markdown files of the code tree AT commit ``sha``.

    Two commits with the same signature hold identical code, whatever else changed between them (docs, other threads' files).
    Raises if ``sha`` does not resolve or the tree is empty (an empty tree must never compare equal to anything real)."""
    import subprocess

    r = subprocess.run(["git", "--no-optional-locks", "-C", str(repo), "ls-tree", "-r", f"{sha}^{{commit}}", "--", *paths],
                       capture_output=True, text=True, check=False)
    if r.returncode != 0:
        raise ValueError(f"cannot read the code tree at {sha!r}: {r.stderr.strip()[:200]}")
    lines = sorted(x for x in r.stdout.split("\n") if x and not x.endswith(".md"))
    if not lines:
        raise ValueError(f"code tree at {sha!r} is empty for {paths}")
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def assert_code_tree(sha: str | None = None, repo=REPO, paths=CODE_TREE) -> None:
    sha = sha or os.environ.get("SHA")
    if not sha:
        raise SystemExit("STOP: SHA is not set; a stamped stage must name the commit that holds its code")
    problems = code_tree_problems(sha, repo, paths)
    if problems:
        raise SystemExit("STOP: code is not the committed code " + sha[:12] + ": " + "; ".join(problems[:8]))


def stamp_artifact(artifact, *, run_id, inputs, producer, started_at, role=None, root=REPO, check_code=True, denylist_path=None,
                   code_paths=CODE_TREE, reuse=None) -> dict:
    """v42 stamp of one artifact after proving the code tree equals the commit (SHA env); ``root`` = namespace and git root.

    ``reuse`` = {repo-relative path: sha256} of the named, hash-pinned reuse allowlist entries this stamp may take as inputs (default none)."""
    sys.path.insert(0, str(REPO / "scripts"))
    import v42_stamp as W

    if check_code:
        assert_code_tree(repo=root, paths=code_paths)
    kw = {} if denylist_path is None else {"denylist_path": denylist_path}
    if reuse:
        kw.update(reuse_paths=frozenset(reuse), reuse_hashes=dict(reuse))
    return W.write_stamp(Path(artifact), run_id=run_id, protocol_sha12="none", inputs=inputs, producer=producer, started_at=started_at,
                         repo_root=root, git_root=root, require_clean_tree=False, stage_role=role, **kw)


def cmd_check_code(args) -> None:
    problems = code_tree_problems(args.sha, paths=CODE_TREE_STRESS if getattr(args, "stress", False) else CODE_TREE)
    for p in problems:
        print("CODE:", p)
    if problems:
        raise SystemExit(1)
    print("code tree equals", args.sha)


def cmd_complete(args) -> None:
    MARKER_DIR.mkdir(parents=True, exist_ok=True)
    f = MARKER_DIR / f"{args.stage}.json"
    if f.exists():
        raise FileExistsError(f)
    f.write_text(json.dumps({"stage": args.stage, "sha": os.environ.get("SHA"), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                             "run_id": x5_run_id(X5_DDOF1_SUFFIX) if args.stage == "x5_ddof1" else run_ids().get({"x1": "X1", "x2": "X2", "x3": "X3", "x4": "X4", "x5": "X5", "x6": "X6", "b15": "X2", "x9": "X9"}[args.stage])}, indent=2))


# ----------------------------------------------------------------------------- X9
def reuse_entries(manifest_path=LINEAGE_MANIFEST) -> list:
    return list((yaml.safe_load(Path(manifest_path).read_text()) or {}).get("reuse_allowlist") or [])


def confirm_reused_inputs(entries, repo_root=REPO) -> dict:
    """sha256 of every reuse-allowlist file vs the hash the lineage manifest pins. An empty list is not a pass."""
    sys.path.insert(0, str(REPO / "scripts"))
    import check_v42_lineage as L

    changed, missing, ok_names = [], [], []
    for e in entries:
        f = Path(repo_root) / e["path"]
        if not f.is_file():
            missing.append({"name": e["name"], "path": e["path"]})
        elif L.sha256_file(f) != e["sha256"]:
            changed.append({"name": e["name"], "path": e["path"], "pinned": e["sha256"], "on_disk": L.sha256_file(f)})
        else:
            ok_names.append(e["name"])
    return {"n": len(entries), "ok": bool(entries) and not changed and not missing, "changed": changed, "missing": missing, "unchanged": ok_names}


def build_reused_pin(repo_root=REPO, groups=REUSED_GROUPS) -> dict:
    sys.path.insert(0, str(REPO / "scripts"))
    import check_v42_lineage as L

    out = {}
    for name, (d, pat) in groups.items():
        files = files_under(Path(repo_root) / d, pat)
        out[name] = {"dir": d, "glob": pat, "files": {str(f.relative_to(repo_root)): L.sha256_file(f) for f in files}}
    return out


def confirm_pinned(pin: dict, repo_root=REPO) -> dict:
    """Every pinned file exists with its pinned hash, and no file matching a group's glob is missing from the pin. Empty is not a pass."""
    sys.path.insert(0, str(REPO / "scripts"))
    import check_v42_lineage as L

    changed, missing, added = [], [], []
    n = 0
    for name, g in pin.items():
        listed = g["files"]
        for rel, want in listed.items():
            n += 1
            f = Path(repo_root) / rel
            if not f.is_file():
                missing.append(rel)
            elif L.sha256_file(f) != want:
                changed.append(rel)
        d = Path(repo_root) / g["dir"]
        now = {str(p.relative_to(repo_root)) for p in d.rglob(g["glob"]) if p.is_file()} if d.is_dir() else set()
        added += sorted(now - set(listed))
    return {"n_files": n, "ok": n > 0 and not (changed or missing or added), "changed": changed, "missing": missing, "added": added}


def cmd_pin_reused(args) -> None:
    out = Path(args.out) if args.out else REUSED_PIN
    refuse_existing(out)
    pin = build_reused_pin()
    out.write_text(json.dumps({"schema": 1, "groups": pin}, indent=1, sort_keys=True))
    print("pinned", {k: len(v["files"]) for k, v in pin.items()})


def cmd_x9(args) -> None:
    out = Path(args.out_dir)
    refuse_existing(out)
    res = confirm_reused_inputs(reuse_entries())
    pinned = confirm_pinned(json.loads(REUSED_PIN.read_text())["groups"])
    out.mkdir(parents=True)
    (out / "x9_reused_inputs.json").write_text(json.dumps({"manifest_allowlist": res, "pinned_stage_inputs": pinned}, indent=2))
    print("X9 allowlist", {k: res[k] for k in ("n", "ok")}, "pinned", {k: pinned[k] for k in ("n_files", "ok")},
          "changed", res["changed"] + pinned["changed"], "missing", res["missing"] + pinned["missing"], "added", pinned["added"], flush=True)
    if not (res["ok"] and pinned["ok"]):
        raise SystemExit("STOP: a reused input differs from its pinned hash, is missing, or a new file appeared in a pinned dir")


# ----------------------------------------------------------------------------- acceptance checks of finished stages
CELL_COUNT_JSON = HERE / "runs" / "cell_count_20261101_1424" / "cell_count.json"
SC_V2_JSON = HERE / "runs" / "sc_v2_20261101_1515" / "sc_v2.json"
X4_TOL = 1e-6
X6_TOL = 1e-3


def cmd_x2_check(_args) -> None:
    """X2 strict scoring vs the experiment: per_type_mean_r of the frozen row, 6t and 7t, within 1e-6 (a miss is stop-and-ask)."""
    v = json.loads(VARIANTS_JSON.read_text())
    bad = []
    for n in (6, 7):
        want = next(r["headline"] for r in v["rows"]
                    if r["variant"] == "a_V3_multinomial" and r["row"] == "package default (frozen)" and r["n_types"] == n)
        f = results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
        got = json.loads(f.read_text())["summary"]["per_type_mean_r"]
        print("X2", n, "types", got, want, abs(got - want), flush=True)
        if abs(got - want) > HEADLINE_TOL:
            bad.append((n, got, want))
    if bad:
        raise SystemExit(f"STOP: X2 differs from variants.json beyond {HEADLINE_TOL}: {bad}")


def cmd_x4_check(args) -> None:
    """X4 fair scoring vs cell_count.json (v2 NO_ROUND, competitors): per-type RMSE within 1e-6; own = S2, matched = S1."""
    from Benchmarking.xenium_benchmarking.rna_refinement import fair_gex_check as fc

    cc = json.loads(CELL_COUNT_JSON.read_text())
    out = Path(args.fair_dir)
    diffs = {mode: fc.compare(out / mode / "xenium_gex_comparison_7type.csv", cc[s], {"CITEgeist_v2_SACE": "v2__NO_ROUND"})
             for mode, s in (("own", "S2"), ("matched", "S1"))}
    worst = max(v for d in diffs.values() for m in d.values() for v in m.values())
    report = {"worst_abs_rmse_diff": worst, "tol": X4_TOL, "per_mode_method_type": diffs}
    (out / ("x4_regression.json" if worst <= X4_TOL else "_FAILED.json")).write_text(json.dumps(report, indent=2, default=float))
    print("X4 worst rmse diff", worst, flush=True)
    if not worst <= X4_TOL:
        raise SystemExit(f"STOP: X4 differs from cell_count.json beyond {X4_TOL}: worst {worst}")


def cmd_x6_check(args) -> None:
    got = json.loads(Path(args.eval_json).read_text())
    ref = json.loads(SC_V2_JSON.read_text())
    d1 = abs(got["current_hungarian"]["pooled_acc"] - ref["pooled_acc"]["v2"])
    d2 = abs(got["random"]["pooled_acc"] - ref["random"]["v2"]["pooled_acc"])
    print("X6 hungarian diff", d1, "random diff", d2, flush=True)
    if not (d1 <= X6_TOL and d2 <= X6_TOL):
        raise SystemExit(f"STOP: X6 differs from sc_v2.json beyond {X6_TOL}: hungarian {d1}, random {d2}")


def cmd_paths(args) -> None:
    ids = run_ids()
    items = {
        "PRED6": pred_root(6), "PRED7": pred_root(7), "SACE": sace_root(), "RES_RNAPROT": results_dir("rnaprot"),
        "RES_RNAGT": results_dir("rnagt"), "FAIR": fair_dir(), "SC": sc_dir(), "VSCORES": variants_scores_dir(),
    }
    for k, v in items.items():
        print(f"{k}={v}")
    for st in ("X1", "X2", "X3", "X4", "X5", "X6", "X9"):
        print(f"RID_{st}={ids[st]}")
    if getattr(args, "x5_ddof1", False):  # extra lines only on request: the default output is unchanged
        print(f"X5_OUT={x5_dir(X5_DDOF1_SUFFIX)}")
        print(f"RID_X5_DDOF1={x5_run_id(X5_DDOF1_SUFFIX)}")


# ----------------------------------------------------------------------------- generic stamp
def _rel(p) -> str:
    p = Path(p)
    p = p if p.is_absolute() else REPO / p
    return str(Path(os.path.normpath(p)).relative_to(os.path.normpath(REPO)))


def gather_stamp_inputs(root, inputs=(), input_glob=(), inputs_from_json=()) -> list:
    """Stamp inputs: explicit files, (dir, glob) pairs, and the repo-relative keys of FILE[KEY] maps (the producer's own record of what it read)."""
    root = Path(root)
    out = [Path(p) if Path(p).is_absolute() else root / p for p in inputs]
    for d, pat in input_glob or ():
        out += files_under(d if Path(d).is_absolute() else root / d, pat)
    for f, key in inputs_from_json or ():
        listed = json.loads(Path(f).read_text()).get(key)
        if not isinstance(listed, dict) or not listed:
            raise SystemExit(f"STOP: {f} has no non-empty {key!r} map to take stamp inputs from")
        out += [root / k for k in listed]
    return out


def cmd_stamp(args) -> None:
    sys.path.insert(0, str(REPO / "scripts"))
    import v42_stamp as W

    inputs = gather_stamp_inputs(REPO, args.input, args.input_glob, args.inputs_from_json)
    rid = args.run_id or run_id(args.stage)
    stamp_artifact(Path(args.artifact), run_id=rid, inputs=inputs, producer=args.producer,
                   started_at=args.started or W.utc_now(), role=args.role)
    print("stamped", args.artifact, "as", rid, flush=True)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("x1").set_defaults(fn=cmd_x1)
    sub.add_parser("x3-manifest").set_defaults(fn=cmd_x3_manifest)
    sub.add_parser("x2-check").set_defaults(fn=cmd_x2_check)
    pc = sub.add_parser("check-code")
    pc.add_argument("sha")
    pc.add_argument("--stress", action="store_true", help="also prove Benchmarking/protein_sensitivity (X7 stages)")
    pc.set_defaults(fn=cmd_check_code)
    pm = sub.add_parser("complete")
    pm.add_argument("stage", choices=["x1", "x2", "x3", "x4", "x5", "x5_ddof1", "x6", "b15", "x9"])
    pm.set_defaults(fn=cmd_complete)
    pr = sub.add_parser("pin-reused")
    pr.add_argument("--out", default=None)
    pr.set_defaults(fn=cmd_pin_reused)
    pp = sub.add_parser("paths")
    pp.add_argument("--x5-ddof1", action="store_true", help="also print X5_OUT / RID_X5_DDOF1 of the corrected X5 rerun")
    pp.set_defaults(fn=cmd_paths)
    p4 = sub.add_parser("x4-check")
    p4.add_argument("fair_dir")
    p4.set_defaults(fn=cmd_x4_check)
    p6 = sub.add_parser("x6-check")
    p6.add_argument("eval_json")
    p6.set_defaults(fn=cmd_x6_check)
    p9 = sub.add_parser("x9")
    p9.add_argument("out_dir")
    p9.set_defaults(fn=cmd_x9)
    ps = sub.add_parser("stamp")
    ps.add_argument("--stage", required=True)
    ps.add_argument("--producer", required=True, help="repo-relative tracked script that produced the artifact")
    ps.add_argument("--artifact", required=True)
    ps.add_argument("--input", action="append", default=[])
    ps.add_argument("--input-glob", nargs=2, action="append", metavar=("DIR", "PATTERN"))
    ps.add_argument("--inputs-from-json", nargs=2, action="append", metavar=("FILE", "KEY"),
                    help="take the repo-relative input paths from the keys of FILE[KEY] (the producer's own record of what it read)")
    ps.add_argument("--run-id", default=None, help="stamp this run id instead of the one of --stage (the X5 ddof1 rerun)")
    ps.add_argument("--role", default=None)
    ps.add_argument("--started", default=None)
    ps.set_defaults(fn=cmd_stamp)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
