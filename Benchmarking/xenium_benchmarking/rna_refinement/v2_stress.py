"""X7 (v42 plan Task 9): protein-input stress test of frozen v2 on the Xenium benchmark (6- and 7-type, strict rnaprot scoring).

The protein channels are perturbed exactly as in the v41p sweep (``protein_sensitivity.conditions.build`` + ``run_sensitivity.apply_condition``,
same seeds), then the two calls of the official X1 stage are made on the perturbed matrix:
  p0 = marker share (``reproduce_spike.ms_share``) and  p = ``refine_proportions(Y, n, p0, **FROZEN)``  (raw RNA and nuclei untouched).
Two arms come out of one pass: ``citegeist`` = frozen v2 (the CITEgeist estimate), ``spots`` = p0 (plain protein marker share, the
schema slot the manifest keys read: ``spots_vs_citegeist.json``). No SACE here: the sentinel GEX of the v41p chain used 5-cell
rounding, which the fair standard retires, and no X7 manifest key is a GEX key.

Commands (CITEgeist_env python, repo root, inside an sbatch job; SHA = the commit that holds the code, proved by v2_official.assert_code_tree):
  sweep TASK_ID        one unit = one (n_types, region) pair, 10 units (TASK_ID 0-9: 6-type regions 0-4, then 7-type regions 0-4)
  score                run check, baseline gate against X1/X2, strict scoring of both arms per seed, aggregation, v42 stamps
  paths                shell assignments (ROOT, RID_X7)
  units                the 10 units
Output root ``Benchmarking/protein_sensitivity/results/sensitivity_v2_rnaprot`` is NEW (the manifest names it). Unit dirs are
``citegeist/xenium_<n>_<r>`` (v2) and ``spots/spots_<n>_<r>`` (p0), the layout ``strict_score.unit_dir`` reads.
Acceptance: the baseline condition reproduces the stamped X1 predictions (<= 1e-9, checked in each unit before anything else is
computed) and the stamped X2 per_type_mean_r (<= 1e-6, at scoring); a miss is stop-and-ask, never tuned.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for _p in (str(REPO), str(REPO / "Benchmarking"), str(REPO / "Benchmarking" / "xenium_benchmarking"), str(REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from Benchmarking.xenium_benchmarking.rna_refinement import v2_official as vo  # noqa: E402

ROOT_REL = "Benchmarking/protein_sensitivity/results/sensitivity_v2_rnaprot"
ROOT = REPO / ROOT_REL
N_REGIONS = 5
UNITS = [(n, r) for n in (6, 7) for r in range(N_REGIONS)]
EXPECTED_N_CONDITIONS = {6: 50, 7: 52}  # fixed constants (as protein_sensitivity.run_check), never derived from the outputs
ARMS = ("citegeist", "spots")  # citegeist = frozen v2; spots = plain marker share p0 (schema slot of spots_vs_citegeist.json)
PRED_TOL = 1e-9
METRIC_TOL = 1e-6
FAIL_NO_MARKERS = "failed_no_markers"
FAIL_LT2 = "failed_fewer_than_2_types"
ALLOWED_FAILURES = {("ablation", "zero"), ("type_removal", "zero")}
VARIANTS_JSON = vo.VARIANTS_JSON
P0_REF_DEFAULT = VARIANTS_JSON.parent / "artifacts" / "a_V3_multinomial"  # the experiment's p0 artifact (marker share, unperturbed)
V41P_ROOT = REPO / "Benchmarking" / "protein_sensitivity" / "results" / "sensitivity_v41p_rnaprot"  # independent record of the v41p perturbation hashes


# ----------------------------------------------------------------------------- layout
def unit_dir(root, arm: str, n: int, r: int) -> Path:
    from protein_sensitivity import strict_score as SS

    return SS.unit_dir(root, "citegeist" if arm == "citegeist" else "spots", n, r)


def unit_of(task_id: int) -> tuple:
    if not 0 <= task_id < len(UNITS):
        raise SystemExit(f"STOP: task id {task_id} is not a unit (0-{len(UNITS) - 1})")
    return UNITS[task_id]


def conditions_for(n: int) -> list:
    """The sweep design of the v41p run (same builder, ms_v41 dictionary, brightness off); the fixed condition count is asserted."""
    from protein_sensitivity import conditions as C
    from protein_sensitivity import variant as V

    conds = json.loads(C.build("sweep", "xenium", n, V.from_args("ms_v41", 0)))
    if len(conds) != EXPECTED_N_CONDITIONS[n]:
        raise SystemExit(f"STOP: {n}-type sweep has {len(conds)} conditions, expected {EXPECTED_N_CONDITIONS[n]}")
    if conds[0]["axis"] != "baseline":
        raise SystemExit("STOP: the first condition must be the baseline (fail-fast acceptance relies on it)")
    return conds


def ensure_conditions(root, n: int, conds: list) -> None:
    """Write conditions_xenium_<n>.json once (exclusive create); a second writer must see identical text."""
    f = Path(root) / f"conditions_xenium_{n}.json"
    text = json.dumps(conds)
    try:
        with open(f, "x") as fh:
            fh.write(text)
    except FileExistsError:
        if f.read_text() != text:
            raise SystemExit(f"STOP: {f} exists with a different condition list")


def _dense(x):
    return np.asarray(x.toarray() if hasattr(x, "toarray") else x, float)


# ----------------------------------------------------------------------------- one condition (pure)
def estimate_arms(X, names, profile, Y, nc, frozen=None):
    """(p0, p, diag) for one perturbed protein matrix: the two calls of ``v2_official.compute_frozen``, with the condition's profile.

    Fewer than 2 types cannot be refined (``refine_proportions`` needs >= 2 supported types): the v2 arm then has no result
    (p is None) while p0 stays defined."""
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    p0 = rs.ms_share(X, names, profile)
    if len(profile) < 2:
        return p0, None, None
    p, diag = refine_proportions(Y, nc, p0, n_bootstrap=0, **(frozen or vo.FROZEN))
    return p0, p, diag


def _write_pred(df: pd.DataFrame, d: Path, cid: str) -> None:
    (d / "predictions").mkdir(parents=True, exist_ok=True)
    f = d / "predictions" / f"{cid}.csv.gz"
    if f.exists():
        raise FileExistsError(f)
    df.to_csv(f)


def run_unit(n: int, region: int, root, conds: list, x1_dir=None, p0_ref_dir=P0_REF_DEFAULT) -> dict:
    """All conditions of one unit into root/{citegeist,spots}/...; stops before the second condition if the baseline misses X1."""
    from protein_sensitivity import baseline as B
    from protein_sensitivity import run_sensitivity as RS
    from protein_sensitivity import variant as V

    root = Path(root)
    dirs = {a: unit_dir(root, a, n, region) for a in ARMS}
    for d in dirs.values():
        vo.refuse_existing(d)
    for d in dirs.values():
        d.mkdir(parents=True)
    variant = V.from_args("ms_v41", 0)
    profile_full = V.profile_for(n, variant)
    gex, cite = B.load_xenium(region)
    if not (cite.obs_names == gex.obs_names).all():
        raise SystemExit("STOP: CITE and GEX spot ids differ")
    Y = _dense(gex.X)
    nc = cite.obs["n_cells"].values.astype(float)
    obs = cite.obs_names
    metas = {a: [] for a in ARMS}
    unexpected = []  # exceptions that are NOT a recorded per-condition status (see below); they fail the unit
    for cond in conds:
        cid = RS.condition_id(cond)
        rng = np.random.default_rng(cond.get("seed", 0))
        c2, prof, meta = RS.apply_condition(cite.copy(), profile_full, cond, rng)
        status_v2 = status_p0 = meta["status"]
        sha = V.perturbed_sha256(c2) if meta["status"] == "ok" else None
        fallback, em = None, None
        if meta["status"] == "ok":
            try:
                types = list(prof)
                p0, p, diag = estimate_arms(_dense(c2.X), list(c2.var_names), prof, Y, nc)
                _write_pred(pd.DataFrame(p0, index=obs, columns=types), dirs["spots"], cid)
                if p is None:
                    status_v2 = FAIL_LT2
                else:
                    _write_pred(pd.DataFrame(p, index=obs, columns=types), dirs["citegeist"], cid)
                    fallback, em = diag["fallback"], list(map(bool, diag["em_converged"]))
            except Exception as e:
                # A ValueError from the fit (e.g. "Fewer than 2 supported types") is a recorded per-condition status, as in v41p.
                # ANY other exception (MemoryError, OSError, KeyError, ...) is a defect, not a result: it is recorded in the unit's
                # completion record and the unit exits nonzero (see cmd_sweep), instead of waiting for scoring to notice.
                msg = str(e).split(".")[0][:60].replace(" ", "_")
                status_v2 = status_p0 = "failed_" + (msg if isinstance(e, ValueError) else f"{type(e).__name__}_{msg}")
                if not isinstance(e, ValueError):
                    unexpected.append({"condition_id": cid, "exception": type(e).__name__, "message": str(e)[:300]})
        for a, st in (("citegeist", status_v2), ("spots", status_p0)):
            metas[a].append({
                "condition_id": cid, "status": st, "dropped": list(meta["dropped"]),
                "active_profile": {t: {"Major": list(sp.get("Major", []))} for t, sp in prof.items()},
                "perturbed_sha256": sha, "variant": V.meta(variant, n),
                **({"v2_fallback": fallback, "v2_em_converged": em} if a == "citegeist" else {}),
            })
        if cond["axis"] == "baseline":
            check_baseline_unit(n, region, dirs, metas, x1_dir, p0_ref_dir)
    for a in ARMS:
        (dirs[a] / "row_meta.jsonl").write_text("".join(json.dumps(m) + "\n" for m in metas[a]))
    cond_text = json.dumps(conds)
    run_meta = {
        "engine": "rna_refined", "frozen": vo.FROZEN, "n_types": n, "region": region, "n_conditions": len(conds),
        "conditions_sha256": hashlib.sha256(cond_text.encode()).hexdigest(), "sha": os.environ.get("SHA"),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "profile": V.meta(variant, n),
    }
    for a in ARMS:
        (dirs[a] / "run_meta.json").write_text(json.dumps({**run_meta, "arm": a}, indent=2))
    status_counts = {}
    for m in metas["citegeist"]:
        status_counts[m["status"]] = status_counts.get(m["status"], 0) + 1
    rec = {"unit": [n, region], "n_conditions": len(conds), "sha": run_meta["sha"], "unexpected_failures": unexpected, "status_counts": status_counts}
    # _UNIT_COMPLETE.json exists only for a unit with no unexpected exception; otherwise _UNIT_FAILED.json (the stage scripts and
    # check_run require the former), and cmd_sweep exits nonzero.
    (dirs["citegeist"] / ("_UNIT_FAILED.json" if unexpected else "_UNIT_COMPLETE.json")).write_text(json.dumps(rec, indent=2))
    return {**run_meta, "unexpected_failures": unexpected}


def check_baseline_unit(n, region, dirs, metas, x1_dir=None, p0_ref_dir=P0_REF_DEFAULT):
    """Fail-fast acceptance: baseline v2 equals the stamped X1 file; baseline p0 equals the experiment's p0 artifact. Raises SystemExit."""
    bm = {a: metas[a][-1] for a in ARMS}
    if bm["citegeist"]["status"] != "ok" or bm["citegeist"].get("v2_fallback"):
        raise SystemExit(f"STOP: baseline v2 did not run cleanly: {bm['citegeist']['status']} fallback={bm['citegeist'].get('v2_fallback')}")
    x1 = Path(x1_dir) if x1_dir else vo.pred_root(n)
    want = pd.read_csv(x1 / f"Xenium_region_{region}" / "predictions" / "baseline.csv.gz", index_col=0)
    got = pd.read_csv(dirs["citegeist"] / "predictions" / "baseline.csv.gz", index_col=0)
    d = vo.max_abs_diff_frames(got, want)
    print(f"unit {n}t r{region}: baseline v2 vs X1 max abs diff {d:.3e}", flush=True)
    if not d <= PRED_TOL:
        raise SystemExit(f"STOP: baseline v2 differs from the stamped X1 predictions by {d} (> {PRED_TOL}); report, do not tune")
    p0_ref = None if p0_ref_dir is None else Path(p0_ref_dir) / f"{n}type" / f"region_{region}" / "p0.csv"
    if p0_ref is not None:
        if not p0_ref.is_file():
            raise SystemExit(f"STOP: the experiment's p0 artifact {p0_ref} is missing")
        d0 = vo.max_abs_diff_frames(pd.read_csv(dirs["spots"] / "predictions" / "baseline.csv.gz", index_col=0), pd.read_csv(p0_ref, index_col=0))
        print(f"unit {n}t r{region}: baseline p0 vs experiment p0 artifact max abs diff {d0:.3e}", flush=True)
        if not d0 <= PRED_TOL:
            raise SystemExit(f"STOP: baseline p0 differs from the experiment's p0 artifact by {d0}")


# ----------------------------------------------------------------------------- run check (pure on files)
def check_against_v41p(root, n: int, region: int, ref_root) -> list:
    """Perturbation hashes of one unit vs the v41p run's ``row_meta.jsonl`` (an independent record written by another run, not by this
    sweep): every row of both arms must carry the hash v41p recorded for the same condition id, and the plain arm the same status.
    A missing reference unit or an unknown condition id is an error (never a silent pass)."""
    from protein_sensitivity import strict_score as SS

    ref_dir = SS.unit_dir(ref_root, "citegeist", n, region)
    rp = ref_dir / "row_meta.jsonl"
    if not rp.is_file():
        return [f"{n}-type region {region}: v41p reference {rp} is missing"]
    ref = {}
    for x in (json.loads(l) for l in rp.read_text().splitlines() if l.strip()):
        ref[x["condition_id"]] = (x.get("status"), x.get("perturbed_sha256"))
    if not ref:
        return [f"{n}-type region {region}: v41p reference {rp} has no rows"]
    errs = []
    for a in ARMS:
        rows = [json.loads(l) for l in (unit_dir(root, a, n, region) / "row_meta.jsonl").read_text().splitlines() if l.strip()]
        for x in rows:
            cid = x["condition_id"]
            if cid not in ref:
                errs.append(f"{n}-type region {region} {a}: condition {cid} is not in the v41p reference")
                continue
            if x.get("perturbed_sha256") != ref[cid][1]:
                errs.append(f"{n}-type region {region} {a}: {cid} perturbed_sha256 differs from the v41p record")
            if a == "spots" and x["status"] != ref[cid][0]:
                errs.append(f"{n}-type region {region} {a}: {cid} status {x['status']!r} differs from the v41p record {ref[cid][0]!r}")
    return errs


def check_run(root, n_types_list=(6, 7), expected=EXPECTED_N_CONDITIONS, ref_root=V41P_ROOT, score_sha=None, require_sha=False,
              repo=REPO, tree_paths=vo.CODE_TREE_STRESS) -> list:
    """Completeness of a finished sweep root; returns error strings (empty = clean). Counts are fixed constants.

    ref_root: the v41p run whose recorded perturbation hashes every unit must match (None only in synthetic tests).
    Commits: units are compared by CODE-TREE CONTENT (blob ids of the tracked non-markdown files of ``tree_paths`` at each unit's
    SHA), not by raw commit id, so a commit outside the code tree between stages is harmless and one inside it is an error. The
    scoring job's own SHA (``score_sha``) is compared the same way, so sweep code and score code can never differ unnoticed.
    require_sha: a unit without a recorded SHA is an error (the real stages always record it)."""
    from protein_sensitivity import run_sensitivity as RS
    from protein_sensitivity import strict_score as SS

    root, errs, shas = Path(root), [], set()
    for n in n_types_list:
        cf = root / f"conditions_xenium_{n}.json"
        if not cf.is_file():
            errs.append(f"missing {cf}")
            continue
        conds = json.loads(cf.read_text())
        ids = [RS.condition_id(c) for c in conds]
        if len(ids) != expected[n]:
            errs.append(f"{n}-type conditions file has {len(ids)} conditions, expected {expected[n]}")
        if len(set(ids)) != len(ids):
            errs.append(f"{n}-type conditions file has duplicate condition ids")
        want_sha = hashlib.sha256(cf.read_text().encode()).hexdigest()
        for r in range(N_REGIONS):
            metas = {}
            for a in ARMS:
                d = unit_dir(root, a, n, r)
                tag = f"{a}/{d.name}"
                p = d / "row_meta.jsonl"
                if not p.is_file() or not (d / "run_meta.json").is_file():
                    errs.append(f"{tag}: missing row_meta.jsonl or run_meta.json")
                    continue
                rows = [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
                metas[a] = rows
                if [x["condition_id"] for x in rows] != ids:
                    errs.append(f"{tag}: row_meta condition ids differ from the conditions file (missing, duplicate or reordered)")
                rm = json.loads((d / "run_meta.json").read_text())
                if rm.get("engine") != "rna_refined" or rm.get("frozen") != vo.FROZEN:
                    errs.append(f"{tag}: run_meta engine/frozen settings differ from frozen v2")
                if rm.get("conditions_sha256") != want_sha:
                    errs.append(f"{tag}: run_meta conditions_sha256 differs from the conditions file")
                shas.add(rm.get("sha"))
                if require_sha and not rm.get("sha"):
                    errs.append(f"{tag}: run_meta has no sha")
                uc = unit_dir(root, "citegeist", n, r) / "_UNIT_COMPLETE.json"
                if a == "citegeist":
                    if not uc.is_file():
                        errs.append(f"{tag}: _UNIT_COMPLETE.json is missing (the unit did not finish cleanly)")
                    elif json.loads(uc.read_text()).get("unexpected_failures") != []:
                        errs.append(f"{tag}: the unit recorded unexpected failures: {json.loads(uc.read_text()).get('unexpected_failures')}")
                for x in rows:
                    cid, st = x["condition_id"], x["status"]
                    cond = conds[ids.index(cid)] if cid in ids else {"axis": "?"}
                    if st == "ok":
                        if not (d / "predictions" / f"{cid}.csv.gz").is_file():
                            errs.append(f"{tag}: {cid} is ok but its prediction file is missing")
                        if not x.get("perturbed_sha256"):
                            errs.append(f"{tag}: {cid} is ok with a null perturbed_sha256")
                    elif not (
                        (st == FAIL_NO_MARKERS and (cond["axis"], cond.get("panel")) in ALLOWED_FAILURES)
                        or (st == FAIL_LT2 and a == "citegeist" and len(cond.get("profile") or {}) < 2)
                    ):
                        errs.append(f"{tag}: {cid} has unexpected status {st!r}")
                if a == "citegeist":
                    for x in rows:
                        if x["status"] == "ok" and x["condition_id"] == "baseline" and x.get("v2_fallback"):
                            errs.append(f"{tag}: baseline v2 fell back")
            if set(metas) == set(ARMS):
                try:
                    bad = SS.check_shared_perturbation(unit_dir(root, "citegeist", n, r), unit_dir(root, "spots", n, r))
                except (ValueError, OSError) as e:
                    errs.append(f"{n}-type region {r}: shared perturbation not checkable: {e}")
                else:
                    errs += [f"{n}-type region {r}: perturbation differs: {b}" for b in bad]
                if ref_root is not None:
                    errs += check_against_v41p(root, n, r, ref_root)
    shas.discard(None)
    sigs = {}
    for s_ in sorted(shas | ({score_sha} if score_sha else set())):
        try:
            sigs[s_] = vo.code_tree_signature(s_, repo, tree_paths)
        except ValueError as e:
            errs.append(f"commit {s_}: {e}")
    if len(set(sigs.values())) > 1:
        errs.append(f"units / the scoring job ran on different code (code-tree content differs between commits {sorted(sigs)})")
    return errs


# ----------------------------------------------------------------------------- gate (baseline acceptance at scoring)
def x2_reference(n: int) -> dict:
    f = vo.results_dir("rnaprot") / f"method_comparison_singler_{n}type" / f"CITEgeist_{n}type_results.json"
    d = json.loads(f.read_text())
    return {"path": f, "per_type_mean_r": float(d["summary"]["per_type_mean_r"]), "mask_report": d["mask_report"]}


def mask_matches_x2(ctx, ref_mask_report: dict) -> list:
    """Per-region mask counts of the scoring context vs the mask_report of the official X2 result."""
    errs = []
    for r in range(N_REGIONS):
        want = ref_mask_report[str(r)]
        got = len(ctx.masks[r])
        if got != int(want["n_mask_spots"]):
            errs.append(f"region {r}: mask has {got} spots, X2 mask_report says {want['n_mask_spots']}")
    return errs


def baseline_gate(ctx, root, n: int) -> dict:
    """v2 baseline vs stamped X2 per_type_mean_r (<= 1e-6) and X1 predictions (<= 1e-9); plain-share baseline vs the experiment's
    v0 headline (<= 1e-6). ok is False on any miss."""
    from protein_sensitivity import strict_score as SS

    ref = x2_reference(n)
    out = {"n_types": n, "x2_file": str(ref["path"].relative_to(REPO))}
    mask_errs = mask_matches_x2(ctx, ref["mask_report"])
    sc = {a: SS.score_condition(ctx, {r: unit_dir(root, a, n, r) for r in range(N_REGIONS)}, "baseline") for a in ARMS}
    v2_r, p0_r = float(sc["citegeist"]["summary"]["per_type_mean_r"]), float(sc["spots"]["summary"]["per_type_mean_r"])
    p0_want = float(json.loads(VARIANTS_JSON.read_text())["v0_marker_share"][str(n)])
    pred = {}
    for r in range(N_REGIONS):
        got = pd.read_csv(unit_dir(root, "citegeist", n, r) / "predictions" / "baseline.csv.gz", index_col=0)
        want = pd.read_csv(vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz", index_col=0)
        pred[str(r)] = vo.max_abs_diff_frames(got, want)
    out.update(
        mask_errors=mask_errs, v2_per_type_mean_r=v2_r, x2_per_type_mean_r=ref["per_type_mean_r"], v2_abs_diff=abs(v2_r - ref["per_type_mean_r"]),
        plain_per_type_mean_r=p0_r, plain_want=p0_want, plain_abs_diff=abs(p0_r - p0_want), pred_max_abs_diff_vs_x1=pred,
        tolerances={"metric": METRIC_TOL, "pred": PRED_TOL},
    )
    out["ok"] = bool(not mask_errs and out["v2_abs_diff"] <= METRIC_TOL and out["plain_abs_diff"] <= METRIC_TOL and max(pred.values()) <= PRED_TOL)
    return out


# ----------------------------------------------------------------------------- aggregation
def v2_allowed_failure(status, cond) -> bool:
    return (status == FAIL_NO_MARKERS and (cond["axis"], cond.get("panel")) in ALLOWED_FAILURES) or (
        status == FAIL_LT2 and len(cond.get("profile") or {}) < 2
    )


def comparison_rows(seed_level: pd.DataFrame, root, n_types_list=(6, 7)) -> dict:
    """spots_vs_citegeist.json rows (v41p schema; ``citegeist_*`` = frozen v2, ``spots_*`` = plain marker share) plus the v2 SEM and the
    number of region fits that fell back to a single pass."""
    import math

    from protein_sensitivity import aggregate_v41p as AGG

    out = AGG.build_spots_comparison(seed_level)
    fb = {}
    for n in n_types_list:
        for r in range(N_REGIONS):
            for x in (json.loads(l) for l in (unit_dir(root, "citegeist", n, r) / "row_meta.jsonl").read_text().splitlines() if l.strip()):
                if x["status"] == "ok":
                    fb.setdefault((n, x["condition_id"]), []).append(bool(x.get("v2_fallback")))
    cid_key = {}
    for n in n_types_list:
        sub = seed_level[(seed_level["method"] == "citegeist") & (seed_level["n_types"] == n)]
        for _, row in sub.iterrows():
            cid_key[(n, row["condition_id"])] = (row["axis"], AGG._cell_key(row["axis"], row["level"], row["panel"]))
    for row in out["rows"]:
        g = AGG._spots_group(seed_level, "citegeist", row["gt"]).get((row["condition"], row["level"]), {"ok": {}})
        k = len(g["ok"])
        row["citegeist_r_sem"] = (AGG._sd(g["ok"].values()) / math.sqrt(k)) if k > 1 else None
        flags = [f for (n, cid), v in fb.items() if cid_key.get((n, cid)) == (row["condition"], row["level"]) and n == row["gt"] for f in v]
        row["n_v2_region_fits"] = len(flags)
        row["n_v2_fallback_region_fits"] = int(sum(flags))
    out["description"] = (
        "Frozen v2 (citegeist_*) vs plain protein marker share (spots_*, the p0 of the same run) under identical protein perturbations "
        "(Xenium RCC, strict per-type mean r, rnaprot GT, one shared mask, ms_v41 dictionary). delta = v2 - plain. Units are (condition, "
        "seed) joint 5-region scores; matched = both arms scored ok for that seed. citegeist_r_sem = SD / sqrt(n seed units), null for "
        "deterministic conditions. n_v2_fallback_region_fits counts region fits where v2 returned its single-pass fallback."
    )
    return out


def score_root(root, n_list=(6, 7), expected=None, ref_root=V41P_ROOT, score_sha=None, require_sha=False) -> dict:
    """Run check, baseline gate (STOP on a miss), strict scoring of both arms per seed, aggregation; writes gate/ mask/ scores/ and the
    three JSON files into ``root`` (all must be new). No code proof, no stamps (cmd_score adds both)."""
    from protein_sensitivity import aggregate_v41p as AGG
    from protein_sensitivity import strict_score as SS

    root = Path(root)
    for o in (root / "gate", root / "mask", root / "scores", root / "sensitivity_summary.json", root / "spots_vs_citegeist.json", root / "inputs_record.json"):
        vo.refuse_existing(o)
    errs = check_run(root, n_list, expected or EXPECTED_N_CONDITIONS, ref_root=ref_root, score_sha=score_sha, require_sha=require_sha)
    if errs:
        raise SystemExit("STOP: run check is not clean:\n" + "\n".join(errs[:30]))
    (root / "gate").mkdir()
    (root / "mask").mkdir()
    ctxs, gates = {}, {}
    for n in n_list:
        ctx = SS.reconstruct_mask(n)
        ctxs[n] = ctx
        (root / "mask" / f"mask_{n}type.json").write_text(json.dumps(SS.mask_record(ctx)))
        gates[n] = baseline_gate(ctx, root, n)
        (root / "gate" / f"baseline_gate_{n}type.json").write_text(json.dumps(gates[n], indent=2, default=float))
        print("gate", n, {k: gates[n][k] for k in ("ok", "v2_abs_diff", "plain_abs_diff")}, flush=True)
    if not all(g["ok"] for g in gates.values()):
        raise SystemExit("STOP: baseline gate missed (see gate/); report the numbers, do not tune")
    (root / "scores").mkdir()
    seed_all, rt_all = [], []
    for n in n_list:
        conds = json.loads((root / f"conditions_xenium_{n}.json").read_text())
        for method in ARMS:
            s, rt = AGG.score_all(root, method, n, ctxs[n], conds, gate_fn=lambda *_a: None, allowed_failure=v2_allowed_failure)
            seed_all.append(s)
            rt_all.append(rt)
    seed_level = pd.concat(seed_all, ignore_index=True)
    for (m, n), g in seed_level.groupby(["method", "n_types"]):
        if not (g["status"] == "ok").any():
            raise SystemExit(f"STOP: {m} {n}-type: no condition was scored")
    region_type = pd.concat([x for x in rt_all if len(x)], ignore_index=True)
    means = AGG.condition_means(seed_level, region_type)
    seed_level.to_csv(root / "scores" / "seed_level.csv", index=False)
    region_type.to_csv(root / "scores" / "region_type.csv", index=False)
    flat = means.copy()
    for c in ("region_type_mean", "per_type"):
        flat[c] = flat[c].map(lambda x: json.dumps(x, sort_keys=True))
    flat.to_csv(root / "scores" / "strict_scores.csv", index=False)
    summary = {"_metric": AGG.METRIC_LABEL, "_arm": "frozen v2 (rna_refined), perturbed protein, raw RNA and nuclei unchanged"}
    for n in n_list:
        summary.update(AGG.summarize_strict(means, n))
    (root / "sensitivity_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True, default=float))
    (root / "spots_vs_citegeist.json").write_text(json.dumps(comparison_rows(seed_level, root, n_list), indent=2, sort_keys=True, default=float))
    rec = {"gate": {str(n): {k: g[k] for k in ("v2_per_type_mean_r", "x2_per_type_mean_r", "v2_abs_diff", "plain_per_type_mean_r", "plain_want")} for n, g in gates.items()},
           "scored_predictions": {f"{r['method']}/{r['n_types']}/{r['condition_id']}": json.loads(r["pred_sha256"]) for _, r in seed_level[seed_level["status"] == "ok"].iterrows()},
           "conditions": {str(n): vo_sha(root / f"conditions_xenium_{n}.json") for n in n_list}}
    (root / "inputs_record.json").write_text(json.dumps(rec, indent=2, sort_keys=True))
    return {"gates": gates, "seed_level": seed_level, "means": means, "summary": summary}


def stamp_inputs(root, n_list=(6, 7)) -> list:
    """Everything the two X7 artifacts were computed from: X1 predictions (the v2 upstream), the X2 results, the conditions, the input record, the h5ad inputs."""
    root = Path(root)
    ins = [vo.pred_root(n) / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz" for n in n_list for r in range(N_REGIONS)]
    ins += [x2_reference(n)["path"] for n in n_list]
    ins += [root / f"conditions_xenium_{n}.json" for n in n_list] + [root / "inputs_record.json"]
    ins += [vo.H5AD_DIR / f"Xenium_region_{r}_{k}.h5ad" for r in range(N_REGIONS) for k in ("CITE", "GEX")]
    return ins


def design_reuse(root=ROOT_REL, n_list=(6, 7), manifest_path=None) -> dict:
    """{path: sha256} of the conditions files X7 reuses from the v41p design: read from the lineage manifest's reuse allowlist
    (category ``design``), the one place the pins live. Refuses a missing entry, so a stamp can never name a reuse the manifest does not."""
    sys.path.insert(0, str(REPO / "scripts"))
    import check_v42_lineage as L

    allow = {a["path"]: a for a in L.load_manifest(manifest_path or L.DEFAULT_MANIFEST).get("reuse_allowlist") or []}
    out = {}
    for n in n_list:
        rel = f"{root}/conditions_xenium_{n}.json"
        a = allow.get(rel)
        if a is None or a.get("category") != "design":
            raise SystemExit(f"STOP: {rel} is not a design entry of the reuse allowlist in the lineage manifest")
        out[rel] = a["sha256"]
    return out


def cmd_score(_args) -> None:
    sys.path.insert(0, str(REPO / "scripts"))
    import v42_stamp as W

    started = W.utc_now()
    vo.assert_code_tree(paths=vo.CODE_TREE_STRESS)
    root = ROOT
    score_root(root, score_sha=os.environ.get("SHA"), require_sha=True)
    ins = stamp_inputs(root)
    for art, role in ((root / "spots_vs_citegeist.json", "x7_protein_stress_v2_vs_plain"), (root / "sensitivity_summary.json", "x7_protein_stress_v2_summary")):
        vo.stamp_artifact(art, run_id=vo.run_id("X7"), inputs=ins, producer=vo._rel(__file__), started_at=started, role=role, code_paths=vo.CODE_TREE_STRESS, reuse=design_reuse())
    complete("x7score")
    print("X7 scoring done", flush=True)


def vo_sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def complete(stage: str) -> None:
    vo.MARKER_DIR.mkdir(parents=True, exist_ok=True)
    f = vo.MARKER_DIR / f"{stage}.json"
    if f.exists():
        raise FileExistsError(f)
    f.write_text(json.dumps({"stage": stage, "sha": os.environ.get("SHA"), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                             "run_id": vo.run_ids().get({"x7": "X7", "x7score": "X7", "x8": "X8"}.get(stage, ""))}, indent=2))


# ----------------------------------------------------------------------------- CLI
def cmd_sweep(args) -> None:
    n, r = unit_of(int(args.task_id))
    vo.assert_code_tree(paths=vo.CODE_TREE_STRESS)
    ROOT.mkdir(parents=True, exist_ok=True)
    conds = conditions_for(n)
    ensure_conditions(ROOT, n, conds)
    print(f"X7 unit {n}-type region {r}: {len(conds)} conditions -> {ROOT_REL}", flush=True)
    meta = run_unit(n, r, ROOT, conds)
    if meta["unexpected_failures"]:
        raise SystemExit(f"STOP: X7 unit {n}-type region {r} finished with unexpected exceptions (see _UNIT_FAILED.json): "
                         f"{meta['unexpected_failures'][:3]}")
    print("X7 unit done", n, r, flush=True)


def first_stage_sha(kind: str, wh_out=None, root=ROOT) -> str:
    """The commit a multi-stage analysis was started at, read from the FIRST stage's own record (never HEAD): x7 = unit 0's run_meta,
    wh = region 0's region_checks file. Later stages (x7rest, x7score, whrest, whscore) run at this commit."""
    import re

    if kind == "x7":
        f = Path(root) / "citegeist" / "xenium_6_0" / "run_meta.json"
    elif kind == "wh":
        f = Path(wh_out) / "region_checks" / "region_0.json"
    else:
        raise ValueError(kind)
    if not f.is_file():
        raise SystemExit(f"STOP: first-stage record {f} does not exist")
    sha = json.loads(f.read_text()).get("sha")
    if not (isinstance(sha, str) and re.fullmatch(r"[0-9a-f]{40}", sha)):
        raise SystemExit(f"STOP: first-stage record {f} has no full commit sha ({sha!r})")
    vo.code_tree_signature(sha)  # raises if it does not resolve
    return sha


def cmd_first_sha(args) -> None:
    print(first_stage_sha(args.kind, args.wh_out))


def cmd_paths(_args) -> None:
    print(f"X7_ROOT={ROOT}")
    print(f"RID_X7={vo.run_id('X7')}")
    print(f"RID_X10={vo.run_id('X10')}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    ps = sub.add_parser("sweep")
    ps.add_argument("task_id")
    ps.set_defaults(fn=cmd_sweep)
    sub.add_parser("score").set_defaults(fn=cmd_score)
    sub.add_parser("paths").set_defaults(fn=cmd_paths)
    pf = sub.add_parser("first-sha")
    pf.add_argument("kind", choices=["x7", "wh"])
    pf.add_argument("wh_out", nargs="?")
    pf.set_defaults(fn=cmd_first_sha)
    pc = sub.add_parser("complete")
    pc.add_argument("stage", choices=["x7", "x7score", "x8"])
    pc.set_defaults(fn=lambda a: complete(a.stage))
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
