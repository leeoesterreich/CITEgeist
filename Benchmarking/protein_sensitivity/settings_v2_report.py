"""Aggregate settings round 2: OAT detection settings (production + benchmark-only families), joint factorial LORO.

CPU only. FAIL-CLOSED: every prerequisite gate must pass before ANY delta, LORO selection, headline or verdict is
computed; on failure only settings_v2_gate_failure.json is written. Gates: frozen-protocol run provenance (protocol
md5, authorized implementation commit, conditions hash, no conflicting duplicate candidates, never mtime), exact
row count and condition-id set, baseline == A11 gate (1e-4), all R1 duplicate cells present and equal (1e-4).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from protein_sensitivity import aggregate as A  # noqa: E402
from protein_sensitivity import loro as L  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402
from protein_sensitivity import settings_conditions as SC  # noqa: E402

RESULTS = Path(__file__).resolve().parent / "results"
MASK_AFFECTING = {"fusion_mode", "gex_detection_k", "gex_detection_min_corr", *SC.DETECTION_VARIANT_PARAMS}
PANELS = (6, 7)
HEX40 = re.compile(r"[0-9a-f]{40}")


def _cid(rec: dict) -> str:
    return "baseline" if rec["axis"] == "baseline" else R.condition_id(rec)


def expected_ids(proto: dict) -> list[str]:
    return sorted(_cid(c) for c in SC.build_list_v2(proto))


def conditions_sha256(proto: dict) -> str:
    """Hash of the conditions file the submit script writes (json.dump of build_list_v2)."""
    return hashlib.sha256(json.dumps(SC.build_list_v2(proto)).encode()).hexdigest()


def _sha_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ------------------------------------------------------------------ run selection (never by mtime)
def select_runs(root, regions, protocol_md5: str, authorized_commits, cond_sha: str, panels=PANELS):
    """Return ({(n_types, region): dir}, problems). Every selected unit must carry the frozen protocol md5, an
    authorized implementation commit and the expected conditions hash; several COMPLETE candidates for one unit are
    accepted only when their metrics.csv are byte-identical. Incomplete dirs (no run_meta/metrics) are ignored."""
    problems = []
    groups: dict[str, list[Path]] = {}
    for d in sorted(Path(root).iterdir()) if Path(root).exists() else []:
        if d.is_dir() and (d / "run_meta.json").exists() and (d / "metrics.csv").exists():
            groups.setdefault(re.sub(r"_race_.*$", "", d.name), []).append(d)
    if not authorized_commits:
        problems.append("protocol lists no authorized_implementation_commits (fail-closed)")
    bad_auth = [a for a in authorized_commits if not (isinstance(a, str) and HEX40.fullmatch(a))]
    if bad_auth:  # short prefixes and empty strings are never accepted
        problems.append(f"authorized_implementation_commits must be 40 lowercase hex chars; invalid: {bad_auth}")
        authorized_commits = [a for a in authorized_commits if a not in bad_auth]
    chosen = {}
    for nt in panels:
        for u in map(str, regions):
            cands = groups.get(f"xenium_{nt}_{u}", [])
            if not cands:
                problems.append(f"{nt}-type region {u}: no complete run directory")
                continue
            if len({_sha_file(c / "metrics.csv") for c in cands}) > 1:
                problems.append(f"{nt}-type region {u}: {len(cands)} complete candidates with DIFFERENT metrics")
                continue
            meta = json.loads((cands[0] / "run_meta.json").read_text())
            if meta.get("settings_protocol_v2_md5") != protocol_md5:
                problems.append(f"{nt}-type region {u}: protocol md5 {meta.get('settings_protocol_v2_md5')} != frozen")
            if str(meta.get("git_hash", "")) not in set(authorized_commits):
                problems.append(f"{nt}-type region {u}: git_hash {meta.get('git_hash')} not authorized")
            if meta.get("conditions_sha256") != cond_sha:
                problems.append(f"{nt}-type region {u}: conditions hash mismatch")
            chosen[(nt, int(u))] = cands[0]
    return chosen, problems


def load_diag(chosen: dict) -> dict:
    out = {}
    for (nt, reg), d in chosen.items():
        dd = d / "diagnostics"
        for f in dd.glob("*.json") if dd.exists() else []:
            out[(nt, reg, f.stem)] = json.loads(f.read_text())
    return out


# ------------------------------------------------------------------ frozen reference inputs
def load_references(proto: dict, ref_root=SC.REPO):
    """Load EXACTLY the A11 gate files and the R1 metrics file named (with sha256) in the protocol.
    Returns (gate_df, r1_df, problems). No directory scan, no mtime."""
    spec = proto.get("reference_files", {})
    problems = []
    root = Path(ref_root)

    def _load(rel, want):
        f = root / rel
        if not want or not re.fullmatch(r"[0-9a-f]{64}", str(want)):
            problems.append(f"reference {rel}: protocol sha256 missing or malformed")
            return None
        if not f.exists():
            problems.append(f"reference {rel}: file missing")
            return None
        if _sha_file(f) != want:
            problems.append(f"reference {rel}: sha256 mismatch")
            return None
        return pd.read_csv(f, dtype={"unit": str})

    gate_files = spec.get("a11_gate", {}).get("files", {})
    if not gate_files:
        problems.append("protocol lists no a11_gate reference files")
    frames = [x for rel, h in gate_files.items() if (x := _load(rel, h)) is not None]
    r1 = None
    r1_spec = spec.get("r1_metrics", {})
    if r1_spec:
        r1 = _load(r1_spec.get("path", ""), r1_spec.get("sha256"))
    else:
        problems.append("protocol lists no r1_metrics reference file")
    gate = pd.concat(frames, ignore_index=True) if frames and len(frames) == len(gate_files) else None
    return gate, r1, problems


# ------------------------------------------------------------------ gates
def r1_equivalence(df: pd.DataFrame, r1: pd.DataFrame, proto: dict, tol: float) -> dict:
    """ALL expected duplicate comparisons (count anchored to the protocol) must be present, ok, finite and equal."""
    spec = proto["r1_duplicates"]
    regions = [str(r) for r in proto["regions"]]
    base = proto["joint_block"]["baseline_cell"]
    items = [(nt, reg, p, lv) for nt in PANELS for reg in regions for p, ls in spec["levels"].items() for lv in ls]
    problems, mism = [], []
    if len(items) != spec["expected_comparisons"]:
        problems.append(f"protocol inconsistent: {len(items)} comparisons vs literal {spec['expected_comparisons']}")
    for nt, reg, p, lvl in items:
        slug = SC.joint_level(*({**base, p: lvl}[k] for k in SC.JOINT_KEYS))
        m2 = df[(df["axis"] == "joint") & (df["level"] == slug) & (df["n_types"] == nt) & (df["unit"].astype(str) == reg)]
        m1 = r1[
            (r1["axis"] == "setting")
            & (r1["param"] == p)
            & (r1["n_types"] == nt)
            & (r1["unit"].astype(str) == reg)
            & (r1["level"].map(SC.norm_level) == SC.norm_level(lvl))
        ]
        ok = len(m2) == 1 and len(m1) == 1
        if ok:
            ok = all(x["status"].iloc[0] == "ok" and np.isfinite(float(x["pearson_r"].iloc[0])) for x in (m1, m2))
        if not ok:
            problems.append(f"missing/failed/non-finite comparison {(nt, reg, p, SC.norm_level(lvl))}")
        elif abs(float(m1["pearson_r"].iloc[0]) - float(m2["pearson_r"].iloc[0])) > tol:
            mism.append((nt, reg, p, SC.norm_level(lvl)))
    return {
        "pass": not problems and not mism,
        "tol": tol,
        "n_expected": len(items),
        "mismatches": mism,
        "problems": problems,
    }


def gate_problems(df, gate, r1, proto, tol) -> tuple[list[str], dict]:
    problems, detail = [], {}
    want = proto["expected_conditions"]["total_per_region_per_ntype"]
    exp = expected_ids(proto)
    regions = [str(r) for r in proto["regions"]]
    for nt in PANELS:
        sub = df[df["n_types"] == nt]
        if len(sub) != want * len(regions):
            problems.append(f"{nt}-type: expected {want * len(regions)} rows, found {len(sub)}")
        for u in regions:
            ids = sorted(_cid(r) for r in sub[sub["unit"].astype(str) == u].to_dict("records"))
            if ids != exp:
                problems.append(f"{nt}-type region {u}: condition-id set mismatch ({len(set(ids) ^ set(exp))} ids)")
        g = gate[(gate["dataset"] == "xenium") & (gate["n_types"] == nt) & (gate["axis"] == "baseline")]
        if g["unit"].astype(str).duplicated().any() or len(g) != len(regions):
            problems.append(f"{nt}-type A11 reference: need exactly one baseline row per region, found {len(g)}")
        ref = {
            str(u): float(r) if s == "ok" else float("nan")
            for u, r, s in zip(g["unit"], g["pearson_r"], g["status"])
        }
        chk = A.check_baseline_gate(df.assign(dataset="xenium"), nt, ref, tol, regions)
        detail[f"baseline_gate_{nt}"] = chk
        if not chk["pass"]:
            problems.append(f"{nt}-type baseline gate: {chk['problems']} {chk['mismatches']}")
    eq = r1_equivalence(df, r1, proto, tol)
    detail["r1_equivalence"] = eq
    if not eq["pass"]:
        problems.append(f"R1 duplicate equivalence: {eq['problems'][:5]} {eq['mismatches'][:5]}")
    return problems, detail


# ------------------------------------------------------------------ summaries (only after all gates pass)
def _eff_thresholds(diag, nt, reg, cid):
    d = diag.get((nt, int(reg), cid), {})
    return sorted({round(c[2], 4) for c in d.get("effective_thresholds", [])})


def _deltas(rows, base, regions):
    per, bad = {}, []
    for reg in regions:
        r = rows[rows["unit"].astype(str) == reg]
        if r.empty or r["status"].iloc[0] != "ok" or not np.isfinite(r["pearson_r"].iloc[0]):
            bad.append(reg)
        else:
            per[reg] = float(r["pearson_r"].iloc[0]) - base[reg]
    return per, bad


def _robust(per, rule):
    d = list(per.values())
    return bool(abs(np.mean(d)) <= rule["mean_abs_threshold"] and max(abs(x) for x in d) <= rule["region_abs_threshold"])


def summarize_oat(df, diag, n_types, proto):
    """Per OAT detection level: deltas, inert / not-estimable status, robust yes/no (only when estimable)."""
    rule = proto["robustness_rule"]
    regions = [str(r) for r in proto["regions"]]
    sub = df[(df["n_types"] == n_types)]
    bb = sub[sub["axis"] == "baseline"]
    base = {str(u): float(r) for u, r in zip(bb["unit"], bb["pearson_r"])}
    out = []
    for param, spec in proto["settings"].items():
        for lvl in spec["levels_non_baseline"]:
            cid = R.condition_id({"axis": "setting", "param": param, "level": lvl})
            rows = sub[
                (sub["axis"] == "setting")
                & (sub["param"] == param)
                & (sub["level"].map(SC.norm_level) == SC.norm_level(lvl))
            ]
            per, bad = _deltas(rows, base, regions)
            status, robust = "ok", None
            if bad:
                status = "not_estimable: failed/missing regions " + ",".join(bad)
            elif param in MASK_AFFECTING:
                hs = [
                    (
                        diag.get((n_types, int(reg), cid), {}).get("mask_hash_fused"),
                        diag.get((n_types, int(reg), "baseline"), {}).get("mask_hash_fused"),
                    )
                    for reg in regions
                ]
                if any(a is None or b is None for a, b in hs):
                    status = "not_estimable: detection-mask hash missing"
                elif all(a == b for a, b in hs):
                    status = "inert"
            if status == "ok":
                robust = _robust(per, rule)
            d = list(per.values())
            out.append(
                {
                    "param": param,
                    "level": lvl,
                    "family": spec["family"],
                    "status": status,
                    "robust": robust,
                    "mean_delta_r": float(np.mean(d)) if d else None,
                    "worst_region_delta": float(max(d, key=abs)) if d else None,
                    "per_region_delta": per,
                    "effective_thresholds_by_region": {reg: _eff_thresholds(diag, n_types, reg, cid) for reg in regions},
                }
            )
    return out


def summarize_joint_oat(df, n_types, proto):
    """OAT levels of the four joint axes read from the joint block (the other three at baseline)."""
    rule = proto["robustness_rule"]
    jb = proto["joint_block"]
    regions = [str(r) for r in proto["regions"]]
    base = {k: SC.norm_level(jb["baseline_cell"][k]) for k in SC.JOINT_KEYS}
    sub = df[df["n_types"] == n_types]
    bb = sub[sub["axis"] == "baseline"]
    b = {str(u): float(r) for u, r in zip(bb["unit"], bb["pearson_r"])}
    out = []
    for key in SC.JOINT_KEYS:
        for lvl in jb["levels"][key]:
            if SC.norm_level(lvl) == base[key]:
                continue
            combo = {k: (lvl if k == key else jb["baseline_cell"][k]) for k in SC.JOINT_KEYS}
            slug = SC.joint_level(*(combo[k] for k in SC.JOINT_KEYS))
            per, bad = _deltas(sub[(sub["axis"] == "joint") & (sub["level"] == slug)], b, regions)
            d = list(per.values())
            out.append(
                {
                    "param": key,
                    "level": lvl,
                    "family": "production",
                    "status": "ok" if not bad else "not_estimable: failed/missing regions " + ",".join(bad),
                    "robust": _robust(per, rule) if not bad else None,
                    "mean_delta_r": float(np.mean(d)) if d else None,
                    "worst_region_delta": float(max(d, key=abs)) if d else None,
                    "per_region_delta": per,
                }
            )
    return out


def headline_counts(levels: list[dict]) -> dict:
    """'robust X of Y estimable production-setting LEVELS'; benchmark-only variants reported separately."""
    prod = [x for x in levels if x["family"] == "production" and x["status"] == "ok"]
    return {
        "phrase": "robust X of Y estimable production-setting levels",
        "production_estimable_levels": len(prod),
        "production_robust_levels": sum(bool(x["robust"]) for x in prod),
        "production_inert_or_not_estimable_levels": sum(
            x["family"] == "production" and x["status"] != "ok" for x in levels
        ),
        "benchmark_only_levels_reported_separately": sum(x["family"] != "production" for x in levels),
    }


def joint_grid(df: pd.DataFrame, proto: dict):
    jb = proto["joint_block"]
    base_slug = SC.joint_level(*(jb["baseline_cell"][k] for k in SC.JOINT_KEYS))
    rows = []
    for _, r in df[df["axis"] == "baseline"].iterrows():
        rows.append((r["n_types"], int(r["unit"]), base_slug, r["pearson_r"]))
    for _, r in df[df["axis"] == "joint"].iterrows():
        rows.append((r["n_types"], int(r["unit"]), str(r["level"]), r["pearson_r"] if r["status"] == "ok" else np.nan))
    return pd.DataFrame(rows, columns=["n_types", "region", "combo", "pearson_r"]), base_slug


def in_sample_best(grid: pd.DataFrame, base_slug: str) -> list[dict]:
    """DESCRIPTIVE only (not part of the rule): best combo on all regions and its gain over the baseline."""
    out = []
    for nt, g in grid.groupby("n_types"):
        piv = g.pivot_table(index="combo", columns="region", values="pearson_r", aggfunc="first").dropna()
        means = piv.mean(axis=1)
        best = base_slug
        for c, m in means.items():
            if m > means[best] + L.TIE_TOL:
                best = c
        out.append(
            {
                "n_types": int(nt),
                "best_combo": best,
                "mean_r": float(means[best]),
                "baseline_mean_r": float(means[base_slug]),
                "gain": float(means[best] - means[base_slug]),
            }
        )
    return out


def write_gate_failure(out: Path, problems, detail=None):
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "gate": "FAILED",
        "problems": problems,
        "detail": detail or {},
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    (out / "settings_v2_gate_failure.json").write_text(json.dumps(payload, indent=2, default=str))


def run(protocol_path, results, out, ref_root=SC.REPO) -> int:
    out = Path(out)
    if out.exists() and any(out.iterdir()):  # never write next to stale summaries
        print(f"REFUSING: output dir {out} is not empty", file=sys.stderr)
        return 2
    protocol_path = Path(protocol_path)
    if not protocol_path.exists():
        write_gate_failure(out, [f"frozen protocol {protocol_path} does not exist"])
        return 1
    proto = SC.load_protocol_v2(protocol_path)
    tol = proto["baseline"]["baseline_equals_gate_tolerance"]
    chosen, problems = select_runs(
        results,
        proto["regions"],
        SC.protocol_md5(protocol_path),
        proto.get("authorized_implementation_commits", []),
        conditions_sha256(proto),
    )
    if problems:
        write_gate_failure(out, problems)
        return 1
    df = pd.concat([pd.read_csv(d / "metrics.csv", dtype={"unit": str}) for d in chosen.values()], ignore_index=True)
    df["param"] = df["param"].fillna("")
    pat = r"a(?P<alpha>[^_]+)_r(?P<lambda_reg>[^_]+)_l(?P<lambda_laplacian>[^_]+)_g(?P<detection_gate_ub>.+)"
    for k in SC.JOINT_KEYS:
        df[k] = np.nan
    j = df["axis"] == "joint"
    parsed = df.loc[j, "level"].astype(str).str.extract(pat)
    for k in SC.JOINT_KEYS:
        df.loc[j, k] = parsed[k].astype(float)
    gate, r1, ref_problems = load_references(proto, ref_root)
    if ref_problems:
        write_gate_failure(out, ref_problems)
        return 1
    problems, detail = gate_problems(df, gate, r1, proto, tol)
    if problems:  # NOTHING else is computed or written
        write_gate_failure(out, problems, detail)
        return 1
    diag = load_diag(chosen)
    levels = {str(nt): summarize_oat(df, diag, nt, proto) + summarize_joint_oat(df, nt, proto) for nt in PANELS}
    grid, base_slug = joint_grid(df, proto)
    sel = L.select_joint(grid, base_slug)
    rule = proto["joint_block"]["decision_rule"]
    dec = L.joint_decision(sel, rule["mean_gain_max"], rule["region_gain_max"])
    decs = {int(r["n_types"]): r for r in json.loads(dec.to_json(orient="records"))}
    git_hash = R._git_hash()
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "settings_v2_metrics.csv", index=False)
    sel.to_csv(out / "settings_v2_loro_joint.csv", index=False)
    summary = {
        "stamp": {
            "report_git_hash": git_hash,
            "protocol_v2_md5": SC.protocol_md5(protocol_path),
            "conditions_sha256": conditions_sha256(proto),
            "run_git_hashes": sorted({json.loads((d / "run_meta.json").read_text())["git_hash"] for d in chosen.values()}),
            "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "n_rows": int(len(df)),
        },
        "gates": {"passed": True, **detail},
        "oat_levels": levels,
        "headline": {nt: headline_counts(v) for nt, v in levels.items()},
        "joint_loro": json.loads(sel.to_json(orient="records")),
        "primary_panel": 6,
        "primary_decision": decs[6],
        "secondary_decision_7type": decs[7],
        "in_sample_best_descriptive": in_sample_best(grid, base_slug),
        "joint_decision_rule": rule,
        "scope_statement": proto["joint_block"]["scope_statement"],
    }
    (out / "settings_v2_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"wrote {len(df)} rows; GATES PASSED")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--protocol", default=str(SC.PROTOCOL_V2))
    ap.add_argument("--results", default=str(RESULTS / "settings_v2"))
    ap.add_argument("--out", required=True, help="fresh or empty output directory")
    a = ap.parse_args(argv)
    return run(a.protocol, a.results, a.out)


if __name__ == "__main__":
    sys.exit(main())
