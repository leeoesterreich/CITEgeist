"""Old-vs-new gate for the RNA-total covariate reruns (protocol deviation 2026-11-06, plan rev 2 / rev 3 gate lists).

Each pair is (superseded per-cell run, rerun with rna="total"): B7 v2b7/v2b71, B8 v2b8/v2b81, B18 v2b18/v2b181, B12 v2b121/v2b122, B23 v2b23/v2b231, B21 v2b21/v2b211
(B221 against v2b22 is gated inside the stage; B24 has no old run). FAILS (exit 1) unless every field that must not move is identical, then prints what moved (old -> new)
for the adjusted rows. Fields that may move are listed per stage in ``MUST_NOT_MOVE`` / ``MAY_MOVE`` below; anything that differs and is in neither rule fails closed.

    python compare_rna_total_reruns.py [--root <namespace root>]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SHA = "fc43dbea67bf"
HEADLINE = "HCC22-088-P4-S2_1i_rep"
RES = f"Benchmarking/midkine_remediation/results/{SHA}"
PC = f"output/positive_control_module3_7type/{SHA}"
PAIRS = {"B7": ("v2b7", "v2b71"), "B8": ("v2b8", "v2b81"), "B18": ("v2b18", "v2b181"), "B12": ("v2b121", "v2b122"), "B23": ("v2b23", "v2b231"), "B21": ("v2b21", "v2b211")}
# run metadata: differs by construction in every rerun (run id, commit, payload hash, paths, deviation record, covariate mode and names)
META = {"_code_commit", "_code_commit_source", "_code_dirty", "_payload_sha256", "_produced_by", "run_id", "analysis_id", "git_commit", "hurdle_reference", "hurdle_run", "keratin_run",
        "protocol_deviation", "rna_covariate", "covariates", "covariates_partial", "coverage_source", "sidecar_files_sha256", "files_sha256", "blas_threads", "adjustments"}
MUST_NOT_MOVE = {
    "B7": "label, n_injected_spots, n_lr_pairs, n_eligible, n_compartments, fold, dis_thr, pair (baseline, injected, report, value may move)",
    "B8": "input stamp hashes except the covariate-dependent inputs (B7 positive control, B24 coverage), per specimen (12 h1 files + the A01-A07 sweep): n_eligible, n_blocks, label, "
          "status / reason_code / n_active / n_active_high / n_active_low / n_mixed_blocks of every run and variant; per_specimen n_eligible and status; H4 file",
    "B18": "everything after dropping run id, paths, rna_covariate, commit and stamp fields; the csv byte for byte",
    "B12": "unadjusted primary and secondary (per specimen, cohort, focal, descriptive); only the partial_* rows may move",
    "B23": "the adjustments whose label has no RNA covariate (unadjusted, log_spot_umi, qp_prop, log_spot_umi+qp_prop): status, rho, p, n_surrogates",
    "B21": "every primary row and table row (check 7 and 8), every activity / score field; only the adjusted rows, the A07 sweep copy and the cohort lines built from them may move",
}


def load(p):
    return json.loads(Path(p).read_text())


def leaves(o, p=""):
    if isinstance(o, dict):
        for k in o:
            yield from leaves(o[k], f"{p}/{k}")
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, f"{p}[{i}]")
    else:
        yield p, o


def same(x, y) -> bool:
    if isinstance(x, float) and isinstance(y, float):
        return (math.isnan(x) and math.isnan(y)) or x == y
    return x == y


def diffs(a, b) -> list:
    """[(path, old, new)] over the union of leaf paths (a path present on one side only counts as a difference)."""
    la, lb = dict(leaves(a)), dict(leaves(b))
    return [(p, la.get(p, "<absent>"), lb.get(p, "<absent>")) for p in sorted(set(la) | set(lb)) if p not in la or p not in lb or not same(la[p], lb[p])]


def is_meta(path: str) -> bool:
    return any(seg.split("[")[0] in META for seg in path.split("/") if seg)


class Gate:
    def __init__(self):
        self.problems, self.moved, self.n_checked = [], [], 0

    def identical(self, what, a, b):
        self.n_checked += sum(1 for _ in leaves(a))
        d = diffs(a, b)
        self.problems += [f"{what}{p}: {x!r} -> {y!r}" for p, x, y in d[:5]] + ([f"{what}: ... {len(d)} differences" ] if len(d) > 5 else [])

    def pick(self, o, keys):
        return {k: o.get(k) for k in keys}


def b7(root, g):
    o, n = (load(root / PC / v / HEADLINE / "positive_control.json") for v in PAIRS["B7"])
    keys = ("label", "n_injected_spots", "n_lr_pairs", "n_eligible", "n_compartments", "fold", "dis_thr", "pair", "control", "estimator", "gates_branch")
    g.identical("B7 ", g.pick(o, keys), g.pick(n, keys))
    for part in ("baseline", "injected"):
        g.moved.append((f"B7 {part} gmr", o[part].get("gmr"), n[part].get("gmr")))
    g.moved += [("B7 value (delta log GMR)", o.get("value"), n.get("value")), ("B7 status", o.get("status"), n.get("status"))]
    assert n.get("rna_covariate") == "total", "B7 rerun does not record rna_covariate=total"


RUN_KEYS = ("status", "reason_code", "n_active", "n_active_high", "n_active_low", "n_mixed_blocks")


def _est(o):
    return {r: {part: {k: o["runs"][r][part].get(k) for k in RUN_KEYS} for part in ("magnitude", "participation")} for r in o["runs"]}


def b8(root, g):
    ro, rn = root / RES / PAIRS["B8"][0], root / RES / PAIRS["B8"][1]
    so, sn = load(f"{ro}/commot_compartment_hurdle.json.stamp.json")["input_sha256"], load(f"{rn}/commot_compartment_hurdle.json.stamp.json")["input_sha256"]
    dep = lambda p: "positive_control" in p or "prod_coverage" in p   # noqa: E731  (the covariate-dependent inputs)
    so, sn = {k: v for k, v in so.items() if not dep(k)}, {k: v for k, v in sn.items() if not dep(k)}
    g.identical("B8 input stamp hashes ", so, sn)
    if not so:
        g.problems.append("B8: no comparable stamp inputs")
    specs = sorted(p.name for p in ro.iterdir() if p.is_dir())
    assert len(specs) == 12, specs
    for s in specs:
        ho, hn = load(f"{ro}/{s}/h1/h1_{s}.json"), load(f"{rn}/{s}/h1/h1_{s}.json")
        g.identical(f"B8 {s} h1 ", {"n_eligible": ho["n_eligible"], "n_blocks": ho["n_blocks"], "label": ho["label"], "runs": _est(ho)},
                    {"n_eligible": hn["n_eligible"], "n_blocks": hn["n_blocks"], "label": hn["label"], "runs": _est(hn)})
        if s == HEADLINE:
            so_, sn_ = load(f"{ro}/{s}/sweep/sweep_{s}.json"), load(f"{rn}/{s}/sweep/sweep_{s}.json")
            for v in so_["variants"]:
                f = lambda d: {"n_eligible": d.get("n_eligible"), "n_blocks": d.get("n_blocks"), "label": d.get("label"),   # noqa: E731
                               **{part: {k: d[part].get(k) for k in RUN_KEYS} for part in ("magnitude", "participation")}}
                g.identical(f"B8 sweep {v} ", f(so_["variants"][v]), f(sn_["variants"][v]))
            g.moved += [(f"B8 sweep {v} gmr", so_["variants"][v]["magnitude"].get("gmr"), sn_["variants"][v]["magnitude"].get("gmr")) for v in ("A01", "A07")]
            mo, mn = ho["runs"]["R00"], hn["runs"]["R00"]
            g.moved += [("B8 P4-S2 R00 GMR", mo["magnitude"].get("gmr"), mn["magnitude"].get("gmr")),
                        ("B8 P4-S2 R00 GMR interval lo", mo["magnitude"]["interval"].get("lo"), mn["magnitude"]["interval"].get("lo")),
                        ("B8 P4-S2 R00 GMR interval hi", mo["magnitude"]["interval"].get("hi"), mn["magnitude"]["interval"].get("hi")),
                        ("B8 P4-S2 R00 participation OR", mo["participation"].get("odds_ratio"), mn["participation"].get("odds_ratio")),
                        ("B8 P4-S2 R11 GMR", ho["runs"]["R11"]["magnitude"].get("gmr"), hn["runs"]["R11"]["magnitude"].get("gmr"))]
    co, cn = load(f"{ro}/commot_compartment_hurdle.json"), load(f"{rn}/commot_compartment_hurdle.json")
    g.identical("B8 per_specimen ", co["per_specimen"], cn["per_specimen"])
    ho4, hn4 = load(f"{ro}/h4_all.json"), load(f"{rn}/h4_all.json")
    g.identical("B8 H4 ", {k: v for k, v in ho4.items() if not is_meta("/" + k)}, {k: v for k, v in hn4.items() if not is_meta("/" + k)})
    g.moved += [("B8 claim tree magnitude branch", (co.get("claim_tree") or {}).get("magnitude", {}).get("branch"), (cn.get("claim_tree") or {}).get("magnitude", {}).get("branch")),
                ("B8 claim tree participation branch", (co.get("claim_tree") or {}).get("participation", {}).get("branch"), (cn.get("claim_tree") or {}).get("participation", {}).get("branch")),
                ("B8 negative control observed GMR", co["controls"]["negative_control"].get("observed", {}).get("gmr"), cn["controls"]["negative_control"].get("observed", {}).get("gmr")),
                ("B8 negative control median", (co["controls"]["negative_control"].get("summary") or {}).get("median"), (cn["controls"]["negative_control"].get("summary") or {}).get("median"))]
    assert cn.get("rna_covariate") == "total" and cn.get("coverage_source", {}).get("stage") == "B24", "B8 rerun does not record total mode and the B24 coverage"


def b18(root, g):
    ro, rn = root / RES / PAIRS["B18"][0], root / RES / PAIRS["B18"][1]
    o, n = load(ro / "celltype_map.json"), load(rn / "celltype_map.json")
    f = lambda d: {k: v for k, v in d.items() if not is_meta("/" + k)}   # noqa: E731
    g.identical("B18 ", f(o), f(n))
    if (ro / "celltype_map.csv").read_bytes() != (rn / "celltype_map.csv").read_bytes():
        g.problems.append("B18 celltype_map.csv differs")
    assert n.get("rna_covariate") == "total", "B18 rerun does not record the B8 covariate mode"


def b12(root, g):
    o, n = (load(root / RES / v / "keratin_secretory_summary.json") for v in PAIRS["B12"])
    d = diffs(o, n)
    for p, x, y in d:
        if is_meta(p) or any(seg.split("[")[0].startswith("partial") for seg in p.split("/")):
            continue
        g.problems.append(f"B12 {p}: {x!r} -> {y!r}")
    g.n_checked += sum(1 for _ in leaves(o))
    P = HEADLINE
    for t in ("partial_primary", "partial_secondary"):
        for k in ("rho", "p_two_sided", "q"):
            g.moved.append((f"B12 P4-S2 {t} {k}", o["focal_P4_S2"][t].get(k), n["focal_P4_S2"][t].get(k)))
    assert P in o["per_specimen"]["primary"]
    for t in ("primary", "secondary"):   # explicit non-vacuity: the unadjusted rows exist and are compared
        assert o["per_specimen"][t] and not [1 for p_, _x, _y in d if p_.startswith(f"/per_specimen/{t}")]
    assert n["config"].get("rna_covariate") == "total"


def _rename(label: str) -> str:
    return label.replace("log_rna_per_cell_equiv", "log_rna_total")


def b23(root, g):
    o, n = (load(root / RES / v / "keratin_covariate_results.json") for v in PAIRS["B23"])
    assert n["config"].get("rna_covariate") == "total"
    n_same = 0
    for s, ro in o["per_specimen"].items():
        rn = n["per_specimen"][s]
        g.identical(f"B23 {s} head ", {k: ro.get(k) for k in ("status", "n", "reason")}, {k: rn.get(k) for k in ("status", "n", "reason")})
        if ro.get("status") != "ok":
            continue
        ao, an = ro["adjustments"], rn["adjustments"]
        if sorted(map(_rename, ao)) != sorted(an):
            g.problems.append(f"B23 {s}: adjustment labels differ after renaming")
            continue
        for lab, rec in ao.items():
            new = an[_rename(lab)]
            if "log_rna" in lab:
                if s == HEADLINE:
                    g.moved += [(f"B23 P4-S2 {_rename(lab)} {k}", rec.get(k), new.get(k)) for k in ("rho", "p_two_sided")]
            else:
                g.identical(f"B23 {s} {lab} ", {k: rec.get(k) for k in ("status", "rho", "p_two_sided", "n_surrogates")}, {k: new.get(k) for k in ("status", "rho", "p_two_sided", "n_surrogates")})
                n_same += 1
    if n_same < 4 * 10:
        g.problems.append(f"B23: only {n_same} unadjusted rows compared (expected 4 labels x at least 10 specimens)")
    g.moved.append(("B23 leading single covariate", o["summary"]["leading_single_covariate"], n["summary"]["leading_single_covariate"]))


def b21(root, g):
    o, n = (load(root / RES / v / "e2_covariation_results.json") for v in PAIRS["B21"])
    assert n["config"].get("rna_covariate") == "total"
    for p, x, y in diffs(o, n):
        if is_meta(p) or "adjusted" in p or "rank_by_rho" in p or "/a07_as_stored/" in p or "/per_specimen_table[" in p:   # table rows: gated row by row below (primary identical, adjusted may move)
            continue
        g.problems.append(f"B21 {p}: {x!r} -> {y!r}")
    g.n_checked += sum(1 for _ in leaves(o))
    for c in ("check7", "check8"):
        ro, rn = o[c]["per_specimen_table"], n[c]["per_specimen_table"]
        if len(ro) != len(rn) or not ro:
            g.problems.append(f"B21 {c} table length {len(ro)} vs {len(rn)}")
            continue
        for x, y in zip(ro, rn):
            key = lambda r: (r["specimen"], r["test"], r["form"], r.get("variant", r.get("part")))   # noqa: E731
            if key(x) != key(y):
                g.problems.append(f"B21 {c}: row order/key differs {key(x)} vs {key(y)}")
            elif x["form"] != "adjusted":
                g.identical(f"B21 {c} {key(x)} ", x, y)
            elif x["specimen"] == HEADLINE and x["test"] in ("K2", "E1", "P1", "L1", "main_mdk", "kpres_vs_mdk_np", "B_kpres_vs_mdk_np", "mdk_vs_esr1") or (c == "check8" and x["specimen"] == HEADLINE and x.get("variant") == "collectri"):
                g.moved.append((f"B21 {c} P4-S2 {x['test']} adjusted rho", x.get("rho"), y.get("rho")))
    g.moved.append(("B21 A07 sweep GMR (copy of the B8 sweep)", o["check7"]["a07_as_stored"]["magnitude"].get("gmr"), n["check7"]["a07_as_stored"]["magnitude"].get("gmr")))


STAGES = {"B7": b7, "B8": b8, "B18": b18, "B12": b12, "B23": b23, "B21": b21}


def fmt(v):
    return f"{v:.4g}" if isinstance(v, float) else str(v)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=REPO)
    ap.add_argument("--stage", action="append", choices=sorted(STAGES), help="default: all")
    a = ap.parse_args(argv)
    rc = 0
    for st in a.stage or list(STAGES):
        g = Gate()
        try:
            STAGES[st](a.root, g)
        except (AssertionError, KeyError, FileNotFoundError) as e:
            g.problems.append(f"{st}: {type(e).__name__}: {e}")
        old, new = PAIRS[st]
        print(f"== {st} {old} -> {new}: must not move: {MUST_NOT_MOVE[st]}")
        print(f"   compared {g.n_checked} values; {'OK, nothing that must stay put moved' if not g.problems else 'FAILED'}")
        for p in g.problems[:20]:
            print("   PROBLEM", p)
        for what, x, y in g.moved:
            print(f"   moved  {what:55s} {fmt(x):>14s} -> {fmt(y)}")
        rc |= bool(g.problems)
    return rc


if __name__ == "__main__":
    sys.exit(main())
