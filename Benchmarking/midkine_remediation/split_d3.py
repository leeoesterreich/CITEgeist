"""B22 part 4 (check 13): which part of check 11's D3 removes the MDK-secretory Moran's I, qp_prop or RNA per cell? Check 11 with a different form list and projector
(forms none, D1, D3, D4..D7). Descriptive only.

Promoted from the exploratory check 13 (``exploratory_v42/check13_split_d3.py``, git history). The analysis functions are that code, unchanged (``make_projector13``, ``label``,
``gate_vs_check11``, ``spearmans``, ``tables``; tests/midkine_remediation/test_moran_robustness.py checks the AST). Changed: the check 11 summary file the exploratory check gated
against is the check 11 result computed earlier IN THE SAME RUN (``c11``, an argument), ``main`` became ``compute`` (check 11's module globals ``FORMS`` and ``make_projector`` are
patched for the duration of the call and restored after, so the checks that run later in the same job see check 11 unchanged), the table header no longer prints a dirty flag.
Rules: ``exploratory_v42/check13_DECISION_RULE.md`` (git history).
"""
from __future__ import annotations

import contextlib
import json
import os
import sys
import tempfile
import types
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import a12_10_real as R  # noqa: E402
import depth_removed_moran as C11  # noqa: E402
import hurdle as H  # noqa: E402
import promoted_common as CM  # noqa: E402

FORMS13 = ("none", "D1", "D3", "D4", "D5", "D6", "D7", "D5_percell")
RNA_COVARIATE = "total"            # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md): D3-D7 use RNA content (rna_total)
# columns of hurdle.covariate_matrix(obs): 0 log_spot_umi, 1 log_rna_total (log_rna_per_cell_equiv for D5_percell), 2 qp_prop
DESIGN = {"D3": (0, 1, 2), "D4": (0, 2), "D5": (0, 1), "D6": (2,), "D7": (1,), "D5_percell": (0, 1)}
PERCELL_FORMS = {"D5_percell": "per_cell"}   # diagnostic arm: the pre-deviation D5 (depth + RNA per cell), kept so the move of D5 is shown beside it (gated against the stored v2b22 D5 in moran_robustness)
S1_SIZE = 0.048217128908550164
_orig_projector = C11.make_projector


def make_projector13(form: str, obs, depth):
    """check 11's projector for none/D1/D2; D3..D7 = intercept + the named columns of hurdle.covariate_matrix(obs, rna=RNA_COVARIATE); D5_percell = the same columns of the per-cell matrix."""
    if form not in DESIGN:
        return _orig_projector(form, obs, depth)
    Z = H.covariate_matrix(obs, rna=PERCELL_FORMS.get(form, RNA_COVARIATE))[:, list(DESIGN[form])]
    return C11.basis(np.column_stack([np.ones(len(depth)), Z]))


def label(m: dict) -> dict:
    """Mechanical reading rule. m = cohort mean real I per form (needs D1, D4, D5)."""
    r4, r5 = 1 - m["D4"] / m["D1"], 1 - m["D5"] / m["D1"]
    hi4, hi5, lo4, lo5 = r4 >= 2 / 3, r5 >= 2 / 3, r4 < 1 / 3, r5 < 1 / 3
    if hi4 and hi5:
        lab = "BOTH"
    elif hi4 and lo5:
        lab = "TUMOUR_CONTENT"
    elif hi5 and lo4:
        lab = "RNA_PER_CELL"
    elif not hi4 and not hi5:
        lab = "NEITHER_ALONE"
    else:
        lab = "UNCLASSIFIED"
    return {"label": lab, "removal_D4": r4, "removal_D5": r5}


def gate_vs_check11(spec: str, res: dict, c11: dict) -> None:
    """D1 and D3 must reproduce check 11's stored I and p (and status) within 1e-9, else fail loudly."""
    for fm in ("D1", "D3"):
        got, exp = res["forms"][fm], c11["per_specimen"][spec]["forms"][fm]
        if got.get("status") != exp.get("status"):
            raise CM.GateFailure(f"GATE FAILED {spec} {fm} status {got.get('status')} vs check11 {exp.get('status')}")
        if exp.get("status") == "ok":
            for k in ("I", "p"):
                CM.gate_equal(f"{spec} {fm} {k} vs check11 stored", got[k], exp[k], rtol=1e-9, atol=1e-9)


def spearmans(obs) -> dict:
    d = {"log_spot_umi": np.log(obs["spot_umi"].to_numpy(float)), "log_rna_total": np.log(H.rna_total(obs)), "qp_prop": obs["qp_prop"].to_numpy(float)}
    pairs = (("qp_prop", "log_rna_total"), ("qp_prop", "log_spot_umi"), ("log_rna_total", "log_spot_umi"))
    return {f"{a}~{b}": float(spearmanr(d[a], d[b])[0]) for a, b in pairs}


def tables(meta, per, coh, b9) -> str:
    f = C11.f
    L = [f"# Check 13 generated tables (n_real {meta['n_real']}, B power {meta['b_power']}, controls {meta['n_ctrl']})", "", f"Output dir `{meta['out']}`; commit {meta['git_commit'][:8]}; gates passed for {meta['n_gates']} specimens.", "",
         "## Reading rule", "", "```", json.dumps(coh["reading"], indent=1, default=float), "```", "",
         "## Part 1: I (MSR p, BH q within form)", "", "| specimen | n | stored I (q) | " + " | ".join(FORMS13) + " |", "|---|---|---|" + "---|" * len(FORMS13)]
    for s, r in per.items():
        if r.get("status") != "ok":
            L.append(f"| {s} | {r.get('n')} | not estimable: {r.get('reason')} |" + "|" * len(FORMS13))
            continue
        cell = lambda v: "not est." if v.get("status") != "ok" else f"{f(v['I'])} (p {f(v['p'], 3)}, q {f(v.get('q_bh'))})"  # noqa: E731
        L.append(f"| {s} | {r['n']} | {f(b9['per_specimen'][s]['I'])} ({f(b9['bh'][s]['q'])}) | " + " | ".join(cell(r["forms"][fm]) for fm in FORMS13) + " |")
    L += ["", "## Spearman among covariates, per specimen", "", "| specimen | qp~rna_total | qp~log_umi | rna_total~log_umi |", "|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") == "ok":
            sp = r["spearman"]
            L.append(f"| {s} | {f(sp['qp_prop~log_rna_total'], 3)} | {f(sp['qp_prop~log_spot_umi'], 3)} | {f(sp['log_rna_total~log_spot_umi'], 3)} |")
    L += ["", "## Part 2: controls per form", "", "| specimen | form | real I | a mean | a p95 | a pct | b mean | b p95 | b pct | >p95 a | >p95 b |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in FORMS13:
            v = r["forms"][fm]
            if v.get("status") == "ok":
                c = v["controls"]
                L.append(f"| {s} | {fm} | {f(v['I'])} | {f(c['a']['mean'])} | {f(c['a']['p95'])} | {f(c['a']['pct'], 1)} | {f(c['b']['mean'])} | {f(c['b']['p95'])} | {f(c['b']['pct'], 1)} | {c['a']['gt95']} | {c['b']['gt95']} |")
    L += ["", "## Part 1b: depth after removal", "", "| specimen | form | r(depth,x) | I(depth,x) | q | r(depth,y) | I(depth,y) | q |", "|---|---|---|---|---|---|---|---|"]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for fm in FORMS13:
            v = r["forms"][fm]
            if v.get("status") == "ok":
                dx, dy = v["depth"]["x"], v["depth"]["y"]
                L.append(f"| {s} | {fm} | {f(dx['pearson_r'])} | {f(dx.get('I'))} | {f(dx.get('q_bh'))} | {f(dy['pearson_r'])} | {f(dy.get('I'))} | {f(dy.get('q_bh'))} |")
    L += ["", "## Power: mean I [combined p] (frac p<0.05)", "", "| specimen | scenario | " + " | ".join(FORMS13) + " |", "|---|---|" + "---|" * len(FORMS13)]
    for s, r in per.items():
        if r.get("status") != "ok":
            continue
        for sc in C11.SCEN:
            blk = r["power"]["scenarios"][sc]
            if "forms" not in blk:
                L.append(f"| {s} | {sc} | {blk.get('status')} |" + "|" * (len(FORMS13) - 1))
                continue
            cells = [f"not est. ({blk['forms'][fm].get('n_valid', 0)} valid)" if blk["forms"][fm]["status"] != "ok" else
                     f"{f(blk['forms'][fm]['mean'])} [{f(blk['forms'][fm]['p_combined'], 3)}] ({f(blk['forms'][fm]['frac_p_lt_0.05'], 2)})" for fm in FORMS13]
            L.append(f"| {s} | {sc} | " + " | ".join(cells) + " |")
    L += ["", "## Cohort", "", "```", json.dumps(coh, indent=1, default=float, sort_keys=True), "```", ""]
    return "\n".join(L)


@contextlib.contextmanager
def patched_check11():
    """check 11's run_specimen / cohort / csv writers read the module globals ``FORMS`` and ``make_projector``: swap them for check 13's and restore on exit."""
    keep = (C11.FORMS, C11.make_projector)
    C11.FORMS, C11.make_projector = FORMS13, make_projector13
    try:
        yield
    finally:
        C11.FORMS, C11.make_projector = keep


def compute(comp_run, b9: dict, b17: dict, c11: dict, specs, out, meta_common: dict, n_real: int = C11.N_REAL, n_ctrl: int = C11.N_CTRL, b_power: int = R.B) -> tuple:
    """Check 13 (was ``main``). ``c11`` = the check 11 summary computed earlier in the same run (``{"per_specimen": ...}``); D1 and D3 are gated against it (I and p to 1e-9).
    Writes check13_{part1,controls,depthcheck,power}.csv and check13_tables.md into ``out``. Returns ``(summary, files)``; raises ``GateFailure``."""
    a = types.SimpleNamespace(n_real=n_real, n_ctrl=n_ctrl, b_power=b_power, comp_run=comp_run)
    per = {}
    with patched_check11():
        for s in specs:
            per[s] = C11.run_specimen(s, a, b9, b17)
            if per[s].get("status") == "ok":
                gate_vs_check11(s, per[s], c11)
                per[s]["spearman"] = spearmans(R.load_cgeist_eligible(comp_run, s)["obs"])
        coh = C11.cohort(per, b9, b17)
        ok = {fm: {s: r["forms"][fm]["I"] for s, r in per.items() if r.get("status") == "ok" and r["forms"][fm].get("status") == "ok"} for fm in FORMS13}
        m = {fm: float(np.mean([ok[fm][s] for s in ok["D1"] if s in ok[fm]])) for fm in FORMS13 if ok[fm]}
        coh["reading"] = {**label(m), "cohort_mean_I_over_specimens_estimable_in_D1": m,
                          "removal_by_form": {fm: 1 - m[fm] / m["D1"] for fm in FORMS13 if fm in m},
                          "side_removal_D6_D7": {fm: 1 - m[fm] / m["D1"] for fm in ("D6", "D7") if fm in m}, "s1_size": S1_SIZE}
        sp = np.array([[r["spearman"][k] for k in r["spearman"]] for r in per.values() if r.get("status") == "ok"])
        coh["spearman"] = {k: {"mean": float(sp[:, i].mean()), "min": float(sp[:, i].min()), "max": float(sp[:, i].max())} for i, k in enumerate(next(r["spearman"] for r in per.values() if r.get("status") == "ok"))}
        meta = {"check": "check13_split_d3", "n_real": n_real, "n_ctrl": n_ctrl, "b_power": b_power, "seed": C11.SEED11, "k": R.K,
                "n_gates": sum(r.get("status") == "ok" for r in per.values()), **meta_common}
        out = Path(out)
        with tempfile.TemporaryDirectory(dir=out) as tmp:                                  # check 11's writer names its files check11_*
            C11.write_csvs(Path(tmp), per)
            for nm in ("part1", "controls", "depthcheck", "power"):
                os.replace(Path(tmp) / f"check11_{nm}.csv", out / f"check13_{nm}.csv")
        (out / "check13_tables.md").write_text(tables(meta, per, coh, b9))
    return {**meta, "per_specimen": per, "cohort": coh}, [f"check13_{nm}.csv" for nm in ("part1", "controls", "depthcheck", "power")] + ["check13_tables.md"]
