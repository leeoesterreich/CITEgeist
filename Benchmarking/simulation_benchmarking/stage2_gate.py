"""Stage 2 gate for a chosen brightness estimator on the calibrated sim (plan Task 20).

Per condition, MEDIAN across replicates: Spearman(b, b*) >= 0.5 across the sim's markers AND per-type mean
proportion r >= plain - 0.01. Every condition must pass; any failure (or undefined/NaN value) => {"outcome": "drop"}.
This is a recovery/stress test on a sim built from the same model family, NOT independent validation.

Plain-only path (Stage 1 outcome = drop, no chosen estimator): the gate is NOT evaluated; the report carries
plain_only=True and says no chosen arm exists.

per_rep columns: condition, replicate, spearman_b_bstar, props_r (chosen arm), props_r_plain.
CLI: stage2_gate.py --out FILE [--method M --chosen-results-dir D --plain-results-dir D --sim-root R]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts"))

SPEARMAN_MIN = 0.5
PLAIN_MARGIN = 0.01
COLUMNS = ["condition", "replicate", "spearman_b_bstar", "props_r", "props_r_plain"]
NOTE_PLAIN_ONLY = (
    "Stage 1 selected no estimator (outcome drop): no chosen arm exists, so the Stage 2 estimator gate was not evaluated. "
    "Only the plain marker-share arm was run on the calibrated sim."
)
NOTE_GATED = (
    "Stage 2 is a recovery/stress test on a sim built from the same model family, not independent validation."
)


def _median(values) -> float:
    """Median across replicates; NaN counts as failing (-inf) so an undefined replicate can never help a condition pass."""
    v = np.asarray(values, dtype=float)
    return float(np.median(np.where(np.isfinite(v), v, -np.inf)))


def stage2_gate(per_rep: pd.DataFrame) -> dict:
    missing = [c for c in COLUMNS if c not in per_rep.columns]
    if missing:
        raise ValueError(f"per_rep is missing columns {missing}")
    if len(per_rep) == 0:
        raise ValueError("per_rep is empty; refusing a vacuous pass")
    per_condition = {}
    for cond, g in per_rep.groupby("condition"):
        sp = _median(g["spearman_b_bstar"])
        r = _median(g["props_r"])
        plain = _median(g["props_r_plain"])
        reasons = []
        if not np.isfinite(sp):
            reasons.append("median Spearman undefined")
        if not np.isfinite(plain):
            reasons.append("median plain props r undefined")
        if not np.isfinite(r):
            reasons.append("median candidate props r undefined")
        ok_sp = bool(np.isfinite(sp) and sp >= SPEARMAN_MIN)
        ok_r = bool(np.isfinite(r) and np.isfinite(plain) and r >= plain - PLAIN_MARGIN)
        per_condition[str(cond)] = {
            "n_replicates": int(len(g)),
            "median_spearman_b_bstar": sp,
            "median_props_r": r,
            "median_props_r_plain": plain,
            "spearman_pass": ok_sp,
            "props_pass": ok_r,
            "pass": ok_sp and ok_r,
            "fail_reasons": reasons,
        }
    ok = all(v["pass"] for v in per_condition.values())
    return {
        "outcome": "pass" if ok else "drop",
        "thresholds": {"spearman_min": SPEARMAN_MIN, "plain_margin": PLAIN_MARGIN},
        "per_condition": per_condition,
    }


def build_report(per_rep, chosen_method: str) -> dict:
    """chosen_method empty => plain-only report (no gate). Otherwise per_rep is required and the gate is embedded."""
    if not chosen_method:
        return {"plain_only": True, "chosen_method": None, "note": NOTE_PLAIN_ONLY}
    if per_rep is None:
        raise ValueError("a chosen estimator needs the per-replicate table to evaluate the gate")
    return {"plain_only": False, "chosen_method": chosen_method, "stage2_gate": stage2_gate(per_rep), "note": NOTE_GATED}


def write_report(path, report: dict) -> dict:
    import provenance

    return provenance.write_stamped_json(path, report, "Benchmarking/simulation_benchmarking/stage2_gate.py")


def _per_type_mean_r(pred_csv: Path, gt_csv: Path) -> float:
    """One replicate's mean over the defined types of the per-type r. Per-type r is score_calsim_per_type.per_type_r (the
    owner: GT constant = undefined, prediction constant = 0). This per-replicate scalar is the gate's own quantity (this
    gate is for the retired brightness estimator); the reported per-type r aggregates type-first
    (score_calsim_per_type.aggregate_per_type_r)."""
    from score_calsim_per_type import _mean_defined, per_type_r

    gt = pd.read_csv(gt_csv, index_col=0).drop(columns=["spot_x", "spot_y"], errors="ignore")
    pr = pd.read_csv(pred_csv, index_col=0).drop(columns=["spot_x", "spot_y", "recon_error"], errors="ignore").loc[gt.index]
    return _mean_defined(per_type_r(pr, gt[[t for t in gt.columns if t in pr.columns]]))


def _spearman(a: pd.Series, b: pd.Series) -> float:
    from scipy.stats import spearmanr

    common = a.index.intersection(b.index)
    if len(common) < 3 or a[common].nunique() < 2 or b[common].nunique() < 2:
        return float("nan")
    return float(spearmanr(a[common], b[common])[0])


def per_rep_table(sim_root: Path, chosen_dir: Path, plain_dir: Path) -> pd.DataFrame:
    """Build the gate table from run outputs. b = applied brightness from the chosen arm's diagnostics (used_brightness if present),
    b* = the sim's per-marker b* from replicates_calsim/calsim_reports. Proportion r = per-type mean Pearson r vs GT."""
    rows = []
    for cond in ("high_seg", "mixed"):
        for rep in range(5):
            gt = sim_root / cond / "ST_sim" / f"Wu_ST_{rep}_prop.csv"
            unit = f"sim_6_{cond}_{rep}"
            props = f"{cond}/Wu_rep_{rep}/props.csv"
            bstar = json.loads((sim_root / "calsim_reports" / f"bstar_{cond}_{rep}.json").read_text())["report"]["sims"][f"{cond}_{rep}"]["bstar"]
            diag = json.loads((chosen_dir / unit / "diagnostics" / "baseline.json").read_text())["brightness"]
            b = {m: (d.get("used_brightness", d.get("applied_b")) if isinstance(d, dict) else d) for m, d in diag.items()}
            b = pd.Series({m: float(v) for m, v in b.items() if v is not None})
            rows.append({
                "condition": cond, "replicate": rep,
                "spearman_b_bstar": _spearman(b, pd.Series(bstar, dtype=float)),
                "props_r": _per_type_mean_r(chosen_dir / unit / props, gt),
                "props_r_plain": _per_type_mean_r(plain_dir / unit / props, gt),
            })
    return pd.DataFrame(rows, columns=COLUMNS)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--method", default="", help="chosen estimator; empty = plain-only report (gate not evaluated)")
    ap.add_argument("--sim-root")
    ap.add_argument("--chosen-results-dir")
    ap.add_argument("--plain-results-dir")
    a = ap.parse_args(argv)
    out = Path(a.out)
    if out.exists():
        raise SystemExit(f"{out} exists; refusing to overwrite")
    table = None
    if a.method:
        table = per_rep_table(Path(a.sim_root), Path(a.chosen_results_dir), Path(a.plain_results_dir))
    rep = build_report(table, a.method)
    if table is not None:
        rep["per_rep"] = table.to_dict(orient="records")
    write_report(out, rep)
    print(json.dumps({k: v for k, v in rep.items() if k != "per_rep"}, indent=1))


if __name__ == "__main__":
    main()
