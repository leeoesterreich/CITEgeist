"""B18: which cell types hold the MDK-NCL sender and receiver COMMOT score. Reported in Section 2.4 (who sends, who receives). Descriptive only.

Promoted from the exploratory check 3 (``exploratory_v42/check3_celltype_map.py``, git history). ``split_type``, ``gate_p4s2``, ``stats`` and
``process`` are the exploratory functions unchanged; ``load`` only takes the COMMOT run dir as an argument instead of a module constant.
The cohort table code is the exploratory one. Removed: the comparison with the earlier protocols' G5 status (v41 runs; they are not inputs of a v42 stage).
No p-values, no ranking claims, no cross-type tests, no cell-type label comparison and no neighbour analysis.

How compartments are encoded: rows_R00.csv.gz column `id` is "<spot barcode>|<cell type>" (e.g. AACACGTGCATCGCAC-1|Cancer_Cells);
the cell type is the text after the last "|". Rows of rows_R00.csv.gz align one-to-one, in order, with R00/sender.csv.gz (checked).
R00 = marker-share proportion >= 0.05 (p_construct 0.05). Score columns s-MDK-<R> / r-MDK-<R> are the COMMOT sent / received
signal of that compartment for the pair MDK-<R>. Active = score > 0. No seed is involved.

CAVEAT: compartments of different types in the same spot share counts and are NOT independent. Shares and counts below describe
where the score sits; they are not tests and do not rank cell types.

GATE (fail closed): P4-S2 R00, eligible_cancer rows, s-MDK-NCL: n_active must equal the B8 hurdle headline_h1.runs.R00 n_active, and the
whole-column score total must equal run.json G5 row_sum_total. Failure writes gate_failure.json, no result file, no stamp, exit 2.

Inputs come from the CLI (set by ``.commands/submit-v2-celltype-map.sh`` through ``chain_preflight.py``): the B6 COMMOT run (``--commot-run``), the B8 run
(``--hurdle-run``) and ``--out-root``.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE))
import chain_guard as PF  # noqa: E402
import promoted_common as CM  # noqa: E402
import protocol_io as PIO  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/celltype_map.py"
STAGE = "B18"
ARTIFACT = "celltype_map.json"
HEADLINE = "HCC22-088-P4-S2_1i_rep"
PAIRS = ["NCL", "LRP1", "SDC2", "SDC4", "ITGA6_ITGB1"]
MIN_ACTIVE = 30


def split_type(ids: pd.Series) -> pd.Series:
    return ids.str.rsplit("|", n=1).str[1]


def load(run_root, spec: str):
    d = Path(run_root) / spec
    rows = pd.read_csv(d / "rows_R00.csv.gz")
    snd = pd.read_csv(d / "R00" / "sender.csv.gz", index_col=0)
    if len(rows) != len(snd) or not (rows["id"].to_numpy() == snd.index.to_numpy()).all():
        raise CM.GateFailure(f"{spec}: rows_R00 ids do not align with sender.csv.gz")
    return rows, snd, json.loads((d / "R00" / "run.json").read_text())


def gate_p4s2(n_active: int, score_total: float, expected_n: int, g5_total: float) -> list:
    return [CM.gate_equal("P4-S2 R00 eligible cancer s-MDK-NCL n_active", n_active, expected_n, rtol=0, atol=0),
            CM.gate_equal("P4-S2 R00 s-MDK-NCL column total vs run.json G5 row_sum_total", score_total, g5_total, rtol=1e-9, atol=1e-9)]


def stats(score: np.ndarray, types: np.ndarray, tvals) -> list:
    tot = float(score.sum())
    out = []
    for t in tvals:
        v = score[types == t]
        a = v[v > 0]
        out.append({"cell_type": t, "n_compartments": int(len(v)), "n_active": int(len(a)),
                    "frac_active": float(len(a) / len(v)) if len(v) else None,
                    "median_active": float(np.median(a)) if len(a) else None,
                    "mean_active": float(a.mean()) if len(a) else None,
                    "sum": float(v.sum()), "share_of_total": float(v.sum() / tot) if tot > 0 else None})
    return out


def process(run_root, spec: str) -> dict:
    rows, snd, run = load(run_root, spec)
    types = split_type(rows["id"]).to_numpy()
    tvals = sorted(set(types))
    g5 = run.get("G5", {})
    res = {"g5_status": g5.get("status"), "n_rows": int(len(rows)), "pairs": {}}
    for p in PAIRS:
        gp = g5.get("pairs", {}).get(f"s-MDK-{p}", {})
        row = {"g5_pair_ok": gp.get("ok"), "g5_row_sum_total": gp.get("row_sum_total")}
        for side in ("s", "r"):
            sc = snd[f"{side}-MDK-{p}"].to_numpy(float)
            row[side] = stats(sc, types, tvals)
            row[f"{side}_total"] = float(sc.sum())
        # flagged: G5 failed for the s- pair (G5 lists s- pairs only), or a score column is all zero
        row["flagged"] = bool(gp.get("ok") is False or row["s_total"] == 0 or row["r_total"] == 0)
        res["pairs"][p] = row
    return res


NOTE = ("Descriptive only. No p-values, no ranking claims, no cross-type tests. Compartments of different cell types in the same "
        "spot share counts and are not independent. Active = score > 0. Compartments = R00 rows (marker-share proportion >= 0.05). "
        "Cell type = text after the last '|' in the compartment id. Pairs flagged by run.json G5 (or with an all-zero score column) are "
        "marked flagged=true and must not be read as a measured score; G5 lists s- pairs only, the r- column of a flagged pair is flagged with it. "
        "Cohort medians include flagged specimens for the unflagged pairs; MDK-NCL is not flagged in any specimen. P4-S2 reported separately.")


def cohort_table(res: dict) -> list:
    """Cohort table, MDK-NCL (the exploratory check 3 code, unchanged)."""
    types = sorted({c["cell_type"] for r in res.values() for c in r["pairs"]["NCL"]["s"]})
    cohort = []
    for t in types:
        e = {"cell_type": t}
        for side in ("s", "r"):
            sh, act, nsp = [], 0, 0
            for s, r in res.items():
                c = next((c for c in r["pairs"]["NCL"][side] if c["cell_type"] == t), None)
                if c is None:
                    continue
                nsp += 1
                if c["share_of_total"] is not None:
                    sh.append(c["share_of_total"])
                act += int(c["n_active"] >= MIN_ACTIVE)
            p4 = next((c for c in res[HEADLINE]["pairs"]["NCL"][side] if c["cell_type"] == t), None)
            e[side] = {"median_share_across_specimens": float(np.median(sh)) if sh else None, "n_specimens_with_type": nsp,
                       f"n_specimens_ge{MIN_ACTIVE}_active": act,
                       "P4-S2_share": p4["share_of_total"] if p4 else None, "P4-S2_n_active": p4["n_active"] if p4 else None}
        cohort.append(e)
    return cohort


def run(commot_run, hurdle_run, out_root, specimens, run_id: str = "") -> int:
    """Process every specimen (P4-S2 first), gate P4-S2, write ``celltype_map.json``, the flat CSV and ``_COMPLETE.json``. 0 ok, 2 gate failure."""
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    code = HERE.parents[1]
    git_commit = subprocess.run(["git", "-C", str(code), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    specs = sorted(specimens, key=lambda s: (s != HEADLINE, s))      # P4-S2 first
    if HEADLINE not in specs:
        raise SystemExit("the headline specimen is not in the specimen list: the P4-S2 gate cannot be evaluated")
    hj = json.loads((Path(hurdle_run) / "commot_compartment_hurdle.json").read_text())
    expected_n = hj["headline_h1"]["runs"]["R00"]["magnitude"]["n_active"]
    meta = {"stage": STAGE, "run_id": run_id, "produced_by": PRODUCED_BY, "descriptive_only": True, "seed": None, "note": NOTE, "git_commit": git_commit,
            "commot_run": str(commot_run), "hurdle_reference": str(Path(hurdle_run) / "commot_compartment_hurdle.json"), "pairs": PAIRS, "min_active": MIN_ACTIVE,
            "headline_specimen": HEADLINE,
            "rna_covariate": hj.get("rna_covariate")}      # protocol deviation 2026-11-06: B18 reads no covariate itself; this is the mode of the B8 hurdle it is gated against ('total' for v2b81)
    res = {}
    for spec in specs:
        res[spec] = process(commot_run, spec)
        if spec == HEADLINE:
            rows, snd, run_json = load(commot_run, spec)
            el = rows["eligible_cancer"].to_numpy(bool) & (split_type(rows["id"]).to_numpy() == "Cancer_Cells")
            sc = snd["s-MDK-NCL"].to_numpy(float)
            try:
                gates = gate_p4s2(int((sc[el] > 0).sum()), float(sc.sum()), expected_n, run_json["G5"]["pairs"]["s-MDK-NCL"]["row_sum_total"])
            except CM.GateFailure as e:
                PIO.atomic_write_json(out / "gate_failure.json", {**meta, "failure": str(e), "specimen": spec})
                print(e, file=sys.stderr)
                return 2
            res[spec]["gate"] = gates
            print("gate ok", gates, flush=True)
    g5 = {s: {"status": r["g5_status"], "failed_pairs": [f"s-MDK-{p}" for p, v in r["pairs"].items() if v["g5_pair_ok"] is False],
              "pair_ok": {f"s-MDK-{p}": v["g5_pair_ok"] for p, v in r["pairs"].items()}} for s, r in res.items()}
    flat = []
    for s, r in res.items():
        for p, v in r["pairs"].items():
            for side in ("s", "r"):
                for c in v[side]:
                    flat.append({"specimen": s, "pair": f"MDK-{p}", "side": side, "flagged_G5_or_zero": v["flagged"], "g5_status": r["g5_status"], **c})
    pd.DataFrame(flat).to_csv(out / "celltype_map.csv", index=False)
    PIO.atomic_write_json(out / ARTIFACT, {**meta, "specimens": res, "g5_status": g5, "cohort_mdk_ncl": cohort_table(res)})
    PIO.atomic_write_json(out / "_COMPLETE.json", {"path": str(out), "jobs": [os.environ.get("SLURM_JOB_ID") or "no-slurm"], "sha": git_commit,
                                                  "completed": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "produced_by": PRODUCED_BY})
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--commot-run", type=Path, required=True, help="the B6 COMMOT run dir (id from v42_run_ids.yaml)")
    ap.add_argument("--hurdle-run", type=Path, required=True, help="the B8 run dir (id from v42_run_ids.yaml); its P4-S2 n_active is the gate")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<B18 id>; never overwritten")
    ap.add_argument("--run-id", required=True, help="the B18 id")
    ap.add_argument("--protocol", type=Path, default=HERE / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry(STAGE, a.protocol, run_id=a.run_id, out_root=a.out_root, commot_root=a.commot_run, hurdle_run=a.hurdle_run)
    if a.out_root.exists():
        raise SystemExit(f"refusing: {a.out_root} exists (a run dir is never overwritten)")
    stamp = PF.stage_stamp_begin(STAGE, a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id)  # job-start stamp gate, before any compute
    specimens = json.loads(a.protocol.read_text())["specimens"]["hurdle_specimens_all"]
    rc = run(a.commot_run, a.hurdle_run, a.out_root, specimens, a.run_id)
    if rc == 0 and stamp:
        stamp.finish(commot_root=a.commot_run, hurdle_run=a.hurdle_run)
    print("done" if rc == 0 else "gate failed", a.run_id, a.out_root, flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
