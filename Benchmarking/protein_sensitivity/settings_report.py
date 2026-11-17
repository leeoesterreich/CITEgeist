"""Aggregate the settings-sensitivity sweep (W3) and run leave-one-region-out selection (W3b).

CPU only. Writes <out>/settings_summary.json (stamped), <out>/settings_metrics.csv, <out>/settings_loro.csv.
Gates (all must pass, else exit 1): row count equals the literal expected count, key set equals the
protocol-derived set, and baseline r reproduces the A11 gate within the protocol tolerance.
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from protein_sensitivity import aggregate as A  # noqa: E402
from protein_sensitivity import loro as L  # noqa: E402
from protein_sensitivity import run_sensitivity as R  # noqa: E402
from protein_sensitivity import settings_conditions as SC  # noqa: E402

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
N_TYPES = (6, 7)
EXPECTED_PER_UNIT = 15  # literal: 14 non-baseline settings levels + 1 baseline (spec)
N_REGIONS = 5


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default=str(RESULTS / "settings"))
    ap.add_argument("--gate", default=str(RESULTS / "gate"))
    ap.add_argument("--out", default=str(RESULTS))
    args = ap.parse_args(argv)

    proto = SC.load_protocol()
    tol = proto["baseline"]["baseline_equals_gate_tolerance"]
    df = A.collect(args.results)
    assert len(df) > 0, "no settings rows collected"
    if "param" not in df.columns:
        df["param"] = ""
    df["param"] = df["param"].fillna("")
    gate = A.collect(args.gate)
    problems = []

    # gate 1: literal row count and key set
    for nt in N_TYPES:
        sub = df[(df["dataset"] == "xenium") & (df["n_types"] == nt)]
        want = EXPECTED_PER_UNIT * N_REGIONS
        if len(sub) != want:
            problems.append(f"{nt}-type: expected {want} rows, found {len(sub)}")
        exp_ids = sorted(R.condition_id(c) for c in SC.build_list(nt, proto))
        for u in map(str, proto["regions"]):
            recs = sub[sub["unit"].astype(str) == u].to_dict("records")
            ids = sorted(
                "baseline" if r["axis"] == "baseline" else R.condition_id({"axis": r["axis"], "param": r["param"], "level": r["level"]})
                for r in recs
            )
            if ids != exp_ids:
                problems.append(f"{nt}-type region {u}: condition-id set mismatch: {sorted(set(ids) ^ set(exp_ids))}")

    # gate 2: baseline equals A11 gate before any delta is read
    gate_checks = {}
    for nt in N_TYPES:
        g = gate[(gate["dataset"] == "xenium") & (gate["n_types"] == nt) & (gate["axis"] == "baseline")]
        gate_r = {str(u): float(r) for u, r in zip(g["unit"], g["pearson_r"])}
        gate_checks[nt] = A.check_baseline_gate(df, nt, gate_r, tol)
        if not gate_checks[nt]["pass"]:
            problems.append(f"{nt}-type baseline != gate: {gate_checks[nt]['mismatches']}")

    summaries = {str(nt): A.summarize_settings(df, nt, proto) for nt in N_TYPES}
    grid = L.to_grid(df[df["dataset"] == "xenium"], proto)
    sel = L.select(grid, proto["baseline"]["levels"])
    loro_summary = L.summarize(sel)

    git_hash = R._git_hash()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "settings_metrics.csv", index=False)
    sel.to_csv(out / "settings_loro.csv", index=False)
    summary = {
        "stamp": {
            "git_hash": git_hash,
            "settings_protocol_md5": SC.protocol_md5(),
            "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
            "results_root": str(args.results),
            "n_rows": int(len(df)),
        },
        "baseline_gate": {str(k): v for k, v in gate_checks.items()},
        "problems": problems,
        "sensitivity": summaries,
        "loro": {
            "rule": proto["loro_rule"],
            "selection": json.loads(sel.to_json(orient="records")),
            "summary": json.loads(loro_summary.to_json(orient="records")),
        },
    }
    (out / "settings_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"wrote {len(df)} rows -> {out}")
    if problems:
        print("GATE FAILED:", *problems, sep="\n  ")
        return 1
    print("GATE PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
