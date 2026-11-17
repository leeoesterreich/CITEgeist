"""Collect every stamped G6/G7 synthetic result into one stamped summary table (read-only on the inputs)."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))

import gate_verdicts as V  # noqa: E402
import provenance  # noqa: E402

SYNTH = HERE.parent / "synthetic"


def main():
    rows = []
    for f in sorted(SYNTH.glob("g6_*.json")):
        if f.stem == "g6_summary":
            continue
        r = V.load_stamped(f)
        name = f.stem[3:]
        if "magnitude" in r:
            rows.append({"run": name, "kind": "type1", "case": r["case"], "scheme": r["scheme"], "mutation": r["mutation"], "n_datasets": r["n_datasets"], "B": r["B"],
                         "magnitude": r["magnitude"]["type1"], "magnitude_ci95": r["magnitude"]["ci95"], "participation": r["participation"]["type1"],
                         "participation_ci95": r["participation"]["ci95"], "n_not_estimable": [r["magnitude"]["n_not_estimable"], r["participation"]["n_not_estimable"]]})
        else:
            rows.append({"run": name, "kind": "injection", **{k: r[k] for k in r if k not in ()}})
    out = SYNTH / "g6_summary.json"
    provenance.write_stamped_json(out, {"rows": rows, "gate": {"max_type1": V.MAX_TYPE1, "required_cases": list(V.G6_NULL_CASES)}}, "Benchmarking/midkine_remediation/src/summarize_synthetic.py")
    print(out, len(rows))


if __name__ == "__main__":
    main()
