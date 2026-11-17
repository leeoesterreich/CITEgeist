"""Run one fix-round G6 run -> synthetic/g6_fix/<name>.json. Usage: run_g6fix.py --variant V (--case C | --inj) [--n N] [--procs P]"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
import provenance  # noqa: E402
import synthetic_g6fix as F  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True, choices=sorted(F.VARIANTS))
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--case", choices=F.CASES)
    g.add_argument("--inj", action="store_true")
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--out", type=Path, default=HERE.parent / "synthetic" / "g6_fix")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if a.inj:
        res = F.injection(a.variant, n=a.n or 50, procs=a.procs)
        name = f"inj_{a.variant}"
    else:
        res = F.type1(a.case, a.variant, n=a.n or 1000, procs=a.procs)
        name = f"t1_{a.variant}_{a.case}"
    provenance.write_stamped_json(a.out / f"{name}.json", {"run": name, "result": res, "rule": "g6_fix_round_selection_rule.md"}, "Benchmarking/midkine_remediation/src/run_g6fix.py")
    print(name, res, flush=True)


if __name__ == "__main__":
    main()
