"""Usage: run_prod_coverage.py --case C [--n N] [--procs P] -> synthetic/g6_fix/prod_coverage/prod_<case>_s6.json (criterion in prod_coverage.py)."""
import argparse
import datetime
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2] / "scripts"))
import prod_coverage as PC  # noqa: E402
import provenance  # noqa: E402

FIX = HERE.parent / "synthetic" / "g6_fix"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=PC.CASES)
    ap.add_argument("--n", type=int, default=PC.N_DATASETS)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--heldout", action="store_true", help="held-out seeds (PC.HELDOUT_SEED0); output defaults to prod_coverage_heldout/; needs A12_CRITERION_COMMIT")
    a = ap.parse_args()
    started = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    if a.out is None:
        a.out = FIX / ("prod_coverage_heldout" if a.heldout else "prod_coverage")
    a.out.mkdir(parents=True, exist_ok=True)
    timing = {}
    if a.heldout:
        commit = os.environ.get("A12_CRITERION_COMMIT", "")
        if not commit:
            raise SystemExit("--heldout requires A12_CRITERION_COMMIT (the commit that fixed seeds + criterion)")
        ctime = subprocess.run(["git", "--no-optional-locks", "show", "-s", "--format=%cI", commit], cwd=HERE, capture_output=True, text=True, check=True).stdout.strip()
        timing = {"criterion_commit": subprocess.run(["git", "--no-optional-locks", "rev-parse", commit], cwd=HERE, capture_output=True, text=True, check=True).stdout.strip(),
                  "criterion_commit_time": ctime, "job_id": os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID")),
                  "task_id": os.environ.get("SLURM_ARRAY_TASK_ID"), "task_started_utc": started}
    new = PC.run_case(a.case, a.n, a.procs, seed0=PC.HELDOUT_SEED0 if a.heldout else PC.SEED0)
    stored = PC.stored_record(a.case, FIX)
    payload = {"run": f"prod_{a.case}_s6", "result": new, "stored_historical": stored, "criterion": PC.reproduces(new, stored),
               "criterion_text": PC.__doc__.split("PRE-STATED", 1)[1].strip()}
    if a.heldout:
        payload["heldout"] = {**timing, "binding_values_compared": True}
    provenance.write_stamped_json(a.out / f"prod_{a.case}_s6.json", payload, "Benchmarking/midkine_remediation/src/run_prod_coverage.py")
    print(payload["run"], new["coverage"], payload["criterion"], flush=True)


if __name__ == "__main__":
    main()
