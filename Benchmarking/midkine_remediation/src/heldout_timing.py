"""Usage: heldout_timing.py JOB_ID -> synthetic/g6_fix/prod_coverage_heldout/heldout_timing.json.

Records the SLURM Start of every held-out task (sacct), the criterion commit and its committer time (read from each output json), and
FAILS unless every task Start is strictly later than the criterion commit time (prospective registration proven against SLURM, not commit order).
"""
import datetime
import json
import subprocess
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "synthetic" / "g6_fix" / "prod_coverage_heldout"
NY = datetime.timezone(datetime.timedelta(hours=-4))        # CRC logs in local time (EDT on the run date); commit time carries its own offset


def parse_start(s: str) -> datetime.datetime:
    return datetime.datetime.fromisoformat(s).replace(tzinfo=NY)


def build(job_id: str, sacct_lines: list, commit: str, commit_time: str) -> dict:
    starts = {}
    for ln in sacct_lines:
        jid, start = ln.split("|")[:2]
        if jid.startswith(f"{job_id}_") and start not in ("Unknown", "None", ""):
            starts[jid] = start
    if len(starts) != 4:
        raise SystemExit(f"expected 4 started tasks for {job_id}, got {sorted(starts)}")
    ct = datetime.datetime.fromisoformat(commit_time)
    first = min(parse_start(s) for s in starts.values())
    ok = all(parse_start(s) > ct for s in starts.values())
    return {"job_id": job_id, "criterion_commit": commit, "criterion_commit_time": commit_time, "task_starts": starts,
            "job_start_min": first.isoformat(), "every_task_start_after_commit": ok, "seconds_between_commit_and_first_start": (first - ct).total_seconds()}


def main():
    job = sys.argv[1]
    outs = sorted(OUT.glob("prod_*_s6.json"))
    heldout = [json.loads(p.read_text())["heldout"] for p in outs]
    commits = {(h["criterion_commit"], h["criterion_commit_time"]) for h in heldout}
    if len(outs) != 4 or len(commits) != 1:
        raise SystemExit("need 4 outputs sharing one criterion commit")
    sacct = subprocess.run(["sacct", "-M", "htc", "-j", job, "-X", "-n", "-P", "--format=JobID,Start,End"], capture_output=True, text=True, check=True).stdout.splitlines()
    (commit, ctime), = commits
    rec = build(job, sacct, commit, ctime)
    rec["sacct_raw"] = sacct
    (OUT / "heldout_timing.json").write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=2))
    if not rec["every_task_start_after_commit"]:
        raise SystemExit("PROSPECTIVE ORDER VIOLATED")


if __name__ == "__main__":
    main()
