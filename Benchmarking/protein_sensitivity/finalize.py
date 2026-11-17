"""Finalizer for the v41p sensitivity run: nine gates, then the DVC completion marker.

``finalize`` always returns all nine named checks; ``_COMPLETE.json`` is written only if every one is True.
Expected counts are fixed constants (``run_check``), never derived from the outputs being checked.

CLI: python -m protein_sensitivity.finalize --root R --jobs J [J ...]      (checks, marker, prints the DVC command)
     python -m protein_sensitivity.finalize --root R --snapshot            (write old_dirs_snapshot.json)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (XB, REPO / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import provenance  # noqa: E402
from gt_rebuild import manifest  # noqa: E402
from protein_sensitivity import gate as G  # noqa: E402
from protein_sensitivity import run_check as RC  # noqa: E402
from protein_sensitivity import strict_score as SS  # noqa: E402

CHECKS = (
    "run_check_clean",
    "gate_ok_6_and_7",
    "single_mask_hash",
    "scores_nonempty",
    "fallback_table_present",
    "stamps_clean",
    "old_dirs_unchanged",
    "acceptance_counts",
    "roster_hashes_match",
)
OLD_DIRS = ("sensitivity_marker_share_v41", "sensitivity_gexfix_v2", "spots_sensitivity")
RESULTS_REL = Path("Benchmarking") / "protein_sensitivity" / "results"
SNAPSHOT = "old_dirs_snapshot.json"
SCORE_CSVS = ("strict_scores", "seed_level", "region_type")
FALLBACK_CSVS = ("fallback_table", "fallback_rates")
DVC_ID = "protein_sensitivity.results.v41p_rnaprot"


# ---------------------------------------------------------------- old-dir snapshot (content hashes only)
def _dir_sha(d: Path) -> str:
    h = hashlib.sha256()
    for p in sorted(x for x in d.rglob("*") if x.is_file()):
        h.update(f"{p.relative_to(d).as_posix()}\t{manifest.sha256_file(p)}\n".encode())
    return h.hexdigest()


def snapshot_old_dirs(repo: Path = REPO) -> dict:
    out = {}
    for name in OLD_DIRS:
        d = Path(repo) / RESULTS_REL / name
        out[name] = {"sha256": _dir_sha(d), "n_files": sum(1 for x in d.rglob("*") if x.is_file())} if d.is_dir() else None
    return out


def write_snapshot(root, repo: Path = REPO) -> dict:
    snap = snapshot_old_dirs(repo)
    if any(v is None for v in snap.values()):
        raise FileNotFoundError(f"old dir missing: {[k for k, v in snap.items() if v is None]}")
    (Path(root) / SNAPSHOT).write_text(json.dumps(snap, indent=2, sort_keys=True) + "\n")
    return snap


def verify_old_dirs(root, repo: Path = REPO) -> List[str]:
    p = Path(root) / SNAPSHOT
    if not p.is_file():
        return [f"missing {SNAPSHOT}"]
    saved, now = json.loads(p.read_text()), snapshot_old_dirs(repo)
    if set(saved) != set(OLD_DIRS):
        return [f"snapshot covers {sorted(saved)}, expected {sorted(OLD_DIRS)}"]
    return [f"{k} changed since the snapshot" for k in OLD_DIRS if saved[k] != now[k] or now[k] is None]


# ---------------------------------------------------------------- the nine checks
def _stamped_files(root: Path, n_types_list) -> List[Path]:
    files = [root / "sensitivity_summary.json", root / "spots_vs_citegeist.json", root / "inputs_record.json"]
    for n in n_types_list:
        files += [root / "gate" / f"baseline_gate_{n}type.json", root / "mask" / f"mask_{n}type.json"]
    files += [root / "scores" / f"{c}.csv" for c in SCORE_CSVS + FALLBACK_CSVS]
    return files


def _stamp_of(p: Path) -> dict:
    return json.loads((provenance._sidecar_path(p) if p.suffix == ".csv" else p).read_text())


def _check_run_clean(root, ntl) -> bool:
    return RC.check_run(root, tuple(ntl)) == []


def _check_gates(root, ntl) -> bool:
    try:
        for n in ntl:
            G.require_gate(root, n)
    except (G.GateFailed, OSError, ValueError, KeyError):
        return False
    return True


def _check_mask(root, ntl) -> bool:
    sl = pd.read_csv(root / "scores" / "seed_level.csv")
    for n in ntl:
        sub = sl[sl["n_types"] == n]["mask_sha256"]
        mf = root / "mask" / f"mask_{n}type.json"
        if sub.empty or sub.nunique() != 1 or not mf.is_file():
            return False
        if sub.iloc[0] != json.loads(mf.read_text())["sha256"]:
            return False
    return True


def _check_scores(root, ntl) -> bool:
    for c in SCORE_CSVS:
        p = root / "scores" / f"{c}.csv"
        if not p.is_file() or len(pd.read_csv(p)) == 0:
            return False
    return True


def _check_fallback(root, ntl) -> bool:
    for c in FALLBACK_CSVS:
        p = root / "scores" / f"{c}.csv"
        if not p.is_file() or len(pd.read_csv(p)) == 0:
            return False
    return set(ntl) <= set(pd.read_csv(root / "scores" / "fallback_table.csv")["n_types"])


def _check_stamps(root, ntl) -> bool:
    commits = set()
    for p in _stamped_files(root, ntl):
        try:
            ok, _ = provenance.verify(p)
            st = _stamp_of(p)
        except (OSError, ValueError):
            return False
        if not ok or st.get("_code_dirty") is not False:
            return False
        commits.add(st.get("_code_commit"))
    return len(commits) == 1 and None not in commits


def _check_counts(root, ntl) -> bool:
    sl = pd.read_csv(root / "scores" / "seed_level.csv", keep_default_na=False, na_values=[""])
    for n in ntl:
        conds = json.loads((root / f"conditions_xenium_{n}.json").read_text())
        if len(conds) != RC.EXPECTED_N_CONDITIONS[n]:
            return False
        for method in RC.METHODS:
            sub = sl[(sl["n_types"] == n) & (sl["method"] == method)]
            if sub["condition_id"].nunique() != RC.EXPECTED_N_CONDITIONS[n] or len(sub) != RC.EXPECTED_N_CONDITIONS[n]:
                return False
            bad = sub[sub["status"] != "ok"]
            if not all((r["axis"], r["panel"]) in RC.ALLOWED_FAILURES and r["status"] == "failed_no_markers" for _, r in bad.iterrows()):
                return False
    return True


def _roster_shas(roster: dict) -> dict:
    """method -> ordered sha256 list. Paths are ignored: scoring runs from a worktree whose repo prefix differs."""
    return {m: [e["sha256"] for e in v] for m, v in roster.items()}


def _check_roster(root, ntl) -> bool:
    rec = json.loads((root / "inputs_record.json").read_text()).get("roster") or {}
    for n in ntl:
        if str(n) not in rec or _roster_shas(rec[str(n)]) != _roster_shas(SS.roster_prediction_hashes(n)):
            return False
    return True


def finalize(
    root,
    jobs: Sequence[str],
    n_types_list: Sequence[int] = (6, 7),
    repo: Path = REPO,
) -> Dict[str, bool]:
    root = Path(root)
    runners = {
        "run_check_clean": lambda: _check_run_clean(root, n_types_list),
        "gate_ok_6_and_7": lambda: _check_gates(root, n_types_list),
        "single_mask_hash": lambda: _check_mask(root, n_types_list),
        "scores_nonempty": lambda: _check_scores(root, n_types_list),
        "fallback_table_present": lambda: _check_fallback(root, n_types_list),
        "stamps_clean": lambda: _check_stamps(root, n_types_list),
        "old_dirs_unchanged": lambda: verify_old_dirs(root, repo) == [],
        "acceptance_counts": lambda: _check_counts(root, n_types_list),
        "roster_hashes_match": lambda: _check_roster(root, n_types_list),
    }
    res: Dict[str, bool] = {}
    for name in CHECKS:
        try:
            res[name] = bool(runners[name]())
        except (OSError, ValueError, KeyError, AssertionError, G.GateFailed):
            res[name] = False
    if all(res.values()):
        sha = _stamp_of(root / "sensitivity_summary.json")["_code_commit"][:8]
        manifest.write_complete_marker(root, [str(j) for j in jobs], sha)
    return res


def dvc_command(root) -> str:
    rel = Path(root).resolve().relative_to(REPO)
    return (
        f".commands/dvc-record.sh {rel} --id {DVC_ID} --version v41p_rnaprot "
        "--producer Benchmarking/protein_sensitivity/aggregate_v41p.py"
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--jobs", nargs="+")
    ap.add_argument("--snapshot", action="store_true")
    ap.add_argument("--n-types", type=int, nargs="+", default=[6, 7])
    a = ap.parse_args(argv)
    if a.snapshot:
        write_snapshot(a.root)
        print(f"wrote {Path(a.root) / SNAPSHOT}")
        return 0
    if not a.jobs:
        ap.error("--jobs is required unless --snapshot")
    res = finalize(a.root, a.jobs, tuple(a.n_types))
    for k, v in res.items():
        print(f"{'PASS' if v else 'FAIL'}  {k}")
    if not all(res.values()):
        print("NO marker written", file=sys.stderr)
        return 1
    print(f"wrote {Path(a.root) / '_COMPLETE.json'}")
    print("next (not run):", dvc_command(a.root))
    return 0


if __name__ == "__main__":
    sys.exit(main())
