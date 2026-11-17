"""Preflight for the v41p protein-sensitivity rerun: refuse before anything is written.

Creates nothing. ``python -m protein_sensitivity.preflight --out-root R [--out-root R2] --code-dir D`` prints JSON
and exits 1 on failure; ``--list`` prints the repo-relative required inputs (used to symlink them into a worktree).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import List, Sequence

REPO = Path(__file__).resolve().parents[2]
XB = REPO / "Benchmarking" / "xenium_benchmarking"
for _p in (XB, XB / "evaluation" / "src", REPO / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import benchmark_paths as bp  # noqa: E402
import provenance  # noqa: E402
from gt_rebuild import manifest  # noqa: E402
from gt_rebuild.score_methods import pred_dirs_for  # noqa: E402
from protein_sensitivity import baseline  # noqa: E402

LABELSET = "rnaprot"


class PreflightError(RuntimeError):
    pass


def _rel(p: Path) -> Path:
    return Path(p).relative_to(REPO)


def required_inputs(n_types_list: Sequence[int] = (6, 7)) -> List[Path]:
    """Repo-relative inputs the gate/sweep read (Xenium data, rnaprot GT, canonical v41p + roster predictions)."""
    out: List[Path] = [
        _rel(baseline.XENIUM_DATA / "h5ad_objects"),
        _rel(baseline.XENIUM_DATA / "cell_to_spot_mapping.csv"),
        _rel(baseline.XENIUM_DATA / "cell_type_assignments.csv"),
        _rel(XB / "ground_truth_xenium_rnaprot"),
        _rel(XB / "ground_truth_xenium_clusters"),  # whole dir: the rnaprot sign-off reads manifest.json + reports/
        _rel(bp.eval_results_dir(LABELSET)),
    ]
    for n in n_types_list:
        out.append(_rel(bp.citegeist_v41p_pred_dir(n, LABELSET)))
        out.append(_rel(bp.citegeist_pred_dir(n, LABELSET)))
        out.extend(_rel(d) for d in pred_dirs_for(LABELSET, n).values())
    seen, uniq = set(), []
    for p in out:
        if p not in seen:
            seen.add(p)
            uniq.append(p)
    return uniq


def _git(code_dir: Path, *args) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "--no-optional-locks", *args], cwd=code_dir, capture_output=True, text=True)


def run_preflight(out_roots: Sequence[Path], code_dir: Path) -> dict:
    """Raises PreflightError (all findings listed) unless every check passes. Creates nothing."""
    code_dir = Path(code_dir)
    # (1) sign-off first: nothing else is looked at if it is missing
    try:
        manifest.require_rnagt_signoff(labelset=LABELSET)
    except Exception as e:  # noqa: BLE001 - any sign-off failure is a refusal
        raise PreflightError(f"rnaprot sign-off check failed: {e}") from e
    errors: List[str] = []
    guard = REPO / ".commands" / "_refuse_protected_outdir.sh"
    for root in out_roots:
        r = (subprocess.run([str(guard), str(root)], capture_output=True, text=True) if guard.exists()
             else subprocess.CompletedProcess([], 0, "", ""))
        if r.returncode != 0:
            errors.append(f"protected out dir {root}: {r.stderr.strip()}")
        root = Path(root)
        if root.exists() and (not root.is_dir() or any(root.iterdir())):
            errors.append(f"out root {root} exists and is not empty")
    for rel in required_inputs():
        if not (code_dir / rel).exists():
            errors.append(f"missing required input: {rel}")
    # Without a git checkout (public export) the clean-tree checks cannot run; they are skipped and the result says so.
    git_checked = provenance.in_git_checkout(code_dir)
    if git_checked:
        st = _git(code_dir, "status", "--porcelain", "-uno")
        if st.returncode != 0 or st.stdout.strip():
            errors.append(f"tracked tree not clean in {code_dir}: {st.stdout.strip() or st.stderr.strip()}")
        if _git(code_dir, "diff", "--cached", "--quiet").returncode != 0:
            errors.append(f"staged changes present in {code_dir}")
    if errors:
        raise PreflightError("\n".join(errors))
    try:
        head, head_source = provenance.code_commit_info(code_dir)
    except provenance.CodeProvenanceError as e:
        raise PreflightError(str(e)) from e
    return {
        "ok": True,
        "head": head,
        "head_source": head_source,
        "clean_tree_checked": git_checked,
        "out_roots": [str(r) for r in out_roots],
        "code_dir": str(code_dir),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-root", action="append", default=[])
    ap.add_argument("--code-dir", default=str(REPO))
    ap.add_argument("--list", action="store_true", help="print required inputs (repo-relative), one per line")
    a = ap.parse_args(argv)
    if a.list:
        print("\n".join(str(p) for p in required_inputs()))
        return 0
    try:
        print(json.dumps(run_preflight([Path(r) for r in a.out_root], Path(a.code_dir)), indent=2))
        return 0
    except PreflightError as e:
        print(f"PREFLIGHT FAILED:\n{e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
