"""Re-run the A12.5 post-run validation (amended 2026-10-22) on a kept .FAILED_* rebuild dir and optionally promote it.

Only the validation rules changed (see results/ec0cf5b80715/rebuild_gate_amendment.md); the allocation itself is not
recomputed. The exact per-spot allocation gate result stored in gate1_allocation.json must still be passed.

Default is read-only (prints the new verdict). --promote writes the new rebuild_validation.json (old one kept as
rebuild_validation.original.json), updates manifest.json, and renames the dir to <out-root>/<specimen> (refuses an
existing target).

Usage: revalidate_failed.py --failed-dir <out-root>/.FAILED_<specimen>_<pid> [--promote]
"""
import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rebuild_single_cell as R  # noqa: E402

sys.path.insert(0, str(R.REPO_ROOT / "scripts"))
import provenance  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--failed-dir", type=Path, required=True)
    ap.add_argument("--nuclei-dir", type=Path, required=True, help="the stamped B2 artifact dir (<root>/v2b2)")
    ap.add_argument("--promote", action="store_true")
    a = ap.parse_args(argv)

    d = a.failed_dir
    if not d.name.startswith(".FAILED_"):
        raise SystemExit(f"not a .FAILED_* dir: {d}")
    old = json.loads((d / "rebuild_validation.json").read_text())
    manifest = json.loads((d / "manifest.json").read_text())
    specimen = manifest["specimen"]
    gate = json.loads((d / "gate1_allocation.json").read_text())
    prop_df = pd.read_csv(d / "qp_proportions_used.csv", index_col=0)
    cmap = pd.read_csv(d / "cell_spot_map.csv")
    b2 = R.NA.load_nuclei(a.nuclei_dir, specimen, prop_df.index, protocol_sha12=manifest["protocol_sha12"])
    val = R.rebuild_validation(prop_df, cmap, R.nuclei_spot_barcodes(b2))
    val.update({"specimen": specimen, "run_id": manifest["run_id"], "protocol_sha256": manifest["protocol_sha256"],
                "gate1_allocation_passed": bool(gate["passed"]), "revalidated_from": d.name})
    val["passed"] = bool(val["passed"] and gate["passed"])
    print(json.dumps({k: v for k, v in val.items() if k not in ("fraction_reference", "fraction_observed")}, indent=1, default=str))
    if not val["passed"]:
        print(f"STILL FAILS: {val['reasons']}")
        return 1
    if not a.promote:
        print("would pass; re-run with --promote to write and rename")
        return 0

    target = d.parent / specimen
    if target.exists():
        raise SystemExit(f"target exists, refusing: {target}")
    (d / "rebuild_validation.json").rename(d / "rebuild_validation.original.json")
    (d / "rebuild_validation.json").write_text(json.dumps(val, indent=2))
    manifest["rebuild_validation_passed"] = True
    manifest["revalidated_git_head"] = provenance.code_commit(R.REPO_ROOT)
    manifest["files"] = {p.name: {"sha256": R._hash_file(p, "sha256"), "bytes": p.stat().st_size}
                         for p in sorted(d.iterdir()) if p.is_file() and p.name != "manifest.json"}
    (d / "manifest.json").write_text(json.dumps(manifest, indent=2))
    d.rename(target)
    print(f"promoted {d.name} -> {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
