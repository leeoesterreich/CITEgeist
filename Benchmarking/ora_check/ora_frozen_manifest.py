"""SHA-256 manifest of the frozen ORA files (rule: ora_rebuilt_DECISION_RULE.md). The corrected-mapping work reads these and must not edit them.

  python ora_frozen_manifest.py write    writes ora_frozen_sha256.json (refuses to overwrite)
  python ora_frozen_manifest.py check    exits 1 on any missing or changed file
"""
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "ora_frozen_sha256.json"
# patterns relative to this directory
PATTERNS = (
    "ora_protocol.json",
    "ora_protocol_v41_7type.json",
    "term_mapping.yaml",
    "term_mapping_7type.yaml",
    "term_mapping_xenium7.yaml",
    "cellmarker2024_term_lineage.tsv",
    "libraries/*.gmt",
    "ora_v2_DECISION_RULE.md",
    "ora_sense_DECISION_RULE.md",
    "ora_audit_DECISION_RULE.md",
    "ora_lib.py",
    "ora_sense_lib.py",
    "ora_audit_lib.py",
    "ora_v2_lib.py",
    "run_ora_sense.py",
    "run_ora_audit.py",
    "results_v2_20261101_1501/readout.json",
    "results_v2_20261101_1501/readout.md",
    "results_sense_20261101_1550/readout_sense.json",
    "results_sense_20261101_1550/readout_sense.md",
    "results_sense_20261101_1550/calls/*",
    "results_audit_20261101_1624/mapping_audit.tsv",
    "results_audit_20261101_1624/audit.md",
    "results_audit_20261101_1624/repro.json",
    "results_audit_20261101_1624/sens_*.csv",
    "results_audit_20261101_1624/audit/*/calls.csv",
    "results_audit_20261101_1624/audit/*/terms.csv.gz",
)


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def listing(root: Path = HERE) -> list:
    out = []
    for pat in PATTERNS:
        out += [p.relative_to(root).as_posix() for p in sorted(root.glob(pat)) if p.is_file()]
    return sorted(set(out))


def write(root: Path = HERE, path: Path = MANIFEST) -> dict:
    if path.exists():
        raise SystemExit(f"refusing to overwrite {path}")
    files = listing(root)
    if not files:
        raise SystemExit("no frozen files found")
    m = {f: sha(root / f) for f in files}
    path.write_text(json.dumps(m, indent=1, sort_keys=True) + "\n")
    return m


def check(root: Path = HERE, manifest: dict = None) -> list:
    """List of problems (empty = all frozen files byte-identical to the manifest). An empty manifest is itself a problem."""
    manifest = json.loads(MANIFEST.read_text()) if manifest is None else manifest
    if not manifest:
        return ["manifest is empty"]
    bad = []
    for f, h in manifest.items():
        p = root / f
        if not p.is_file():
            bad.append(f"missing: {f}")
        elif sha(p) != h:
            bad.append(f"changed: {f}")
    return bad


if __name__ == "__main__":
    cmd = sys.argv[1:]
    if cmd == ["write"]:
        print(len(write()), "files hashed")
    elif cmd == ["check"]:
        problems = check()
        print("\n".join(problems) or "all frozen files unchanged")
        sys.exit(1 if problems else 0)
    else:
        raise SystemExit(__doc__)
