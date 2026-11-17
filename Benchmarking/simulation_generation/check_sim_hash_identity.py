"""Hash-identity gate: a new sim root must match the canonical sim byte-for-byte
on the spot mapping, spot RNA GEX h5ad and GT proportions for all 10 units
(high_seg/mixed x replicates). A mismatch stops the run (SystemExit).

Usage: python check_sim_hash_identity.py NEW_ROOT [--canonical ROOT]
"""
import argparse
import hashlib
import sys
from pathlib import Path

CANONICAL_ROOT = (
    Path(__file__).resolve().parents[1] / "simulation_benchmarking" / "replicates"
)
CONDITIONS = ("high_seg", "mixed")
# Independent anchor: 2 conditions x 5 replicates (the README's "3" is wrong).
EXPECTED_UNITS = 10


def _artifacts(rep: int) -> dict:
    return {
        "spot_mapping": f"ST_sim/Wu_ST_{rep}_index.csv",
        "spot_rna_gex": f"h5ad_objects/Wu_rep_{rep}_GEX.h5ad",
        "gt_proportions": f"ST_sim/Wu_ST_{rep}_prop.csv",
    }


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_unit(sim_dir: Path, rep: int) -> dict:
    """SHA-256 of the three artifacts for one condition dir and replicate.

    A missing artifact is recorded as "MISSING" so it can never match.
    """
    out = {}
    for name, rel in _artifacts(rep).items():
        p = Path(sim_dir) / rel
        out[name] = _sha256(p) if p.is_file() else "MISSING"
    return out


def assert_identical(new: dict, canonical: dict, unit: str) -> None:
    """Raise SystemExit naming the unit and artifact on any mismatch."""
    for name in sorted(set(new) | set(canonical)):
        a, b = new.get(name, "MISSING"), canonical.get(name, "MISSING")
        if a != b or a == "MISSING":
            raise SystemExit(
                f"HASH IDENTITY FAILED: unit {unit}, artifact {name}: "
                f"new={a} canonical={b}"
            )


def find_replicates(canonical_root: Path, condition: str) -> list:
    reps = sorted(
        int(p.stem.split("_")[2])
        for p in (Path(canonical_root) / condition / "ST_sim").glob("Wu_ST_*_index.csv")
    )
    return reps


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("new_root", type=Path)
    ap.add_argument("--canonical", type=Path, default=CANONICAL_ROOT)
    ap.add_argument("--expect-units", type=int, default=EXPECTED_UNITS)
    args = ap.parse_args(argv)
    plan = []
    for cond in CONDITIONS:
        reps = find_replicates(args.canonical, cond)
        print(f"{cond}: {len(reps)} replicates in canonical root")
        plan += [(cond, r) for r in reps]
    if len(plan) != args.expect_units:
        raise SystemExit(
            f"UNIT COUNT FAILED: canonical root has {len(plan)} units, "
            f"expected {args.expect_units}"
        )
    n = 0
    for cond, rep in plan:
        unit = f"{cond}/rep{rep}"
        assert_identical(
            hash_unit(args.new_root / cond, rep),
            hash_unit(args.canonical / cond, rep),
            unit,
        )
        n += 1
        print(f"OK {unit}")
    print(f"{n} units hash-identical")
    return 0


if __name__ == "__main__":
    sys.exit(main())
