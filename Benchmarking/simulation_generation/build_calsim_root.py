"""Create the calibrated-sim root: a mirror of the canonical replicates tree where every file is a symlink to the
canonical file, EXCEPT the per-unit CITE h5ad (h5ad_objects/Wu_rep_{i}_CITE.h5ad), which is left absent for the
generation array tasks to write (ruling P14). Refuses if NEW_ROOT exists.

Usage: python build_calsim_root.py NEW_ROOT [CANONICAL_ROOT]
"""
import os
import re
import sys
from pathlib import Path

CANONICAL_ROOT = Path(__file__).resolve().parents[1] / "simulation_benchmarking" / "replicates"
CITE_RE = re.compile(r"^Wu_rep_\d+_CITE\.h5ad$")


def build(new_root: Path, canonical: Path = CANONICAL_ROOT) -> dict:
    new_root, canonical = Path(new_root), Path(canonical)
    if new_root.exists() or new_root.is_symlink():
        raise FileExistsError(f"{new_root} exists; refusing")
    if not canonical.is_dir():
        raise FileNotFoundError(canonical)
    n_links = n_skipped = 0
    for dirpath, dirnames, filenames in os.walk(canonical):
        rel = Path(dirpath).relative_to(canonical)
        (new_root / rel).mkdir(parents=True, exist_ok=False if str(rel) == "." else True)
        for f in filenames:
            if f.startswith("."):  # NFS silly-rename (.nfs*) and other dotfiles are never part of the tree
                continue
            if Path(dirpath).name == "h5ad_objects" and CITE_RE.match(f):
                n_skipped += 1
                continue
            (new_root / rel / f).symlink_to(Path(dirpath) / f)
            n_links += 1
    return {"symlinks": n_links, "cite_left_for_generation": n_skipped}


if __name__ == "__main__":
    print(build(Path(sys.argv[1]), Path(sys.argv[2]) if len(sys.argv) > 2 else CANONICAL_ROOT))
