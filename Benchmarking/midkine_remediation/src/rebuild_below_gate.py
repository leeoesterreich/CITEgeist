"""A12.5b bypass runner: build v4 Hungarian per-cell data for the specimens that are
hard-gated OUT of the normal A12.5 rebuild because their fullres H&E is coarser than
the <=1.0 um/px StarDist segmentation gate.

Per the 2026-10-18 user decision (overriding the 2026-10-14 "gate kept as is"), the
following specimens are built anyway:
    HCC22-088-P3-S2, HCC22-088-P4-S2_1i_rep   (fullres H&E 1.1922 um/px)

This wrapper is a SCOPED, LOCAL-ONLY bypass. It does NOT edit
``rebuild_single_cell.py`` or ``protocol.json``. It shims three gate checks at runtime
for THIS process only, so the normal runner behavior for the 10 supported specimens is
untouched and protocol.json keeps its frozen hash:

  * ``rebuild_single_cell.preflight_specimen``'s ``fullres_he_resolution`` check reads the
    module-level ``MAX_UM_PER_PX``; we raise it above 1.19 so the gate passes here.
  * ``run_specimen`` aborts if the specimen is in
    ``proto["rebuild"]["no_rebuild_below_gate"]["specimens"]``; we wrap
    ``protocol_io.load_protocol`` to return the protocol with that list emptied (the
    hash-verified sha12/sha256 fields are untouched, so the stamped output still carries
    the frozen hash).
  * ``CitegeistModel._validate_he_resolution`` raises on um/px > 1.0 during
    ``run_sace_allocation(single_cell)``; we replace it with a no-op (H&E hard-require).

The rest of ``run_specimen`` (SACE single-cell Hungarian, patchwise-segmentation route,
atomic temp-dir promote, gate1 + rebuild validation, manifest stamping) is byte-identical
to the normal runner.

Usage (mirrors rebuild_single_cell.py):
    python -u rebuild_below_gate.py --specimen <SPEC> --data-root <DATA_ROOT> \
      --module3-dir <stamped B1 dir> --nuclei-dir <stamped B2 dir> --out-root <OUT_ROOT> --run-id <RUN_ID> \
      --protocol Benchmarking/midkine_remediation/protocol.json

    python -u rebuild_below_gate.py --dry-run --specimen <SPEC> ...   (preflight only)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]

# The runner and protocol_io live next to this file; chain_preflight one level up. CITEgeist lives at the repo root.
for p in (str(HERE), str(HERE.parent), str(REPO_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import chain_guard as PF  # noqa: E402
import protocol_io as PIO  # noqa: E402
import rebuild_single_cell as R  # noqa: E402
from CITEgeist.model.citegeist_model import CitegeistModel  # noqa: E402

# UM/px ceiling used for the two gated specimens. The real gate is 1.0; P3-S2/P4-S2_1i_rep
# are 1.1922. 1.25 clears it without relaxing anything else (no "inside" tolerance change).
BYPASS_MAX_UM_PER_PX = 1.25


def _patch_gates() -> None:
    """Defuse the three per-process gates described in the module docstring."""
    # 1. preflight fullres_he_resolution gate (rebuild_single_cell.py:215 reads module global).
    R.MAX_UM_PER_PX = BYPASS_MAX_UM_PER_PX

    # 2. no_rebuild_below_gate list (rebuild_single_cell.py:312). Wrap the hash-verified
    #    loader so the record in the returned dict is emptied; sha12/sha256 untouched.
    _orig_load_protocol = PIO.load_protocol

    def _load_protocol_allow_below_gate(path, allow_draft: bool = False) -> dict:
        proto = _orig_load_protocol(path, allow_draft=allow_draft)
        proto.setdefault("rebuild", {}).setdefault("no_rebuild_below_gate", {}).setdefault("specimens", []).clear()
        return proto

    PIO.load_protocol = _load_protocol_allow_below_gate

    # 3. CitegeistModel._validate_he_resolution hard-raise on um/px > 1.0 (citegeist_model.py:2019).
    CitegeistModel._validate_he_resolution = staticmethod(lambda *a, **k: None)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--specimen", required=True)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--module3-dir", type=Path, required=True, help="the B1 output dir (<specimen>_cell_prop_finetuned_results.csv)")
    ap.add_argument("--out-root", type=Path, required=True)
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--protocol", type=Path, default=R.DEFAULT_PROTOCOL)
    ap.add_argument("--nuclei-dir", type=Path, required=True, help="the stamped B2 artifact dir (<root>/v2b2)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    if a.specimen not in R.SPECIMENS:
        print(f"ERROR: unknown specimen {a.specimen!r}", file=sys.stderr)
        return 2

    PF.guard_entry("B3", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.module3_dir, nuclei_dir=a.nuclei_dir)
    proto = R.load_protocol(a.protocol)
    _patch_gates()  # must be applied before preflight AND before run_specimen

    if a.dry_run:
        res = R.preflight_specimen(a.specimen, a.data_root, a.module3_dir, a.out_root, proto, a.run_id, a.nuclei_dir)
        print(R.format_table([res]))
        for n, c in res["checks"].items():
            if not c["ok"]:
                print(f"  FAIL {a.specimen} {n}: {c['detail']}")
        print(f"overall ok={res['ok']}")
        return 0 if res["ok"] else 1

    R.run_specimen(a.specimen, a.data_root, a.module3_dir, a.out_root, a.run_id, a.protocol, a.nuclei_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
