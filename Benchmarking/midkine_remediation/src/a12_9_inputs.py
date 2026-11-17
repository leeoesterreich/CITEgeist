"""A12.9 input wiring (built in A12.4-A12.8, NOT yet run on real outcomes): load A12.7 compartments + A12.8 COMMOT sender scores into the shapes
that ``descriptive.run_h1_specimen`` and ``sweep.run_sweep`` take, and the A01-A07 sweep entry that writes through ``sweep.write_sweep``.

``scores[run]`` = the MDK-NCL sender score aligned to the ELIGIBLE cancer compartments of R00, or a typed failure dict
({'status': 'failed'|'not_estimable', 'reason': ...}) when the run failed, is missing, or lacks rows.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compartments as C  # noqa: E402

PAIR = "s-MDK-NCL"


def load_compartments(h5ad_path: Path):
    """(obs, xy) of one A12.7 file; obs index = '<barcode>|<type>'."""
    import anndata as ad

    a = ad.read_h5ad(str(h5ad_path), backed="r")
    obs, xy = a.obs.copy(), np.asarray(a.obsm["spatial"], float).copy()
    a.file.close()
    return obs, xy


def typed(status: str, reason: str) -> dict:
    return {"status": status, "reason": reason}


def load_score(commot_specimen_dir: Path, run_id: str, ids) -> object:
    """Sender score of one run aligned to ``ids`` (compartment index labels), or a typed failure dict."""
    d = Path(commot_specimen_dir) / run_id
    if not (d / "run.json").is_file():
        return typed("failed", f"run {run_id} has no run.json (not run)")
    rec = json.loads((d / "run.json").read_text())
    if rec.get("status") != "ok":
        return typed("not_estimable" if rec.get("status") == "not_estimable" else "failed", f"{run_id}: {rec.get('reason', '')}")
    tab = pd.read_csv(d / "sender.csv.gz", index_col=0)
    if PAIR not in tab.columns:
        return typed("failed", f"{run_id}: no {PAIR} column")
    missing = pd.Index(ids).difference(tab.index)
    if len(missing):
        return typed("failed", f"{run_id}: {len(missing)} requested compartments are absent from the run")
    v = tab[PAIR].reindex(list(ids)).to_numpy(float)
    if not np.isfinite(v).all():
        return typed("failed", f"{run_id}: non-finite sender scores")
    return v


def sweep_main(a) -> int:  # pragma: no cover - real data
    import protocol_io as PIO
    import sweep as SW

    proto = PIO.load_protocol(a.protocol)
    obs, xy = load_compartments(Path(a.compartments_run) / f"{a.specimen}.h5ad")
    score = load_score(Path(a.commot_root) / a.specimen, "R00", list(obs.index))
    res = SW.run_sweep(a.specimen, obs, xy, score)
    SW.write_sweep(Path(a.out_root) / a.specimen, res, "Benchmarking/midkine_remediation/src/a12_9_inputs.py", a.protocol, a.run_id)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", action="store_true", help="A01-A07 sweep for the headline specimen (A12.9; not run in A12.4-A12.8)")
    ap.add_argument("--specimen")
    ap.add_argument("--compartments-run", type=Path)
    ap.add_argument("--commot-root", type=Path)
    ap.add_argument("--out-root", type=Path)
    ap.add_argument("--run-id")
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)
    if not a.sweep:
        ap.error("only --sweep is wired here; the H1 driver is the A12.9 step")
    return sweep_main(a)


if __name__ == "__main__":
    sys.exit(main())
