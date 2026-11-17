"""Run Module 3 (marker-share proportions, v41 7-type breast dictionary) on one specimen.

Thin wrapper around the shipped CITEgeist/examples/run_module3_unified.py. The shipped script carries a
placeholder DATA_ROOT and anonymised sample names, so this wrapper injects the lab data root and the real
specimen names and then calls the shipped ``run_module3`` unchanged (marker share, UNIFIED_PROFILE = v41).

Spot set: this applies the shipped min_counts rule (biopsy 100, surgical 25) and nothing else. The v40
``output/module3_7type`` spot set is a strict subset of the one produced here (P4-S2 3370 vs 3931; P4-S1 510 vs
648; P2-S1 1898 vs 2050), from an undocumented extra filter that no min_counts value reproduces. v41 keeps the
shipped, reproducible spot set.

    python run_module3_marker_share.py --sample HCC22-088-P1-S1 --output-dir output/module3_marker_share_7type
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import argparse
import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
DATA_ROOT = Path(str(_cg.root("DATA") / "GSE289326" / "processed_files"))
SPECIMENS = {  # name -> min_counts (biopsy 100, surgical 25): the shipped run_module3_unified.py rule
    "HCC22-088-P1-S1": 100, "HCC22-088-P1-S2": 25, "HCC22-088-P2-S1": 100, "HCC22-088-P2-S2": 25,
    "HCC22-088-P3-S1_A": 100, "HCC22-088-P3-S2": 25, "HCC22-088-P4-S1": 100, "HCC22-088-P4-S2_1i_rep": 25,
    "HCC22-088-P5-S1": 100, "HCC22-088-P5-S2_F_rep": 25, "HCC22-088-P6-S1": 100, "HCC22-088-P6-S2_D": 25,
}


def load_unified():
    path = REPO / "CITEgeist" / "examples" / "run_module3_unified.py"
    spec = importlib.util.spec_from_file_location("_run_module3_unified", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_run_module3_unified"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, choices=sorted(SPECIMENS))
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--skip-gex", action="store_true")
    a = ap.parse_args()
    m = load_unified()
    m.DATA_ROOT = DATA_ROOT
    m.SAMPLE_METADATA = {k: {"min_counts": v} for k, v in SPECIMENS.items()}
    out = Path(a.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    res = m.run_module3(a.sample, out, profile=m.UNIFIED_PROFILE, skip_gex=a.skip_gex, phase1_method="marker_share")
    res["timestamp"] = datetime.now().isoformat()
    res["profile"] = {ct: list(v.get("Major", [])) for ct, v in m.UNIFIED_PROFILE.items()}
    (out / f"{a.sample}_module3_meta.json").write_text(json.dumps(res, indent=2))
    print(f"done {a.sample}: {res['n_spots']} spots, types={res['cell_types']}")


if __name__ == "__main__":
    main()
