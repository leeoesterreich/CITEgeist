"""Module 1 (marker interest) retained-marker count for one clinical Visium specimen folder (v42 recount).

Same loader and parameters as examples/scripts/run_module12_discovery.py (morans_k 8, 99 permutations, seed 42).
One SLURM array task = one specimen folder (index into the sorted folder list, 14 folders: 12 analysis specimens plus
the two superseded originals P4-S2 and P5-S2, whose replicates P4-S2_1i_rep and P5-S2_F_rep replaced them).
Writes <out_dir>/<specimen>.json with the retained markers. The manuscript counts only the 12 analysis specimens
(the builder drops the superseded folders by an explicit list, see the manuscript data builder's marker_inventory.json).

Usage (via run_module1_recount.sbatch):  run_module1_recount.py OUT_DIR   (env SLURM_ARRAY_TASK_ID = 0..13)
"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
import _paths as _cg  # noqa: E402

import json
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "examples" / "scripts"))
sys.path.insert(0, str(REPO))
import run_module12_discovery as R  # noqa: E402

R.DATA_ROOT = Path(
    os.environ.get(
        "M1_DATA_ROOT", str(_cg.root("DATA") / "GSE289326" / "processed_files")
    )
)
names = sorted(p.name for p in R.DATA_ROOT.iterdir() if (p / "outs").is_dir())
assert len(names) == 14, names
specimen = names[int(os.environ["SLURM_ARRAY_TASK_ID"])]
X, coords, markers, n = R.load_patient_antibody_data(specimen)
m = R.identify_interesting_markers(
    X=X, coords=coords, marker_names=markers, morans_k=8, morans_n_perm=99, seed=42, verbose=True
)
out = dict(
    sample=specimen,
    n_spots=int(n),
    n_markers=len(markers),
    n_retained=len(m.interesting_markers),
    retained=list(m.interesting_markers),
    commit=subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
)
Path(sys.argv[1], specimen + ".json").write_text(json.dumps(out, indent=1))
print("DONE", specimen, out["n_retained"], "of", out["n_markers"])
