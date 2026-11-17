"""A12.11 mutation check: the Cancer_Cells non-empty assertion is NON-VACUOUS.

Builds real A12.6-format layer dirs on disk, then drives the real ``run_specimen`` (which calls the real
``discover_programs_from_layers``). Verifies:

  * empty Cancer_Cells (n_active < 20 -> the module skips it)  -> AssertionError (the gate fires)
  * non-empty Cancer_Cells (40 active spots)                  -> status ok, cancer present
  * missing layers dir                                        -> not_estimable / layers_missing
  * missing proportion csv                                    -> not_estimable / proportions_missing
  * missing a layer npz                                       -> not_estimable / layer_type_mismatch
  * types ("Cancer_Basal", "Cancer_Luminal") that do NOT exist -> the exact-name select still fires (cancer absent)

Run via a12_11_nmf_mutation.sbatch (htc). Exit 0 iff every check behaves as expected.
"""
from __future__ import annotations

import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nmf_module4 as M  # noqa: E402

RNG = np.random.default_rng(7)
TYPES = M.TYPES
N_SPOTS = 40
N_GENES = 60


def _write_layer_dir(root: Path, type_names, cancer_active: bool, omit_manifest=False, drop_type=None) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    genes = [f"g{i}" for i in range(N_GENES)]
    barcodes = [f"s{i}" for i in range(N_SPOTS)]
    xy = RNG.uniform(0, 100, (N_SPOTS, 2))
    for t in type_names:
        if t == drop_type:
            continue
        if t == "Cancer_Cells" and not cancer_active:
            X = sp.csr_matrix((N_SPOTS, N_GENES), dtype=np.float64)
        else:
            X = sp.random(N_SPOTS, N_GENES, density=0.5, format="csr", random_state=int(RNG.integers(1e6)))
            X.data = np.abs(X.data) + 1.0
        sp.save_npz(root / f"layer_{t}.npz", X.astype(np.float64))
    (root / "genes.txt").write_text("\n".join(genes) + "\n")
    (root / "barcodes.txt").write_text("\n".join(barcodes) + "\n")
    pd.DataFrame({"barcode": barcodes, "x": xy[:, 0], "y": xy[:, 1]}).to_csv(root / "coords.csv", index=False)
    if not omit_manifest:
        (root / "manifest.json").write_text(json.dumps({
            "specimen": "X", "run_id": "r1", "conservation": {"passed": True, "max_rel_error": 0.0},
            "types": list(type_names), "n_spots": N_SPOTS, "n_genes": N_GENES,
        }))
    return root


def _props_dir(base: Path, name: str, present: bool) -> Path:
    d = base / f"m3_{name}"
    d.mkdir(parents=True, exist_ok=True)
    if present:
        pd.DataFrame({t: RNG.uniform(0.05, 0.95, N_SPOTS) for t in TYPES},
                     index=[f"s{i}" for i in range(N_SPOTS)]) \
            .to_csv(d / f"{name}_cell_prop_finetuned_results.csv")
    return d


def check(label, fn, expect_raise=False, expect_status=None, expect_code=None) -> bool:
    try:
        out = fn()
    except AssertionError as e:
        if expect_raise:
            print(f"PASS {label}: assert fired -> {str(e)[:110]}", flush=True)
            return True
        print(f"FAIL {label}: unexpected AssertionError", flush=True)
        return False
    except Exception as e:  # noqa: BLE001
        print(f"FAIL {label}: unexpected {type(e).__name__}: {e}", flush=True)
        return False
    if expect_raise:
        print(f"FAIL {label}: expected AssertionError but got status={out.get('status')}", flush=True)
        return False
    if expect_status is not None and out.get("status") != expect_status:
        print(f"FAIL {label}: status {out.get('status')} != {expect_status}", flush=True)
        return False
    if expect_code is not None and out.get("reason_code") != expect_code:
        print(f"FAIL {label}: reason_code {out.get('reason_code')} != {expect_code}", flush=True)
        return False
    print(f"PASS {label}: status={out.get('status')} code={out.get('reason_code')}", flush=True)
    return True


def main() -> int:
    ok = True
    with tempfile.TemporaryDirectory() as td:
        base = Path(td)

        # (1) empty cancer -> the non-empty assertion must fire (non-vacuous gate)
        d = _write_layer_dir(base / "ld_empty", TYPES, cancer_active=False)
        m3 = _props_dir(base, "S_empty", present=True)
        ok &= check("empty-cancer fires", lambda: M.run_specimen("S_empty", d, m3), expect_raise=True)

        # (2) non-empty cancer -> ok, cancer present
        d = _write_layer_dir(base / "ld_full", TYPES, cancer_active=True)
        m3 = _props_dir(base, "S_full", present=True)
        r = M.run_specimen("S_full", d, m3)
        ok &= check("non-empty-cancer ok", lambda: r, expect_status="ok")
        ok &= bool(r["cancer"]["present"]) and r["cancer"]["n_spots_used"] == N_SPOTS
        if not ok:
            print("FAIL cancer present 40", flush=True)

        # (3) missing manifest -> layers_missing
        d = _write_layer_dir(base / "ld_nomani", TYPES, cancer_active=True, omit_manifest=True)
        m3 = _props_dir(base, "S_nomani", present=True)
        ok &= check("layers_missing", lambda: M.run_specimen("S_nomani", d, m3),
                    expect_status="not_estimable", expect_code="layers_missing")

        # (4) missing proportions csv -> proportions_missing
        d = _write_layer_dir(base / "ld_noprops", TYPES, cancer_active=True)
        m3 = _props_dir(base, "S_noprops", present=False)
        ok &= check("proportions_missing", lambda: M.run_specimen("S_noprops", d, m3),
                    expect_status="not_estimable", expect_code="proportions_missing")

        # (5) missing one layer npz -> layer_type_mismatch
        d = _write_layer_dir(base / "ld_drop", TYPES, cancer_active=True, drop_type="Fibroblasts")
        m3 = _props_dir(base, "S_drop", present=True)
        ok &= check("layer_type_mismatch", lambda: M.run_specimen("S_drop", d, m3),
                    expect_status="not_estimable", expect_code="layer_type_mismatch")

        # (6) legacy/Luminal/Basal/etc. type names -> rejected at the frozen types gate (exact Cancer_Cells only)
        legacy = [t for t in TYPES if t != "Cancer_Cells"]
        d = _write_layer_dir(base / "ld_legacy", ["Cancer_Basal", "Cancer_Luminal", "Epithelial"] + legacy,
                             cancer_active=False)
        m3 = _props_dir(base, "S_legacy", present=True)
        ok &= check("legacy cancer names rejected", lambda: M.run_specimen("S_legacy", d, m3),
                    expect_status="not_estimable", expect_code="layer_type_mismatch")

    print("MUTATION RESULT:", "ALL PASS" if ok else "FAILURES", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
