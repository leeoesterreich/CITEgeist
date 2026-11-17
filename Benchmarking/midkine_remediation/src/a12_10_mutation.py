"""A12.10 mutation test: each gate turns RED on a perturbed input and GREEN on the real one.

The morans gates are:
  * Definition C estimability (eligible cancer compartments and nonzero-MDK counts);
  * the dense-eigh ceiling (n > 6000 -> not_estimable 'exceeds_eigh');
  * Definition A: the two H&E-resolution-gate specimens are typed not_estimable, and a missing
    v4 h5ad is typed not_estimable;
  * G8 literal ``n_estimable + n_not_estimable == 12`` (fails closed on a dropped specimen).
Pure inference helpers (BH, exact sign-flip, Holm, trajectory, aggregate) are exercised directly.

Run: python -u Benchmarking/midkine_remediation/src/a12_10_mutation.py   (CITEgeist_env)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import a12_10_real as A  # noqa: E402
import compartments as C  # noqa: E402
import protocol_io as PIO  # noqa: E402

SPEC = "HCC22-088-P0-S0"
SEC = A.SECRETORY
GENES = SEC + ["MDK", "KRT5", "KRT6A", "KRT6B", "KRT14", "KRT16", "KRT17"]


def _write_comp(root: Path, n: int, mdk_nz: int):
    """Small synthetic compartment h5ad: n eligible cancer compartments, ``mdk_nz`` nonzero MDK."""
    import anndata as ad

    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    genes = list(GENES)
    ng = len(genes)
    X = np.zeros((n, ng), dtype=np.float32)
    for j, g in enumerate(genes):
        X[:, j] = rng.random(n) + 0.1
    midx = genes.index(C.MDK)
    X[:, midx] = 0.0
    if mdk_nz > 0:
        X[:mdk_nz, midx] = rng.random(mdk_nz) * 4 + 1.0
    side = int(np.ceil(np.sqrt(n)))
    coords = np.array([[i // side, i % side] for i in range(n)], float)
    obs = pd.DataFrame({
        "spot_id": np.arange(n), "barcode": [f"bc{i}" for i in range(n)],
        "cell_type": ["Cancer_Cells"] * n, "qp_prop": 0.5, "total_counts": 500.0,
        "depth_ratio": 1.0,
    }, index=[f"bc{i}|Cancer_Cells" for i in range(n)])
    a = ad.AnnData(X=X.astype(np.float32), obs=obs)
    a.var_names = pd.Index(genes)
    a.obsm["spatial"] = coords
    a.layers["N_P_nf"] = X.astype(np.float32)
    a.write_h5ad(str(root / f"{SPEC}.h5ad"))
    return None


def _write_v4(root: Path, n_cells: int, mdk_nz: int, n_spots: int):
    import anndata as ad

    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(11)
    genes = list(GENES)
    X = np.zeros((n_cells, len(genes)), dtype=np.float32)
    for j in range(len(genes)):
        X[:, j] = rng.random(n_cells) + 0.1
    midx = genes.index(C.MDK)
    X[:, midx] = 0.0
    if mdk_nz > 0:
        X[:mdk_nz, midx] = rng.random(mdk_nz) * 4 + 1.0
    spot_bc = [f"s{i % n_spots}" for i in range(n_cells)]
    coords = np.array([[i % 20, i // 20] for i in range(n_cells)], float)
    obs = pd.DataFrame({"cell_id": list(range(n_cells)), "spot_barcode": spot_bc,
                        "cell_type": ["Cancer_Cells"] * n_cells})
    a = ad.AnnData(X=X.astype(np.float32), obs=obs)
    a.var_names = pd.Index(genes)
    a.obsm["spatial"] = coords
    out = root / "r1" / SPEC / f"{SPEC}_single_cell.h5ad"
    out.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(str(out))
    return None


PASS = []


def check(label, fn, expect_raise=False, expect_status=None, expect_code=None, expect_value=None):
    try:
        out = fn()
    except PIO.GateFailure as e:
        ok = bool(expect_raise)
        PASS.append((label, ok))
        print(f"[{'PASS' if ok else 'FAIL'}] {label}: raised GateFailure: {e}", flush=True)
        if not ok:
            raise SystemExit(f"mutation {label} FAILED (unexpected raise)")
        return
    except Exception as e:  # noqa: BLE001
        ok = bool(expect_raise)
        PASS.append((label, ok))
        print(f"[{'PASS' if ok else 'FAIL'}] {label}: raised {type(e).__name__}: {e}", flush=True)
        if not ok:
            raise SystemExit(f"mutation {label} FAILED (unexpected exception)")
        return
    status = out.get("status") if isinstance(out, dict) else None
    code = out.get("reason_code") if isinstance(out, dict) else None
    if expect_status is not None and status != expect_status:
        raise SystemExit(f"mutation {label} FAILED: status {status} != {expect_status}")
    if expect_code is not None and code != expect_code:
        raise SystemExit(f"mutation {label} FAILED: reason_code {code} != {expect_code}")
    if expect_value is not None:
        got = out if not isinstance(out, dict) else out.get("value", float('nan'))
        if not np.isclose(got, expect_value, atol=1e-9):
            raise SystemExit(f"mutation {label} FAILED: value {got} != {expect_value}")
    ok = True
    PASS.append((label, ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {label}: status={status} code={code}", flush=True)


def main() -> int:
    import tempfile

    tmp = Path(tempfile.mkdtemp(prefix="a12_10"))
    comp = tmp / "comp"
    # --- Definition C gates
    _write_comp(comp, 40, 20)
    check("def C estimable (real)", lambda: A.def_c_result(A.def_c_core(comp, SPEC, {})),
          expect_status="ok", expect_code=None)
    _write_comp(comp, 29, 20)
    check("def C n_eligible=29 -> eligibility_lt30 (perturbed)", lambda: A.def_c_core(comp, SPEC, {}),
          expect_status="not_estimable", expect_code="eligibility_lt30")
    _write_comp(comp, 40, 9)
    check("def C mdk_nonzero=9 -> mdk_nonzero_lt10 (perturbed)", lambda: A.def_c_core(comp, SPEC, {}),
          expect_status="not_estimable", expect_code="mdk_nonzero_lt10")
    _write_comp(comp, 6001, 20)
    check("def C n=6001 -> exceeds_eigh (perturbed)", lambda: A.def_c_core(comp, SPEC, {}),
          expect_status="not_estimable", expect_code="exceeds_eigh")

    # --- Definition A gates
    v4 = tmp / "v4"
    _write_v4(v4, 48, 20, 12)
    check("def A estimable (real)", lambda: A.def_a_specimen(v4, SPEC), expect_status="ok")
    check("def A H&E gate P3-S2 -> he_resolution_gate (perturbed)", lambda: A.def_a_specimen(v4, "HCC22-088-P3-S2"),
          expect_status="not_estimable", expect_code="he_resolution_gate")
    empty_v4 = tmp / "v4empty"
    check("def A missing v4 -> v4_rebuild_missing (perturbed)", lambda: A.def_a_specimen(empty_v4, SPEC),
          expect_status="not_estimable", expect_code="v4_rebuild_missing")

    # --- G8 literal
    def g8_11():
        return A.compute_g8({f"s{i}": {"status": "ok"} for i in range(11)})
    check("G8 with 11 specimens -> GateFailure (perturbed)", g8_11, expect_raise=True)
    def g8_12():
        return A.compute_g8({f"s{i}": {"status": "ok"} for i in range(10)} |
                            {"s10": {"status": "not_estimable"}, "s11": {"status": "not_estimable"}})
    check("G8 with 12 specimens (10 ok + 2 not_estimable) -> ok", lambda: g8_12()["equals_12"],
          expect_value=True)

    # --- pure inference helpers
    check("BH over 10 spec p-values non-vacuous", lambda: len(A._bh([0.001, 0.01, 0.05, 0.5, 0.9, 0.02, 0.03, 0.04, 0.06, 0.2])),
          expect_value=10)
    check("BH smallest p is significant", lambda: A._bh([0.001, 0.9, 0.95, 0.85])[0] < 0.05, expect_value=True)
    check("sign-flip smallest attainable 2/2^6", lambda: A._sign_flip([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])["smallest_attainable"],
          expect_value=2 / 64)
    _sing_flip = A._sign_flip([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
    check("sign-flip all-positive d_p -> p=0.03125 (two-sided minimum)", lambda: _sing_flip["sign_flip_p"],
          expect_value=0.03125)
    check("Holm over 2 keeps order", lambda: A._holm([0.01, 0.05])[0] <= A._holm([0.01, 0.05])[1], expect_value=True)

    # --- trajectory gate (a pair with a not_estimable member is skipped)
    def_c = {"HCC22-088-P1-S1": {"status": "ok", "I": 0.1}, "HCC22-088-P1-S2": {"status": "ok", "I": 0.2},
             "HCC22-088-P2-S1": {"status": "ok", "I": 0.3}, "HCC22-088-P2-S2": {"status": "not_estimable"},
             "HCC22-088-P3-S1_A": {"status": "ok", "I": 0.4}, "HCC22-088-P3-S2": {"status": "not_estimable"}}
    pairs = [["HCC22-088-P1-S1", "HCC22-088-P1-S2"], ["HCC22-088-P2-S1", "HCC22-088-P2-S2"],
             ["HCC22-088-P3-S1_A", "HCC22-088-P3-S2"]]
    tr = A.trajectory(def_c, pairs)
    check("trajectory skips pairs with a not_estimable member", lambda: tr["n_pairs"], expect_value=1)

    fails = sum(1 for _, ok in PASS if not ok)
    print(f"\nmutation summary: {len(PASS)} checks, {fails} failures", flush=True)
    print("ALL CHECKS PASSED" if fails == 0 else f"{fails} FAILURES", flush=True)
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
