"""A12.11 producer: Module 4 (NMF) program counts from the frozen path.

    python nmf_module4.py --out-root output/nmf_counts_module3_7type/<protocol_sha12>/<RUN_ID> --run-id <RUN_ID>
        [--layers-root output/sace_layers_module3_7type/<protocol_sha12>] [--layers-run-id r1]
        [--module3-dir output/module3_7type]

The ONE frozen path (protocol ``nmf.frozen_module4``): rerun Module 4
``discover_programs_from_layers`` (CITEgeist/model/programs/anchored_program_discovery.py:2217) once per
specimen on the A12.6 SACE layers (``output/sace_layers_module3_7type/<protocol_sha12>/<layers_run_id>/<specimen>/``;
each type's slice ``layer_<Type>.npz`` is rebuilt into ``adata.layers['<Type>_genes_pass1']``) with the finetuned
module3_7type proportions (``output/module3_7type/<specimen>_cell_prop_finetuned_results.csv``). Every other
argument is the code default. Legacy program tables are NEVER read; the irreproducible manual constants
490 / 231 / 47 and any allowlist entries for legacy program tables are removed.

Parameters (code defaults, frozen): K_programs 5, lambda_spatial 0.0, lambda_sparsity 0.01,
min_expression_threshold 0.0, top_n_genes 50, random_state 42. NMF init nndsvda, max_iter 200, l1_ratio 0.5.
Moran: kNN k=8 row-normalised binary weights, 199 permutations, two-sided p. The Moran permutations use the
global numpy RNG (anchored_program_discovery.py:1247), so the producer calls ``np.random.seed(42)`` immediately
before each specimen to make the p-values reproducible.

Specimens: exactly the 12 of ``examples/slurm/sbatch_module3_7type.sh`` (never the non-representative
P4-S2 / P5-S2). A program is coherent iff Moran's I > 0.15 AND p < 0.01. Counts are reported for all eligible
types and for ``Cancer_Cells`` only. The cancer selection uses the exact name ``Cancer_Cells`` and the producer
ASSERTS a non-empty selection (never a warning).

Output (namespaced, empty-target-refused, atomic, frozen-protocol-stamped):
    <out-root>/nmf_counts.json   merged per-specimen typed statuses + all-type / cancer-only counts

Typed statuses (specimen ``status`` in {ok, not_estimable, failed}; a specimen ``ok`` with a NaN-derived value
aborts via PIO.validate_typed):
    layers_missing, proportions_missing, layer_type_mismatch, type_skipped_lt20_active, cancer_cells_absent,
    nmf_failed_random_fallback, moran_not_computed, nonfinite_program.
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(HERE.parent))
import chain_guard as PF  # noqa: E402
import protocol_io as PIO  # noqa: E402
import provenance  # noqa: E402
from breast_dict_loader import BREAST_PROFILE_V40, types_for_protocol  # noqa: E402

from CITEgeist.model.programs.anchored_program_discovery import discover_programs_from_layers  # noqa: E402

MODULE4_LOGGER = "CITEgeist.model.programs.anchored_program_discovery"

# Exactly the 12 specimens of examples/slurm/sbatch_module3_7type.sh (order matches the sbatch array).
SPECIMENS = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2", "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2", "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep", "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]
TYPES = list(BREAST_PROFILE_V40)  # v40 default ONLY; callers must pass protocol types (types_for_protocol) to run_specimen

# Frozen code-default settings (anchored_program_discovery.py, discover_programs_from_layers signature).
K_PROGRAMS = 5
LAMBDA_SPATIAL = 0.0
LAMBDA_SPARSITY = 0.01
MIN_EXPRESSION_THRESHOLD = 0.0
TOP_N_GENES = 50
RANDOM_STATE = 42
LAYER_PATTERN = "_genes_pass1"

# Frozen coherence rule (protocol nmf.frozen_module4.coherence_rule).
MORAN_I_MIN = 0.15
MORAN_P_MAX = 0.01

# The A12.6 producer wrote its layers under this sub-namespace.
DEFAULT_LAYERS_RUN_ID = "r1"


class _NmfFailedHandler(logging.Handler):
    """Detects the 'NMF failed' warning (the code then substitutes RANDOM matrices, line :2526)."""

    def __init__(self):
        super().__init__(level=logging.WARNING)
        self.failed = False

    def emit(self, record):
        if "NMF failed" in record.getMessage():
            self.failed = True


def _nmf_fallback_during(fn):
    """Run ``fn`` while capturing the Module 4 logger's 'NMF failed' signal; return (result, did_fallback)."""
    logger = logging.getLogger(MODULE4_LOGGER)
    handler = _NmfFailedHandler()
    logger.addHandler(handler)
    try:
        res = fn()
    finally:
        logger.removeHandler(handler)
    return res, handler.failed


def _active_mask(layer: sp.csr_matrix) -> np.ndarray:
    """Spot rows with total expression > min_expression_threshold (0.0); matches the module's :2332 rule."""
    totals = np.asarray(layer.sum(axis=1)).ravel()
    return totals > MIN_EXPRESSION_THRESHOLD


def _program_status(prog, result, active_mask: np.ndarray, k: int) -> dict:
    """Classify the k-th program: coherent / not_coherent / moran_not_computed / nonfinite_program."""
    I = float(prog.spatial_moran_i)
    p = float(prog.spatial_moran_pvalue)
    W = np.asarray(result.W)
    Hk = np.asarray(result.H[k, :])
    h_active = Hk[active_mask]
    finite = bool(all(np.isfinite(x) for x in (I, p)) and np.all(np.isfinite(W[:, k])) and np.all(np.isfinite(h_active)))
    if not finite:
        status = "nonfinite_program"
    else:
        zero_var = bool(h_active.size == 0 or h_active.min() == h_active.max())
        if zero_var and I == 0.0 and p == 1.0:
            status = "moran_not_computed"
        else:
            status = "ok"
    coherent = bool(status == "ok" and I > MORAN_I_MIN and p < MORAN_P_MAX)
    return {
        "program_id": int(prog.program_id),
        "spatial_moran_i": None if not np.isfinite(I) else I,
        "spatial_moran_pvalue": None if not np.isfinite(p) else p,
        "status": status,
        "coherent": coherent,
        "variance_explained": None if not np.isfinite(prog.variance_explained) else float(prog.variance_explained),
        "top_genes": list(prog.top_genes[:8]),
    }


def _not_estimable(specimen, reason_code, detail="", secs=0.0) -> dict:
    return {"specimen": specimen, "status": "not_estimable", "reason_code": reason_code,
            "reason": detail if detail else reason_code, "seconds": round(secs, 1)}


def run_specimen(specimen: str, layers_dir: Path, module3_dir: Path, types=None) -> dict:
    """Run the frozen Module 4 path for one specimen; returns its typed report (aborts on an empty cancer gate)."""
    t0 = time.time()
    TYPES_ = list(TYPES if types is None else types)

    # ---- proportions (the same finetuned file A12.6 used) -------------------------------
    qp_path = module3_dir / f"{specimen}_cell_prop_finetuned_results.csv"
    if not qp_path.exists():
        return _not_estimable(specimen, "proportions_missing", secs=time.time() - t0)

    # ---- layers + conservation gate (read-only) -----------------------------------------
    manifest_path = layers_dir / "manifest.json"
    if not manifest_path.exists():
        return _not_estimable(specimen, "layers_missing", secs=time.time() - t0)
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get("conservation", {}).get("passed", False):
        return _not_estimable(specimen, "layers_missing", "conservation gate failed", time.time() - t0)
    if manifest.get("types") != TYPES_:
        return _not_estimable(specimen, "layer_type_mismatch", f"layer types {manifest.get('types')}", time.time() - t0)

    for _aux in ("genes.txt", "barcodes.txt", "coords.csv"):
        if not (layers_dir / _aux).is_file():
            return _not_estimable(specimen, "layers_missing", f"missing {_aux}", time.time() - t0)
    genes = [ln for ln in (layers_dir / "genes.txt").read_text().splitlines() if ln.strip()]
    barcodes = [ln for ln in (layers_dir / "barcodes.txt").read_text().splitlines() if ln.strip()]
    coords = pd.read_csv(layers_dir / "coords.csv").set_index("barcode").reindex(barcodes)
    xy = coords[["x", "y"]].to_numpy(float)
    n_spots, n_genes = len(barcodes), len(genes)
    if n_spots == 0 or n_genes == 0:
        return _not_estimable(specimen, "layers_missing", "empty genes/barcodes", time.time() - t0)

    prop_df = pd.read_csv(qp_path, index_col=0)
    if list(prop_df.columns) != TYPES_:
        return _not_estimable(specimen, "layer_type_mismatch", f"proportion columns {list(prop_df.columns)}",
                              time.time() - t0)
    prop_df = prop_df.loc[barcodes]

    # ---- Build the AnnData the frozen path consumes (layers '<Type>_genes_pass1', obsm['spatial']).
    import anndata as ad
    adata = ad.AnnData(X=sp.csr_matrix((n_spots, n_genes), dtype=np.float64),
                       obs=pd.DataFrame(index=pd.Index(barcodes, name="barcodes")),
                       var=pd.DataFrame(index=pd.Index(genes, name="genes")))
    adata.obsm["spatial"] = xy
    for t in TYPES_:
        try:
            layer = sp.load_npz(layers_dir / f"layer_{t}.npz").tocsr()
        except (FileNotFoundError, OSError):
            return _not_estimable(specimen, "layer_type_mismatch", f"missing layer_{t}.npz", time.time() - t0)
        except Exception as e:  # noqa: BLE001 - malformed npz (wrong/corrupt encoding)
            return _not_estimable(specimen, "layer_type_mismatch", f"layer_{t}.npz unreadable: {e}",
                                  time.time() - t0)
        if layer.shape != (n_spots, n_genes):
            return _not_estimable(specimen, "layer_type_mismatch", f"{t} shape {layer.shape}", time.time() - t0)
        adata.layers[f"{t}{LAYER_PATTERN}"] = layer

    active_by_type = {t: _active_mask(adata.layers[f"{t}{LAYER_PATTERN}"]) for t in TYPES_}

    # ---- frozen Module 4 call (Moran reads the global RNG -> seed immediately before each specimen).
    np.random.seed(RANDOM_STATE)
    (result, nmf_fallback) = _nmf_fallback_during(
        lambda: discover_programs_from_layers(
            adata,
            cell_type_proportions=prop_df,
            profile_discovery_result=None,
            protein_adata=None,
            K_programs=K_PROGRAMS,
            lambda_spatial=LAMBDA_SPATIAL,
            lambda_sparsity=LAMBDA_SPARSITY,
            min_expression_threshold=MIN_EXPRESSION_THRESHOLD,
            validate_with_proteins=True,  # code default; protein_adata=None so it is a no-op
            top_n_genes=TOP_N_GENES,
            detect_subpopulations=True,   # code default; not used by any count
            n_subpopulations=3,
            random_state=RANDOM_STATE,
        )
    )

    results_by_anchor = result.results_by_anchor
    types = {}
    for t in TYPES_:
        entry = results_by_anchor.get(t)
        n_active = int(active_by_type[t].sum())
        if entry is None:
            # Skipped: n_active < 20 per the module (:2335), or the type never ran. Excluded from every total.
            types[t] = {"status": "type_skipped_lt20_active", "reason": f"n_active={n_active} < 20 (code skips)",
                        "n_active": n_active, "n_programs": 0, "n_coherent": 0, "programs": []}
            continue
        programs = [_program_status(prog, entry, active_by_type[t], k) for k, prog in enumerate(entry.programs)]
        types[t] = {
            "status": "ok",
            "n_active": n_active,
            "n_spots_used": int(entry.n_spots_used),
            "n_programs": int(len(programs)),
            "n_coherent": int(sum(p["coherent"] for p in programs)),
            "nmf_failed_random_fallback": bool(nmf_fallback),
            "programs": programs,
        }

    # ---- Cancer_Cells selection: exact name, ASSERT a non-empty selection (never a warning).
    cancer = results_by_anchor.get("Cancer_Cells")
    cancer_present = bool(cancer is not None and cancer.n_spots_used >= 20 and len(cancer.programs) > 0)
    if not cancer_present:
        raise AssertionError(
            f"Cancer_Cells selection empty for {specimen}: present={cancer is not None}, "
            f"n_spots_used={int(cancer.n_spots_used) if cancer is not None else 0}, "
            f"n_programs={int(len(cancer.programs)) if cancer is not None else 0}. "
            "The frozen path requires a non-empty Cancer_Cells selection (the legacy "
            "Cancer_Basal/Cancer_Luminal/Epithelial silent-empty failure class). Aborting."
        )

    rep = {
        "specimen": specimen, "status": "ok", "reason_code": "ok", "reason": "",
        "n_spots": n_spots, "n_genes": n_genes, "n_types": len(TYPES_),
        "cancer": {"present": True, "n_spots_used": int(cancer.n_spots_used), "n_programs": len(cancer.programs),
                   "n_coherent": int(sum(p["coherent"] for p in types["Cancer_Cells"]["programs"]))},
        "types": types,
        "seconds": round(time.time() - t0, 1),
    }
    PIO.validate_typed({"status": "ok", "value": float(rep["n_spots"])})
    return rep


def counts(reports, type_filter=None) -> dict:
    """Merged counts over the ok specimens. ``type_filter`` restricts to e.g. {'Cancer_Cells'}."""
    total_emitted = coherent = not_coherent = moran_not_computed = nonfinite = discarded = skipped_type = 0
    per_spec = {}
    n_specimens = 0
    for r in reports:
        if r["status"] != "ok":
            continue
        n_specimens += 1
        sc = {"total": 0, "coherent": 0, "not_coherent": 0, "moran_not_computed": 0, "nonfinite": 0,
              "discarded_nmf_fallback": 0, "skipped_type": 0}
        for t, te in r["types"].items():
            if type_filter is not None and t not in type_filter:
                continue
            if te["status"] != "ok":
                skipped_type += 1
                sc["skipped_type"] += 1
                continue
            if te.get("nmf_failed_random_fallback"):
                discarded += te["n_programs"]
                sc["discarded_nmf_fallback"] += te["n_programs"]
                continue
            for p in te["programs"]:
                total_emitted += 1
                sc["total"] += 1
                if p["status"] == "moran_not_computed":
                    moran_not_computed += 1
                    sc["moran_not_computed"] += 1
                elif p["status"] == "nonfinite_program":
                    nonfinite += 1
                    sc["nonfinite"] += 1
                elif p["coherent"]:
                    coherent += 1
                    sc["coherent"] += 1
                else:
                    not_coherent += 1
                    sc["not_coherent"] += 1
        per_spec[r["specimen"]] = sc
    percent = (coherent / total_emitted * 100) if total_emitted else 0.0
    return {
        "n_specimens": n_specimens,
        "total_programs": total_emitted,
        "coherent": coherent,
        "not_coherent": not_coherent,
        "moran_not_computed": moran_not_computed,
        "nonfinite": nonfinite,
        "discarded_nmf_fallback": discarded,
        "skipped_type": skipped_type,
        "percent_coherent": round(percent, 2),
        "per_specimen": per_spec,
        "scope": "Cancer_Cells only" if type_filter == {"Cancer_Cells"} else "all eligible types",
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layers-root", type=Path, default=REPO_ROOT / "output" / "sace_layers_module3_7type")
    ap.add_argument("--layers-run-id", default=None, help=f"the B4 id; REQUIRED for a v42 protocol, legacy default {DEFAULT_LAYERS_RUN_ID}")
    ap.add_argument("--module3-dir", type=Path, default=None, help="the stamped B1 dir; REQUIRED for a v42 protocol, legacy default output/module3_7type otherwise")
    ap.add_argument("--out-root", type=Path, required=True, help="must end in <protocol_sha12>/<run_id>")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)

    proto = PIO.load_protocol(a.protocol)  # raises on a hash mismatch / non-FROZEN status
    PF.guard_entry("B10", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.module3_dir,
                   layers_run=(Path(a.layers_root) / proto["_protocol_sha12"] / a.layers_run_id) if a.layers_run_id else None)
    a.module3_dir = PF.legacy_or_required(a.module3_dir, a.protocol, "--module3-dir", REPO_ROOT / "output" / "module3_7type")
    a.layers_run_id = PF.legacy_or_required(a.layers_run_id, a.protocol, "--layers-run-id", DEFAULT_LAYERS_RUN_ID)
    if tuple(Path(a.out_root).parts[-2:]) != (proto["_protocol_sha12"], a.run_id):
        raise PIO.OutputTargetError(f"--out-root must end in <{proto['_protocol_sha12']}>/<{a.run_id}>: {a.out_root}")
    target = Path(a.out_root) / "nmf_counts.json"
    if target.exists():
        raise PIO.OutputTargetError(f"refusing non-empty target: {target}")
    PIO.check_output_target(Path(a.out_root), proto["_protocol_sha12"], a.run_id)
    stamp = PF.stage_stamp_begin("B10", a.protocol, target, Path(__file__).resolve(), a.run_id)  # v42: job-start stamp gate, before any compute

    layers_base = Path(a.layers_root) / proto["_protocol_sha12"] / a.layers_run_id
    print(f"protocol={proto['_protocol_sha12']} layers_base={layers_base} specimens={len(SPECIMENS)}", flush=True)

    reports = []
    for i, spec in enumerate(SPECIMENS):
        print(f"[{i + 1}/{len(SPECIMENS)}] {spec} ...", flush=True)
        rep = run_specimen(spec, layers_base / spec, a.module3_dir, types_for_protocol(proto))
        reports.append(rep)
        print(f"    status={rep['status']} reason_code={rep['reason_code']} "
              f"n_spots={rep.get('n_spots', '-')} secs={rep.get('seconds')}", flush=True)

    for r in reports:
        if r["status"] not in ("ok", "not_estimable"):
            raise SystemExit(f"unexpected specimen status {r['status']} for {r['specimen']}")
        if r["status"] != "ok" and not r["reason"]:
            raise PIO.GateFailure(f"{r['specimen']} not_estimable without a reason")
    # Fail closed: every `ok` specimen must carry a finite value and no NaN under any `ok` program.
    for r in reports:
        if r["status"] != "ok":
            continue
        PIO.validate_typed({"status": "ok", "value": float(r["n_spots"])})
        for te in r["types"].values():
            if te["status"] != "ok":
                continue
            for p in te["programs"]:
                if p["status"] != "ok":
                    continue
                if not np.isfinite(p["spatial_moran_i"]) or not np.isfinite(p["spatial_moran_pvalue"]):
                    raise PIO.GateFailure(f"{r['specimen']} {te['programs']} ok with non-finite Moran value")

    all_counts = counts(reports)
    cancer_counts = counts(reports, type_filter={"Cancer_Cells"})

    payload = {
        "protocol_sha256": proto["_protocol_sha256"],
        "protocol_sha12": proto["_protocol_sha12"],
        "run_id": a.run_id,
        "layers_root": str(layers_base),
        "module3_dir": str(Path(a.module3_dir)),
        "frozen_settings": {
            "entry_point": "discover_programs_from_layers (anchored_program_discovery.py:2217)",
            "K_programs": K_PROGRAMS, "lambda_spatial": LAMBDA_SPATIAL, "lambda_sparsity": LAMBDA_SPARSITY,
            "min_expression_threshold": MIN_EXPRESSION_THRESHOLD, "top_n_genes": TOP_N_GENES,
            "random_state": RANDOM_STATE, "layer_pattern": LAYER_PATTERN,
            "nmf": "sklearn NMF init=nndsvda, max_iter 200, alpha_W=lambda_sparsity, l1_ratio 0.5, random_state 42",
            "moran": "kNN k=8 row-normalised binary weights, 199 permutations, two-sided p",
            "permutation_seed": "np.random.seed(42) immediately before each specimen (global RNG)",
            "coherence_rule": {"moran_I_min": MORAN_I_MIN, "p_max": MORAN_P_MAX},
        },
        "specimens": reports,
        "counts": {"all_types": all_counts, "cancer_only": cancer_counts},
        "git_head": provenance.code_commit(REPO_ROOT),
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    PIO.atomic_write_json(target, payload)
    if stamp:
        stamp.finish(module3_dir=a.module3_dir, layers_run=layers_base)
    print(f"wrote {target}", flush=True)
    print(f"ALL: {json.dumps(all_counts)}", flush=True)
    print(f"CANCER: {json.dumps(cancer_counts)}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
