"""A12 program-recovery producer: "MDK-associated NMF program recovery by method" (manuscript Figure 5D).

    python a12_program_recovery.py --out-root Benchmarking/midkine_remediation/results/<sha12>/<run_id> --run-id <run_id>
        --c2l-root output/competitor_gex/cell2location --tangram-root output/competitor_gex/tangram

Recomputes the legacy ``aggregate_comparison.classify_programs`` comparison (Benchmarking/competitor_gex_validation/
aggregate_comparison.py; legacy table ``output/competitor_gex/comparison/failure_mode_summary.csv``, analysis
``programs_mdk``) with the CITEgeist side taken from the A12 rebuild: Module 4 (frozen path, protocol
``nmf.frozen_module4``) rerun on the A12.6 SACE layers, SAME parameters and seeding as A12.11 (``nmf_module4.py``),
but keeping the top-50 genes of every program (A12.11's table keeps only 8).

RECOVERY RULE (verbatim from the legacy classifier, applied per specimen to the CANCER programs of a method):
  full    -- some cancer program has MDK in its top-50 genes AND >= 6 of the 11 secretory genes in its top-50
  partial -- MDK is in the top-50 of some cancer program but no program reaches 6 secretory genes
  absent  -- MDK is in no cancer program's top-50 (including: the method has no cancer programs at all)
  Legacy cancer-type names: {"cancer epithelial", "cancer_epithelial", "epithelial"} (competitors, Wu reference);
  CITEgeist: the exact name ``Cancer_Cells``.

COMPETITORS: the per-cell-type GEX layers of Cell2Location / Tangram on the clinical specimens are NOT changed by A12.
Their programs are recomputed here with the legacy NMF settings (sklearn NMF k=5, alpha_W=0.01, init nndsvda,
max_iter 500, random_state 42, mean cell-type proportion >= 0.1 to be eligible, >= 10 non-empty spots) and, when
``--legacy-programs-root`` points at the legacy ``<method>/programs/<specimen>_programs.json`` files, cross-checked
against them (``matches_legacy_top50``). The legacy rule is a competitor-only per-specimen classification; it never
evaluated CITEgeist. The CITEgeist block therefore reports the same rule applied to CITEgeist's own programs plus
the number of MDK-associated programs (MDK in the top-50; ``n_mdk_coherent`` additionally requires the frozen
coherence rule Moran I > 0.15 and p < 0.01).

Output (namespaced, empty-target-refused, atomic, protocol-hash + provenance stamped):
    <out-root>/program_recovery.json
Typed per-specimen statuses (ok | not_estimable | failed); an ``ok`` record always carries a finite value
(``n_cancer_programs``); a missing input is ``not_estimable`` with a reason_code and is NEVER silently counted absent.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(HERE.parent))
import chain_guard as PF  # noqa: E402
import protocol_io as PIO  # noqa: E402
from breast_dict_loader import types_for_protocol  # noqa: E402

PRODUCED_BY = "Benchmarking/midkine_remediation/src/a12_program_recovery.py"
OUT_NAME = "program_recovery.json"

SPECIMENS = [
    "HCC22-088-P1-S1", "HCC22-088-P1-S2", "HCC22-088-P2-S1", "HCC22-088-P2-S2",
    "HCC22-088-P3-S1_A", "HCC22-088-P3-S2", "HCC22-088-P4-S1", "HCC22-088-P4-S2_1i_rep",
    "HCC22-088-P5-S1", "HCC22-088-P5-S2_F_rep", "HCC22-088-P6-S1", "HCC22-088-P6-S2_D",
]
COMPETITORS = ["cell2location", "tangram"]

# Legacy constants (Benchmarking/competitor_gex_validation/constants.py, aggregate_comparison.py).
SECRETORY_GENES = [
    "HSP90B1", "HSPA5", "CALR", "CANX", "PDIA4", "PDIA6",
    "SEC23A", "SEC61B", "ATF6", "MAN1A1", "XBP1",
]
SECRETORY_OVERLAP_THRESHOLD = 6
LEGACY_CANCER_NAMES = frozenset({"cancer epithelial", "cancer_epithelial", "epithelial"})
CITEGEIST_CANCER = "Cancer_Cells"
TOP_N = 50
CATEGORIES = ("full", "partial", "absent")

# Legacy competitor NMF (run_module4_programs.py).
LEGACY_K = 5
LEGACY_ALPHA = 0.01
LEGACY_MAX_ITER = 500
LEGACY_RANDOM_STATE = 42
LEGACY_MIN_MEAN_PROP = 0.1
LEGACY_MIN_SPOTS = 10
C2L_PROP_PREFIX = "q05cell_abundance_w_sf_means_per_cluster_mu_fg_"
C2L_LAYER_PREFIX = "means_per_cluster_mu_fg_"

# Frozen coherence rule (protocol nmf.frozen_module4.coherence_rule), used only for CITEgeist n_mdk_coherent.
MORAN_I_MIN = 0.15
MORAN_P_MAX = 0.01


# --------------------------------------------------------------------------- the recovery rule
def classify_program(top_genes) -> dict:
    """Per-program flags: MDK in top-50, secretory-gene overlap count. Empty gene list is a hard error."""
    genes = list(top_genes)
    if not genes:
        raise PIO.GateFailure("program without top genes: the recovery rule cannot be evaluated")
    top = genes[:TOP_N]
    n_sec = len(set(SECRETORY_GENES) & set(top))
    return {"mdk_in_top50": "MDK" in top, "n_secretory": int(n_sec)}


def classify_specimen(cancer_programs) -> dict:
    """Legacy rule over a specimen's cancer programs: a list of dicts with ``top_genes`` (and optional ``coherent``).

    Returns category full/partial/absent, the count of MDK-associated programs and the best program's detail.
    """
    best, best_idx, best_nsec = "absent", None, None
    n_mdk = n_mdk_coherent = 0
    per_program = []
    for i, p in enumerate(cancer_programs):
        flags = classify_program(p["top_genes"])
        coherent = bool(p.get("coherent", False))
        per_program.append({"program_idx": int(p.get("program_idx", i)), **flags, "coherent": coherent})
        if not flags["mdk_in_top50"]:
            continue
        n_mdk += 1
        n_mdk_coherent += int(coherent)
        if flags["n_secretory"] >= SECRETORY_OVERLAP_THRESHOLD:
            if best != "full":
                best, best_idx, best_nsec = "full", per_program[-1]["program_idx"], flags["n_secretory"]
        elif best != "full":
            best = "partial"
            if best_nsec is None or flags["n_secretory"] > best_nsec:
                best_idx, best_nsec = per_program[-1]["program_idx"], flags["n_secretory"]
    return {
        "category": best,
        "best_program_idx": best_idx,
        "best_n_secretory": best_nsec,
        "n_cancer_programs": len(per_program),
        "n_mdk_programs": n_mdk,
        "n_mdk_coherent": n_mdk_coherent,
        "programs": per_program,
    }


def _ok_record(specimen, method, cls, extra=None) -> dict:
    rec = {"specimen": specimen, "method": method, "status": "ok", "reason_code": "ok", "reason": "",
           "category": cls["category"], "best_program_idx": cls["best_program_idx"],
           "best_n_secretory": cls["best_n_secretory"], "n_cancer_programs": cls["n_cancer_programs"],
           "n_mdk_programs": cls["n_mdk_programs"], "n_mdk_coherent": cls["n_mdk_coherent"],
           "programs": cls["programs"], "value": float(cls["n_cancer_programs"])}
    if cls["n_cancer_programs"] == 0:
        rec["absent_reason"] = "no_cancer_programs"
    elif cls["category"] == "absent":
        rec["absent_reason"] = "mdk_not_in_any_cancer_top50"
    rec.update(extra or {})
    return rec


def _typed_failure(specimen, method, code, detail="", status="not_estimable") -> dict:
    return {"specimen": specimen, "method": method, "status": status, "reason_code": code,
            "reason": detail or code}


def totals(records) -> dict:
    """Merged totals over the ``ok`` records of ONE method; non-ok records are counted, never folded into a category."""
    out = {c: 0 for c in CATEGORIES}
    n_ok = n_not_est = n_failed = 0
    n_mdk = n_mdk_coh = 0
    for r in records:
        if r["status"] == "ok":
            n_ok += 1
            out[r["category"]] += 1
            n_mdk += r["n_mdk_programs"]
            n_mdk_coh += r["n_mdk_coherent"]
        elif r["status"] == "not_estimable":
            n_not_est += 1
        else:
            n_failed += 1
    out.update({"n_specimens_ok": n_ok, "n_not_estimable": n_not_est, "n_failed": n_failed,
                "n_mdk_programs": n_mdk, "n_mdk_coherent_programs": n_mdk_coh})
    return out


# --------------------------------------------------------------------------- competitor side (legacy NMF)
def _strip_prefix(name: str) -> str:
    for pfx in (C2L_PROP_PREFIX, C2L_LAYER_PREFIX):
        if name.startswith(pfx):
            return name[len(pfx):]
    return name


def legacy_nmf_top50(X, gene_names):
    """Legacy ``run_nmf_for_celltype`` minus the Moran step: list of K top-50 gene lists, or None (< 10 active spots)."""
    from sklearn.decomposition import NMF
    X_dense = X.toarray() if hasattr(X, "toarray") else np.array(X, dtype=np.float64)
    X_dense = np.clip(X_dense, 0, None)
    spot_mask = X_dense.sum(axis=1) > 0
    if spot_mask.sum() < LEGACY_MIN_SPOTS:
        return None
    model = NMF(n_components=LEGACY_K, alpha_W=LEGACY_ALPHA, init="nndsvda", max_iter=LEGACY_MAX_ITER,
                random_state=LEGACY_RANDOM_STATE)
    model.fit_transform(X_dense[spot_mask])
    H = model.components_
    return [[gene_names[i] for i in np.argsort(H[k])[::-1][:TOP_N]] for k in range(LEGACY_K)]


def competitor_cancer_programs(layers: dict, gene_names, prop_df: pd.DataFrame) -> list:
    """Cancer-type programs of one competitor specimen under the legacy eligibility + naming rules.

    ``layers``: {layer_name: (spots x genes) matrix} already row-aligned to ``prop_df``.
    """
    prop_clean = {_strip_prefix(c): c for c in prop_df.columns}
    cancer = []
    for layer_name, X in layers.items():
        ct = _strip_prefix(layer_name.replace("_genes_pass1", ""))
        col = prop_clean.get(ct)
        mean_prop = float(prop_df[col].mean()) if col else 0.0
        if mean_prop < LEGACY_MIN_MEAN_PROP:
            continue
        if ct.lower().strip() not in LEGACY_CANCER_NAMES:
            continue
        progs = legacy_nmf_top50(X, gene_names)
        if progs is None:
            continue
        cancer.extend({"program_idx": k, "top_genes": g} for k, g in enumerate(progs))
    return cancer


def run_competitor(method: str, specimen: str, root: Path, legacy_root=None) -> dict:
    sdir = Path(root) / specimen
    layers_path = next((p for p in (sdir / f"{specimen}_gex_layers.h5ad", sdir / "gex_layers.h5ad",
                                    sdir / "c2l_gex_layers.h5ad") if p.exists()), None)
    prop_path = next((p for p in (sdir / f"{specimen}_proportions.csv", sdir / "proportions.csv") if p.exists()), None)
    if layers_path is None:
        return _typed_failure(specimen, method, "layers_missing", f"no GEX layers h5ad in {sdir}")
    if prop_path is None:
        return _typed_failure(specimen, method, "proportions_missing", f"no proportions.csv in {sdir}")
    import scanpy as sc
    t0 = time.time()
    adata = sc.read_h5ad(layers_path)
    adata.var_names_make_unique()
    prop_df = pd.read_csv(prop_path, index_col=0)
    common = adata.obs_names.intersection(prop_df.index)
    if len(common) == 0:
        return _typed_failure(specimen, method, "barcode_mismatch", "no shared barcodes between layers and proportions")
    if len(common) < adata.n_obs:
        adata = adata[common].copy()
        prop_df = prop_df.loc[common]
    layers = {name: adata.layers[name] for name in adata.layers}
    cancer = competitor_cancer_programs(layers, list(adata.var_names), prop_df)
    cls = classify_specimen(cancer)
    extra = {"seconds": round(time.time() - t0, 1)}
    if legacy_root is not None:
        lp = Path(legacy_root) / method / "programs" / f"{specimen}_programs.json"
        if lp.exists():
            legacy = json.loads(lp.read_text())
            lcancer = [p for ct, ps in legacy.get("programs", {}).items()
                       if ct.lower().strip() in LEGACY_CANCER_NAMES for p in ps]
            extra["matches_legacy_top50"] = bool(
                [p["top_50_genes"] for p in lcancer] == [p["top_genes"] for p in cancer])
        else:
            extra["matches_legacy_top50"] = None
    return _ok_record(specimen, method, cls, extra)


# --------------------------------------------------------------------------- CITEgeist side (frozen Module 4)
def run_citegeist(specimen: str, layers_dir: Path, module3_dir: Path, types=None) -> dict:
    """Frozen Module 4 on the A12.6 layers (same call, parameters and seeding as A12.11), keeping top-50 genes."""
    import nmf_module4 as N
    types = list(N.TYPES if types is None else types)
    t0 = time.time()
    method = "CITEgeist"
    qp_path = Path(module3_dir) / f"{specimen}_cell_prop_finetuned_results.csv"
    if not qp_path.exists():
        return _typed_failure(specimen, method, "proportions_missing", str(qp_path))
    manifest_path = Path(layers_dir) / "manifest.json"
    if not manifest_path.exists():
        return _typed_failure(specimen, method, "layers_missing", str(manifest_path))
    manifest = json.loads(manifest_path.read_text())
    if not manifest.get("conservation", {}).get("passed", False):
        return _typed_failure(specimen, method, "layers_missing", "conservation gate failed")
    if manifest.get("types") != types:
        return _typed_failure(specimen, method, "layer_type_mismatch", f"layer types {manifest.get('types')}")
    for aux in ("genes.txt", "barcodes.txt", "coords.csv"):
        if not (Path(layers_dir) / aux).is_file():
            return _typed_failure(specimen, method, "layers_missing", f"missing {aux}")
    import anndata as ad
    import scipy.sparse as sp
    layers_dir = Path(layers_dir)
    genes = [ln for ln in (layers_dir / "genes.txt").read_text().splitlines() if ln.strip()]
    barcodes = [ln for ln in (layers_dir / "barcodes.txt").read_text().splitlines() if ln.strip()]
    xy = pd.read_csv(layers_dir / "coords.csv").set_index("barcode").reindex(barcodes)[["x", "y"]].to_numpy(float)
    prop_df = pd.read_csv(qp_path, index_col=0)
    if list(prop_df.columns) != types:
        return _typed_failure(specimen, method, "layer_type_mismatch", f"proportion columns {list(prop_df.columns)}")
    prop_df = prop_df.loc[barcodes]
    adata = ad.AnnData(X=sp.csr_matrix((len(barcodes), len(genes)), dtype=np.float64),
                       obs=pd.DataFrame(index=pd.Index(barcodes, name="barcodes")),
                       var=pd.DataFrame(index=pd.Index(genes, name="genes")))
    adata.obsm["spatial"] = xy
    for t in types:
        try:
            layer = sp.load_npz(layers_dir / f"layer_{t}.npz").tocsr()
        except Exception as e:  # noqa: BLE001
            return _typed_failure(specimen, method, "layer_type_mismatch", f"layer_{t}.npz unreadable: {e}")
        if layer.shape != adata.shape:
            return _typed_failure(specimen, method, "layer_type_mismatch", f"{t} shape {layer.shape}")
        adata.layers[f"{t}{N.LAYER_PATTERN}"] = layer
    np.random.seed(N.RANDOM_STATE)
    result, fallback = N._nmf_fallback_during(lambda: N.discover_programs_from_layers(
        adata, cell_type_proportions=prop_df, profile_discovery_result=None, protein_adata=None,
        K_programs=N.K_PROGRAMS, lambda_spatial=N.LAMBDA_SPATIAL, lambda_sparsity=N.LAMBDA_SPARSITY,
        min_expression_threshold=N.MIN_EXPRESSION_THRESHOLD, validate_with_proteins=True,
        top_n_genes=N.TOP_N_GENES, detect_subpopulations=True, n_subpopulations=3, random_state=N.RANDOM_STATE))
    if fallback:
        return _typed_failure(specimen, method, "nmf_failed_random_fallback", "Module 4 substituted random matrices",
                              status="failed")
    entry = result.results_by_anchor.get(CITEGEIST_CANCER)
    if entry is None or entry.n_spots_used < 20 or len(entry.programs) == 0:
        return _typed_failure(specimen, method, "cancer_cells_absent",
                              f"{CITEGEIST_CANCER} missing/skipped (n_spots_used="
                              f"{0 if entry is None else int(entry.n_spots_used)})")
    active = N._active_mask(adata.layers[f"{CITEGEIST_CANCER}{N.LAYER_PATTERN}"])
    cancer = []
    for k, prog in enumerate(entry.programs):
        st = N._program_status(prog, entry, active, k)
        if st["status"] == "nonfinite_program":
            return _typed_failure(specimen, "CITEgeist", "nonfinite_program", f"program {k}", status="failed")
        cancer.append({"program_idx": int(prog.program_id), "top_genes": list(prog.top_genes),
                       "coherent": st["coherent"], "_top8": list(prog.top_genes[:8]),
                       "_moran_i": st["spatial_moran_i"], "_moran_p": st["spatial_moran_pvalue"]})
    cls = classify_specimen(cancer)
    for rp, p in zip(cls["programs"], cancer):
        rp["spatial_moran_i"], rp["spatial_moran_pvalue"] = p["_moran_i"], p["_moran_p"]
    extra = {"seconds": round(time.time() - t0, 1), "_top8": [p["_top8"] for p in cancer]}
    return _ok_record(specimen, method, cls, extra)


def cross_check_a1211(records, counts_path: Path) -> dict:
    """Compare the CITEgeist side with the A12.11 table (same specimens, programs, coherence flags, top-8 genes)."""
    ref = json.loads(Path(counts_path).read_text())
    by = {r["specimen"]: r for r in ref["specimens"]}
    mismatches = []
    for r in records:
        if r["status"] != "ok":
            continue
        a = by.get(r["specimen"])
        if a is None or a["status"] != "ok":
            mismatches.append({"specimen": r["specimen"], "why": "absent in A12.11"})
            continue
        progs = a["types"][CITEGEIST_CANCER]["programs"]
        same = (len(progs) == len(r["programs"]) and all(
            pa["coherent"] == pr["coherent"] and pa["top_genes"] == t8
            for pa, pr, t8 in zip(progs, r["programs"], r["_top8"])))
        if not same:
            mismatches.append({"specimen": r["specimen"], "why": "programs differ from A12.11"})
    return {"a12_11_counts": str(counts_path), "a12_11_run_id": ref.get("run_id"), "n_mismatch": len(mismatches),
            "mismatches": mismatches}


def legacy_summary(methods: dict) -> dict:
    """Per competitor: how many ok specimens' recomputed cancer top-50 lists equal the legacy program JSONs."""
    out = {}
    for m in COMPETITORS:
        recs = methods[m]
        out[m] = {"n_match": sum(1 for r in recs if r.get("matches_legacy_top50") is True),
                  "n_mismatch": sum(1 for r in recs if r.get("matches_legacy_top50") is False),
                  "n_unchecked": sum(1 for r in recs if r["status"] == "ok" and r.get("matches_legacy_top50") is None)}
    return out


def assemble(proto, run_id, args) -> dict:
    layers_base = Path(args.layers_root) / proto["_protocol_sha12"] / args.layers_run_id
    methods = {}
    for method, root in (("cell2location", args.c2l_root), ("tangram", args.tangram_root)):
        recs = []
        for i, spec in enumerate(SPECIMENS):
            print(f"[{method} {i + 1}/{len(SPECIMENS)}] {spec}", flush=True)
            recs.append(run_competitor(method, spec, root, args.legacy_programs_root))
        methods[method] = recs
    cg = []
    for i, spec in enumerate(SPECIMENS):
        print(f"[CITEgeist {i + 1}/{len(SPECIMENS)}] {spec}", flush=True)
        cg.append(run_citegeist(spec, layers_base / spec, args.module3_dir, types_for_protocol(proto)))
    methods["CITEgeist"] = cg
    for recs in methods.values():
        for r in recs:
            if r["status"] == "ok":
                PIO.validate_typed({"status": "ok", "value": r["value"]})
            elif not r.get("reason"):
                raise PIO.GateFailure(f"{r['specimen']}/{r['method']} non-ok without a reason")
    cross = cross_check_a1211(cg, args.a12_11_counts) if args.a12_11_counts else None
    for r in cg:
        r.pop("_top8", None)
    tot = {m: totals(recs) for m, recs in methods.items()}
    legacy_match = legacy_summary(methods)
    return {
        "run_id": run_id,
        "rule": {"mdk_gene": "MDK", "top_n": TOP_N, "secretory_genes": SECRETORY_GENES,
                 "secretory_overlap_threshold": SECRETORY_OVERLAP_THRESHOLD,
                 "competitor_cancer_names": sorted(LEGACY_CANCER_NAMES), "citegeist_cancer_name": CITEGEIST_CANCER,
                 "full": "MDK in top-50 and >=6/11 secretory genes in top-50 of one cancer program",
                 "partial": "MDK in top-50 of a cancer program, no program with >=6 secretory genes",
                 "absent": "MDK in no cancer program top-50 (or no cancer programs)",
                 "coherence_rule_citegeist_only": {"moran_I_min": MORAN_I_MIN, "p_max": MORAN_P_MAX}},
        "legacy_competitor_nmf": {"K": LEGACY_K, "alpha_W": LEGACY_ALPHA, "init": "nndsvda",
                                  "max_iter": LEGACY_MAX_ITER, "random_state": LEGACY_RANDOM_STATE,
                                  "min_mean_proportion": LEGACY_MIN_MEAN_PROP, "min_active_spots": LEGACY_MIN_SPOTS},
        "citegeist_module4": {"K_programs": 5, "lambda_spatial": 0.0, "lambda_sparsity": 0.01, "top_n_genes": 50,
                              "random_state": 42, "permutation_seed": "np.random.seed(42) before each specimen",
                              "layers_root": str(layers_base), "module3_dir": str(args.module3_dir)},
        "specimens": SPECIMENS,
        "per_specimen": methods,
        "totals": tot,
        "n_specimens": len(SPECIMENS),
        "competitor_legacy_cross_check": legacy_match,
        "citegeist_a12_11_cross_check": cross,
        "legacy_reference_totals": {
            "cell2location": {"full": 0, "partial": 3, "absent": 9},
            "tangram": {"full": 0, "partial": 8, "absent": 4},
            "source": "output/competitor_gex/comparison/failure_mode_summary.csv (2026-05-29; analysis=programs_mdk)"},
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layers-root", type=Path, default=REPO_ROOT / "output" / "sace_layers_module3_7type")
    ap.add_argument("--layers-run-id", default=None, help="the B4 id; REQUIRED for a v42 protocol, legacy default r1")
    ap.add_argument("--module3-dir", type=Path, default=None, help="the stamped B1 dir; REQUIRED for a v42 protocol, legacy default output/module3_7type otherwise")
    ap.add_argument("--c2l-root", type=Path, required=True)
    ap.add_argument("--tangram-root", type=Path, required=True)
    ap.add_argument("--legacy-programs-root", type=Path, default=None,
                    help="output/competitor_gex (holds <method>/programs/<specimen>_programs.json) for the cross-check")
    ap.add_argument("--a12-11-counts", type=Path, default=None, help="A12.11 nmf_counts.json for the CITEgeist cross-check")
    ap.add_argument("--out-root", type=Path, required=True, help="results/<sha12>/<run_id> (namespaced, empty)")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    a = ap.parse_args(argv)
    proto = PIO.load_protocol(a.protocol)
    PF.guard_entry("B11", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.module3_dir,
                   layers_run=(Path(a.layers_root) / proto["_protocol_sha12"] / a.layers_run_id) if a.layers_run_id else None,
                   nmf_run=a.a12_11_counts.parent if a.a12_11_counts else None)
    a.module3_dir = PF.legacy_or_required(a.module3_dir, a.protocol, "--module3-dir", REPO_ROOT / "output" / "module3_7type")
    a.layers_run_id = PF.legacy_or_required(a.layers_run_id, a.protocol, "--layers-run-id", "r1")
    PIO.check_output_target(a.out_root, proto["_protocol_sha12"], a.run_id)
    stamp = PF.stage_stamp_begin("B11", a.protocol, Path(a.out_root) / OUT_NAME, Path(__file__).resolve(), a.run_id)  # v42: job-start stamp gate, before any compute
    story = assemble(proto, a.run_id, a)
    import provenance
    stamped = provenance.stamp_payload(PIO.with_protocol_hash(story, proto), PRODUCED_BY)
    target = Path(a.out_root) / OUT_NAME
    PIO.atomic_write_json(target, stamped)
    if stamp:
        stamp.finish(module3_dir=a.module3_dir, layers_run=Path(a.layers_root) / proto["_protocol_sha12"] / a.layers_run_id, nmf_run=a.a12_11_counts.parent)
    print("wrote", target, flush=True)
    print(json.dumps(story["totals"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
