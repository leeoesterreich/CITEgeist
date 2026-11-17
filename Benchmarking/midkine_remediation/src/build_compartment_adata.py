"""A12.7: compartment AnnData for the 12 specimens from the A12.6 SACE layers (protocol ``compartments`` / ``normalisation``).

For each specimen: spot x type allocated counts -> ``compartments.build_compartments`` (p_construct 0.05, n_floor 0.5; N-P primary in X, N_A, N_S,
N_P_nf, E as sparse float32 layers; covariate/label columns in obs; coordinates in obsm['spatial']).
Gates and reports (labels only, no outcomes):
  * G3: rows per type == independent recount on the QP csv restricted to the barcodes of the raw Visium h5 (read here, not taken from the layers),
    after adding back spots dropped for zero depth; plus the sanity bracket 0.99*literal <= rows <= literal from the protocol;
  * eligible cancer compartments (qp_prop >= 0.10 and allocated total >= 200), nuclei-floor / n_floor binding fractions (``heavy_report``);
  * the keratin-high label rule from k_share counts only (``select_label``): rule, threshold, n_high / n_low, typed status.
The pooled depth reference M is the median raw spot UMI over all spots of all 12 specimens present in both the QP csv and the GEX object.
P4-S2 only: sensitivity builds p_construct 0.02/0.10/0.20 and n_floor 0.25/1.0 (separate h5ads, same code path).
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(HERE.parent))
import chain_guard as PF  # noqa: E402
import provenance  # noqa: E402
import compartments as C  # noqa: E402
import protocol_io as PIO  # noqa: E402

SPECIMENS_ORDER = None   # taken from the protocol
HEADLINE = "HCC22-088-P4-S2_1i_rep"
SENSITIVITY = [("p0.02", dict(p_construct=0.02)), ("p0.10", dict(p_construct=0.10)), ("p0.20", dict(p_construct=0.20)),
               ("floor0.25", dict(n_floor=0.25)), ("floor1.0", dict(n_floor=1.0))]


def load_inputs(layer_dir: Path, qp_path: Path):
    """Layers (dense float64 per type), genes, qp (restricted to the layer barcodes), nuclei, umi, coords."""
    genes = (layer_dir / "genes.txt").read_text().split("\n")[:-1]
    barcodes = (layer_dir / "barcodes.txt").read_text().split("\n")[:-1]
    qp_full = pd.read_csv(qp_path, index_col=0)
    qp = qp_full.loc[barcodes]
    layers = {t: sp.load_npz(layer_dir / f"layer_{t}.npz").toarray() for t in qp.columns}
    coords_df = pd.read_csv(layer_dir / "coords.csv", index_col=0).loc[barcodes]
    umi = pd.read_csv(layer_dir / "spot_umi.csv", index_col=0)["spot_umi"].loc[barcodes]
    nuc = pd.read_csv(layer_dir / "nuclei.csv", index_col=0)["nuclei_count_raw"].reindex(barcodes)
    return layers, genes, qp_full, qp, nuc, umi, coords_df[["x", "y"]].to_numpy(float)


def h5_barcodes(h5_path: Path):
    import h5py

    with h5py.File(h5_path, "r") as f:
        return [b.decode() for b in f["matrix"]["barcodes"][:]]


def to_sparse_f32(a: np.ndarray) -> sp.csr_matrix:
    return sp.csr_matrix(a.astype(np.float32))


def write_adata(comp: C.Compartments, path: Path, extra_uns: dict) -> None:
    import anndata as ad

    a = ad.AnnData(X=to_sparse_f32(comp.X["N_P"]), obs=comp.obs.copy())
    a.var_names = pd.Index(comp.genes)
    a.obsm["spatial"] = np.asarray(comp.spatial, float)
    for k in ("N_A", "N_S", "N_P_nf"):
        a.layers[k] = to_sparse_f32(comp.X[k])
    a.layers["E"] = to_sparse_f32(comp.E)
    a.uns["specimen_params"] = {k: float(v) for k, v in comp.specimen.items()}
    a.uns["provenance"] = json.dumps(extra_uns)
    assert "spatial_distance" not in a.obsp
    a.write_h5ad(str(path))


def g3(comp: C.Compartments, qp_full: pd.DataFrame, h5_bc, p_construct: float, literal: dict | None, types, v42: bool = False) -> dict:
    exp = C.recount_rows(qp_full, h5_bc, p_construct)                 # independent recount, barcodes from the raw h5
    got = comp.obs.groupby("cell_type").size().reindex(types).fillna(0).astype(int)
    dropped = comp.specimen["n_dropped_zero_depth"]
    # a zero-depth spot drops at most one row per type it qualifies for: allow got <= exp and exp - got <= dropped
    ok_exact = bool((got <= exp.reindex(types)).all() and ((exp.reindex(types) - got) <= dropped).all())
    rec = {"expected_recount": {t: int(exp[t]) for t in types}, "rows": {t: int(got[t]) for t in types},
           "n_dropped_zero_depth_spots": int(dropped), "passed": ok_exact}
    if dropped == 0:
        rec["passed"] = bool((got == exp.reindex(types)).all())
    if literal is not None:
        c = int(got[types[0]])
        rec["cancer_rows"] = c
        if v42:
            # the protocol's literal counts are the v41 marker-share row counts (protocol.json is frozen): comparing v2 rows to them is meaningless,
            # so the informational bracket is reported as not applicable rather than as a mismatch. The binding G3 gate is the independent recount above.
            rec["literal_bracket"] = None
            rec["bracket_ok"] = None
            rec["bracket_status"] = "not_applicable_v42: protocol literal counts are v41 marker-share row counts"
        else:
            rec["literal_bracket"] = [0.99 * literal, literal]
            rec["bracket_ok"] = bool(0.99 * literal <= c <= literal)
    return rec


def process(specimen, proto, layer_root: Path, qp_dir: Path, comp_root: Path, M: float, data_root: Path, write_sens: bool) -> dict:
    t0 = time.time()
    import rebuild_single_cell as RB

    layers, genes, qp_full, qp, nuc, umi, coords = load_inputs(layer_root / specimen, qp_dir / RB.SPECIMENS[specimen]["qp_csv"])
    types = list(qp.columns)
    bc_h5 = h5_barcodes(data_root / RB.SPECIMENS[specimen]["visium_dir"] / "filtered_feature_bc_matrix.h5")
    comp = C.build_compartments(layers, genes, qp, nuc, umi, coords, M)
    rep = {"specimen": specimen, "n_spots_layers": int(len(qp)), "n_genes": len(genes), "n_compartments": int(len(comp.obs)),
           "rows_per_type": comp.obs.groupby("cell_type").size().reindex(types).fillna(0).astype(int).to_dict(), "specimen_params": comp.specimen,
           "below_gate": specimen in proto["specimens"]["below_gate_specimens"]}
    lit = proto["specimens"]["expected_cancer_compartments_stored_prop"]["ge_0.05"].get(specimen)
    rep["G3"] = g3(comp, qp_full, bc_h5, 0.05, lit, types, v42=PF.is_v42(proto["_protocol_sha12"]))
    if not rep["G3"]["passed"]:
        raise RuntimeError(f"G3 failed for {specimen}: {rep['G3']}")
    el = C.eligible_cancer(comp.obs)
    rep["eligible_cancer"] = int(el.sum())
    rep["heavy"] = C.heavy_report(comp)
    rep["n_eligible_by_type_rule"] = {"qp_prop>=0.10_and_total>=200": int(el.sum())}
    k = comp.obs.loc[el, "k_share"].to_numpy()
    rep["k_share_quantiles"] = {str(q): float(np.quantile(k, q)) for q in (0.1, 0.25, 0.5, 0.7, 0.9, 0.99)} if len(k) else {}
    rep["n_high_at_k0_0.01"] = int((k >= C.K0_PRIMARY).sum()) if len(k) else 0
    if len(k) >= 2 * C.MIN_PER_GROUP:
        lab = C.select_label(k)
        rep["label"] = {kk: vv for kk, vv in lab.items() if kk != "high"}
    else:
        rep["label"] = {"status": "not_estimable", "reason_code": "group_size", "reason": f"{int(el.sum())} eligible compartments < {2 * C.MIN_PER_GROUP}"}
    # label columns in obs (labels only; recomputable from k_share)
    comp.obs["eligible_cancer"] = el.to_numpy()
    comp.obs["keratin_high_primary"] = False
    if rep["label"]["status"] == "ok":
        comp.obs.loc[comp.obs.index[el.to_numpy()], "keratin_high_primary"] = lab["high"]
    comp.obs["below_gate"] = rep["below_gate"]
    prov = {"specimen": specimen, "protocol_sha256": proto["_protocol_sha256"], "M": M}
    comp_root.mkdir(parents=True, exist_ok=True)
    write_adata(comp, comp_root / f"{specimen}.h5ad", prov)
    if write_sens and specimen == HEADLINE:
        rep["sensitivity_files"] = []
        for tag, kw in SENSITIVITY:
            c2 = C.build_compartments(layers, genes, qp, nuc, umi, coords, M, **kw)
            e2 = C.eligible_cancer(c2.obs)
            c2.obs["eligible_cancer"] = e2.to_numpy()
            write_adata(c2, comp_root / f"{specimen}__{tag}.h5ad", {**prov, "variant": tag, **{a: b for a, b in kw.items()}})
            rep["sensitivity_files"].append({"tag": tag, "n_compartments": int(len(c2.obs)), "eligible_cancer": int(e2.sum())})
            del c2
    rep["seconds"] = round(time.time() - t0, 1)
    (comp_root / f"{specimen}_report.json").write_text(json.dumps(rep, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    return rep


def pooled_M(layer_root: Path, qp_dir: Path, specimens) -> float:
    import rebuild_single_cell as RB

    umis, qps = {}, {}
    for s in specimens:
        umis[s] = pd.read_csv(layer_root / s / "spot_umi.csv", index_col=0)["spot_umi"]
        qps[s] = pd.read_csv(qp_dir / RB.SPECIMENS[s]["qp_csv"], index_col=0)
    return C.pooled_depth_reference(umis, qps)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--layers-root", type=Path, required=True, help="output/sace_layers_module3_7type/<sha12>/<layers_run_id>")
    ap.add_argument("--out-root", type=Path, required=True, help="output/compartments_module3_7type/<sha12>/<run_id>")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--qp-dir", type=Path, default=None, help="the stamped B1 dir; REQUIRED for a v42 protocol, legacy default output/module3_7type otherwise")
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "protocol.json")
    ap.add_argument("--specimen", default=None, help="one specimen (default all 12)")
    a = ap.parse_args(argv)
    proto = PIO.load_protocol(a.protocol)
    PF.guard_entry("B5", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.qp_dir, layers_run=a.layers_root)
    a.qp_dir = PF.legacy_or_required(a.qp_dir, a.protocol, "--qp-dir", REPO_ROOT / "output" / "module3_7type")
    if tuple(a.out_root.parts[-2:]) != (proto["_protocol_sha12"], a.run_id):
        raise PIO.OutputTargetError(f"out-root must end in <{proto['_protocol_sha12']}>/<{a.run_id}>")
    PIO.check_output_target(a.out_root, proto["_protocol_sha12"], a.run_id)
    tmp = a.out_root.parent / f".tmp_{a.run_id}_{os.getpid()}"
    tmp.mkdir(parents=True)
    specs = proto["specimens"]["ordered"]
    M = pooled_M(a.layers_root, a.qp_dir, specs)
    (tmp / "M.json").write_text(json.dumps({"M": M, "definition": "median raw spot UMI over all spots of all 12 specimens present in both the QP csv and the GEX object"}))
    reps = {}
    for s in ([a.specimen] if a.specimen else specs):
        reps[s] = process(s, proto, a.layers_root, a.qp_dir, tmp, M, a.data_root, write_sens=True)
        print(s, reps[s]["n_compartments"], "eligible", reps[s]["eligible_cancer"], reps[s]["label"].get("rule"), flush=True)
    summary = {"run_id": a.run_id, "protocol_sha256": proto["_protocol_sha256"], "M": M, "layers_root": str(a.layers_root),
               "git_head": provenance.code_commit(REPO_ROOT),
               "specimens": {s: {"eligible_cancer": r["eligible_cancer"], "label": r["label"], "G3": r["G3"]["passed"]} for s, r in reps.items()}}
    (tmp / "a12_7_summary.json").write_text(json.dumps(summary, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    if a.out_root.exists():
        a.out_root.rmdir()
    os.rename(tmp, a.out_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
