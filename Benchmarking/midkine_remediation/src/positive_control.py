"""Positive-control COMMOT rerun (protocol hurdle.validation_controls.real_data_post_freeze.positive_control_commot; spec 9.4.7).

Run with the COMMOT env python (commot 0.0.3 needs numpy < 2; CITEgeist_env cannot import commot), ONE process, ONE ``import commot``
(about 100-220 s on NFS): the baseline and the injected COMMOT runs share it.

Steps (``run_positive_control``):
  1. build the compartments from the SACE layers (baseline); eligible cancer compartments; the primary label (``select_label`` on k_share);
  2. the ORIGINALLY defined eligible keratin-high spot ids are fixed from these baseline labels;
  3. the LR list is frozen once, from the baseline N-P object (``filter_lr_database``, min_cell_pct 0.01), and reused for the injected run;
  4. baseline COMMOT (R00 settings: N-P input, CellChat human, heteromeric, pathway_sum, ``dis_thr``) -> MDK-NCL sender score;
  5. ``controls.inject_mdk_layers`` multiplies the LINEAR allocated MDK counts of the cancer type by ``fold`` (3.0) in those spots, the
     compartments are rebuilt, and the rows are asserted identical (labels, eligibility and covariates are NOT recomputed; k_share excludes MDK);
  6. injected COMMOT with the same frozen LR list and ``dis_thr``;
  7. the same estimator (``descriptive.fit_magnitude``: GMR among active senders) on the same eligible rows, labels, blocks and covariates;
  8. ``controls.injection_report``: baseline GMR, injected GMR, delta log(GMR). Typed not_estimable if either run is; no gate, no p-value.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Callable, Dict, Optional

import numpy as np
import pandas as pd

import commot_gates as G
import compartments as C
import controls as K
import descriptive as D
import hurdle as H

PAIR = "s-MDK-NCL"
MIN_CELL_PCT = 0.01
DIS_THR_PRIMARY = 500.0
CANCER_TYPE = "Cancer_Cells"
RNA_COVARIATE = "total"        # protocol deviation 2026-11-06 (PROTOCOL_DEVIATION_2026-11-06_rna_total.md)
SenderFn = Callable[..., dict]


# --------------------------------------------------------------------------- COMMOT wrappers (need the COMMOT env)
def _adata(X, spatial, genes, obs_names):
    import anndata as ad
    import scipy.sparse as sp

    a = ad.AnnData(X=sp.csr_matrix(np.asarray(X, dtype=np.float32)), obs=pd.DataFrame(index=list(obs_names)), var=pd.DataFrame(index=list(genes)))
    a.obsm["spatial"] = np.asarray(spatial, float)
    return a


def commot_frozen_lr(X, spatial, genes, obs_names, min_cell_pct: float = MIN_CELL_PCT) -> pd.DataFrame:
    """LR list frozen from the primary N-P object (CellChat, human, heteromeric)."""
    import commot as ct

    a = _adata(X, spatial, genes, obs_names)
    db = ct.pp.ligand_receptor_database(species="human", database="CellChat")
    return ct.pp.filter_lr_database(db, a, heteromeric=True, min_cell_pct=min_cell_pct).reset_index(drop=True)


def commot_sender(X, spatial, genes, obs_names, lr: pd.DataFrame, dis_thr: float) -> dict:
    """One COMMOT run -> typed result. ``ok`` carries the MDK-NCL native sender score (one value per row) after the G5 key-mapping check
    for that pair; a crash (e.g. zero-size reduction when no pair lies inside ``dis_thr``) is a typed ``failed`` result."""
    def go():
        import commot as ct

        a = _adata(X, spatial, genes, obs_names)
        ct.tl.spatial_communication(a, database_name="cellchat", df_ligrec=lr, dis_thr=dis_thr, heteromeric=True, pathway_sum=True)
        sender = a.obsm["commot-cellchat-sum-sender"]
        keys = [k for k in a.obsp.keys() if k.startswith("commot-cellchat-") and k.count("-") >= 3 and not k.endswith("total-total")]
        rs = {k: np.asarray(a.obsp[k].sum(axis=1)).ravel() for k in keys}
        gate = G.key_mapping_check(rs, {c: sender[c].to_numpy() for c in sender.columns}, [tuple(map(str, r[:2])) for r in lr.values])
        if not gate["pairs"][PAIR]["ok"]:                                  # only the pair this control reads has to pass (the other four are descriptive endpoints)
            raise RuntimeError(f"key-mapping gate failed for {PAIR}: {gate['pairs'][PAIR]}")
        return np.asarray(sender[PAIR].to_numpy(), float)

    r = G.run_commot_typed(go)
    return {"status": "ok", "scores": r["result"]} if r["status"] == "ok" else {"status": "failed", "reason": r["reason"]}


# --------------------------------------------------------------------------- positive control
def _mag(run: dict, el, group, cov, blocks) -> dict:
    if run["status"] != "ok":
        return D._nonest("commot_run_failed", f"COMMOT run failed: {run.get('reason', '')}")
    return D.fit_magnitude(np.asarray(run["scores"], float)[el], group, cov, blocks)


def run_positive_control(layers: Dict[str, np.ndarray], genes, qp: pd.DataFrame, nuclei: pd.Series, spot_umi: pd.Series, coords, M: float, *,
                         lr_fn: Callable[..., pd.DataFrame], sender_fn: SenderFn, dis_thr: float = DIS_THR_PRIMARY, fold: float = K.INJECT_FOLD,
                         cancer_type: str = CANCER_TYPE, p_construct: float = 0.05) -> dict:
    """See the module docstring. ``lr_fn(X, spatial, genes, obs_names) -> LR frame`` and ``sender_fn(X, spatial, genes, obs_names, lr, dis_thr) ->
    {'status': 'ok', 'scores': ndarray} | {'status': 'failed', 'reason'}`` are the COMMOT wrappers above (injectable so the logic is testable without commot)."""
    genes = list(genes)
    base = C.build_compartments(layers, genes, qp, nuclei, spot_umi, coords, M, p_construct=p_construct)
    el = C.eligible_cancer(base.obs, cancer_type).to_numpy()
    obs = base.obs[el]
    out = {"control": "positive_commot_mdk_injection", "fold": float(fold), "dis_thr": float(dis_thr), "pair": PAIR, "gates_branch": False,
           "rna_covariate": RNA_COVARIATE, "covariates": list(H.primary_covariates(RNA_COVARIATE)), "protocol_deviation": H.deviation_record(),
           "n_compartments": int(len(base.obs)), "n_eligible": int(el.sum()), "estimator": "descriptive.fit_magnitude (GMR among active senders), point estimate"}
    if el.sum() < D.MIN_ELIGIBLE:
        return {**out, **D._nonest("group_size", f"{int(el.sum())} eligible compartments < {D.MIN_ELIGIBLE}"), "status": "not_estimable"}
    lab = C.select_label(obs["k_share"].to_numpy())
    out["label"] = {k: v for k, v in lab.items() if k != "high"}
    if lab["status"] != "ok":
        return {**out, "status": "not_estimable", "reason_code": lab["reason_code"], "reason": lab["reason"]}
    group = lab["high"].astype(int)
    high_spots = obs["spot_id"].to_numpy()[lab["high"]]                    # fixed from the BASELINE labels
    out["n_injected_spots"] = int(len(high_spots))

    names = list(base.obs.index)
    lr = lr_fn(base.X["N_P"], base.spatial, genes, names)
    out["n_lr_pairs"] = int(len(lr))
    xy = base.spatial[el]
    blocks = H.make_blocks(xy, H.median_pitch(xy), H.BLOCK_SIDE_PITCHES, H.BLOCK_MIN_N)
    cov = H.covariate_matrix(obs, rna=RNA_COVARIATE)

    run0 = sender_fn(base.X["N_P"], base.spatial, genes, names, lr, dis_thr)
    inj_layers = K.inject_mdk_layers(layers, genes, cancer_type, high_spots, fold)
    inj = C.build_compartments(inj_layers, genes, qp, nuclei, spot_umi, coords, M, p_construct=p_construct)
    assert list(inj.obs.index) == names, "injection changed the compartment set"
    inj_ks = inj.obs["k_share"].to_numpy()
    base_ks = base.obs["k_share"].to_numpy()
    # k_share is mathematically invariant under the MDK injection (denom = tot - MDK = kcount + other), but computed
    # as E.sum(axis=1) - MDK it carries ~2 ULP of float64 round-off when MDK is large (real data: max rel diff ~2.6e-16).
    # The binding "labels must not move" invariant is membership (which label each eligible compartment gets), which the
    # second assert checks; the allclose bound just confirms k_share itself did not drift beyond float round-off.
    assert np.allclose(inj_ks, base_ks, rtol=1e-9, atol=0.0), "injection changed k_share beyond float round-off"
    inj_lab = C.select_label(inj.obs[el]["k_share"].to_numpy())
    assert inj_lab["status"] == "ok" and np.array_equal(inj_lab["high"], lab["high"]), "injection changed the k_share label (membership must not move)"
    assert np.allclose(H.covariate_matrix(inj.obs[el], rna=RNA_COVARIATE), cov), "injection changed the covariates"
    run1 = sender_fn(inj.X["N_P"], inj.spatial, genes, names, lr, dis_thr)

    base_mag, inj_mag = _mag(run0, el, group, cov, blocks), _mag(run1, el, group, cov, blocks)
    rep = K.injection_report(base_mag, inj_mag)
    out.update({"baseline": base_mag, "injected": inj_mag, "report": rep, "status": rep["status"], **({"value": rep["delta_log_gmr"]} if rep["status"] == "ok" else {}),
                **({"reason_code": rep["reason_code"], "reason": rep["reason"]} if rep["status"] != "ok" else {})})
    D.assert_no_nan_ok(out)
    return out


# --------------------------------------------------------------------------- toy input (also used by the unit tests)
TOY_GENES = list(C.KERATINS) + ["MDK", "NCL", "SDC2", "SDC4", "LRP1", "ITGA6", "ITGB1", "GA", "GB"]


def toy_inputs(nx: int = 36, ny: int = 30, seed: int = 0):
    """Synthetic specimen: spots on a grid, Cancer_Cells (0.4-0.9) + Other (the rest), keratin-high blobs, MDK/receptors expressed everywhere."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(nx), np.arange(ny))
    coords = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    S = len(coords)
    bc = [f"b{i}" for i in range(S)]
    f = np.sin(coords[:, 0] / 3 + 0.7) + np.sin(coords[:, 1] / 3 + 1.9)
    high = f > 0.4
    pc = rng.uniform(0.4, 0.9, S)                                          # varying qp_prop (a constant covariate would be collinear)
    qp = pd.DataFrame({CANCER_TYPE: pc, "Other": 1.0 - pc}, index=bc)
    g = {t: rng.poisson(40.0, (S, len(TOY_GENES))).astype(float) + 1.0 for t in qp.columns}
    ki = [TOY_GENES.index(k) for k in C.KERATINS]
    g[CANCER_TYPE][:, ki] = np.where(high[:, None], rng.poisson(8.0, (S, len(ki))), rng.poisson(0.05, (S, len(ki)))).astype(float)
    g[CANCER_TYPE][:, TOY_GENES.index("MDK")] = rng.poisson(15.0, S) + 2.0
    nuc = pd.Series(np.full(S, 6.0), index=bc)
    umi = pd.Series(rng.uniform(1500, 2500, S), index=bc)
    return g, TOY_GENES, qp, nuc, umi, coords, float(np.median(umi))


def real_main(a) -> int:  # pragma: no cover - real data, sbatch only (envs/COMMOT)
    """Real-data entry: P4-S2 (or any specimen) from the A12.6 layers and the A12.7 pooled M; COMMOT sees only the LR genes and only the five MDK
    pairs (per-pair independence proved by ``commot_runs --verify-subset``); the frozen LR list is the full filtered list from the baseline object."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import chain_guard as PF

    # v42 job-start checks first, before the heavy imports (commot, squidpy stack)
    PF.guard_entry("B7", a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.qp_dir, layers_run=a.layers_root, comp_run=a.compartments_run)
    import build_compartment_adata as B
    import commot as ct
    import commot_runs as R
    import protocol_io as PIO
    import rebuild_single_cell as RB

    proto = PIO.load_protocol(a.protocol)
    a.qp_dir = PF.legacy_or_required(a.qp_dir, a.protocol, "--qp-dir", Path(__file__).resolve().parents[3] / "output" / "module3_7type")
    if tuple(Path(a.out_root).parts[-2:]) != (proto["_protocol_sha12"], a.run_id):
        raise PIO.OutputTargetError(f"out-root must end in <{proto['_protocol_sha12']}>/<{a.run_id}>")
    target = PIO.check_output_target(Path(a.out_root) / a.specimen, proto["_protocol_sha12"], a.run_id)
    M = json.loads((Path(a.compartments_run) / "M.json").read_text())["M"]
    layers, genes, qp_full, qp, nuc, umi, coords = B.load_inputs(Path(a.layers_root) / a.specimen, Path(a.qp_dir) / RB.SPECIMENS[a.specimen]["qp_csv"])
    need = R.lr_gene_set(ct.pp.ligand_receptor_database(species="human", database="CellChat"))
    gi = [i for i, g in enumerate(genes) if g in need]
    lr_genes = [genes[i] for i in gi]

    def lr_fn(X, spatial, genes_, names):
        return commot_frozen_lr(X[:, gi], spatial, lr_genes, names)

    def sender_fn(X, spatial, genes_, names, lr, dis_thr):
        return commot_sender(X[:, gi], spatial, lr_genes, names, lr, dis_thr)

    res = run_positive_control(layers, genes, qp, nuc, umi, coords, M, lr_fn=lr_fn, sender_fn=sender_fn, dis_thr=a.dis_thr)
    res = {**res, "specimen": a.specimen}
    sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
    import provenance

    stamped = provenance.stamp_payload(PIO.with_protocol_hash(res, proto), "Benchmarking/midkine_remediation/src/positive_control.py")
    PIO.atomic_write_json(target / "positive_control.json", stamped)
    print("positive control", a.specimen, res["status"], res.get("value"), flush=True)
    return 0


def main(argv=None) -> int:
    import argparse

    ap = argparse.ArgumentParser(description="Positive-control COMMOT rerun (COMMOT env). --toy runs the built-in synthetic specimen.")
    ap.add_argument("--toy", action="store_true", help="synthetic specimen")
    ap.add_argument("--dis-thr", type=float, default=None, help="toy default 3 (grid pitch 1); real data default 500 (protocol)")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--specimen", help="real data: specimen id (needs the options below)")
    ap.add_argument("--layers-root", type=Path)
    ap.add_argument("--compartments-run", type=Path, help="A12.7 output dir (holds M.json)")
    ap.add_argument("--out-root", type=Path, help="output/positive_control_module3_7type/<sha12>/<run_id>")
    ap.add_argument("--run-id")
    ap.add_argument("--qp-dir", type=Path, default=None, help="the stamped B1 dir; REQUIRED for a v42 protocol, legacy default output/module3_7type otherwise")
    ap.add_argument("--protocol", type=Path, default=Path(__file__).resolve().parents[1] / "protocol.json")
    args = ap.parse_args(argv)
    if not args.toy:
        for req in ("specimen", "layers_root", "compartments_run", "out_root", "run_id"):
            if getattr(args, req) is None:
                raise SystemExit(f"real-data run needs --{req.replace('_', '-')} (or use --toy)")
        args.dis_thr = DIS_THR_PRIMARY if args.dis_thr is None else args.dis_thr
        return real_main(args)
    args.dis_thr = 3.0 if args.dis_thr is None else args.dis_thr
    layers, genes, qp, nuc, umi, coords, M = toy_inputs()
    res = run_positive_control(layers, genes, qp, nuc, umi, coords, M,
                               lr_fn=commot_frozen_lr, sender_fn=lambda X, sp, ge, nm, lr, dt: commot_sender(X, sp, ge, nm, lr, dt), dis_thr=args.dis_thr)
    txt = json.dumps(res, default=lambda o: o.item() if hasattr(o, "item") else list(o))
    if args.out:
        Path(args.out).write_text(txt)
    print("JSON_BEGIN" + txt + "JSON_END")
    return 0


if __name__ == "__main__":
    sys.exit(main())
