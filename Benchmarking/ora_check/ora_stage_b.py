"""B19: corrected ORA Test B on the rebuilt breast spots (stamped stage; post-hoc, descriptive). Promoted from ``run_ora_rebuilt.py``.

Rule file: ``ora_rebuilt_DECISION_RULE.md`` (unchanged). Question: with the rebuilt spots, do the specificity genes of a cell type, ranked by cross-fit partial
Spearman correlation with the type's share, hit that type's cell-type terms (CellMarker_2024 ORA)? Frozen mapping and three post-hoc corrected mappings, the
frozen "current" call and a "collapsed" call, a Xenium calibration, the null conditions. A verdict label printed by the readout is a VALUE, never a gate.

What is unchanged from ``run_ora_rebuilt.py`` (exploratory run ``results_rebuilt_20261103_1001``): the ORA library code (``ora_lib``, ``ora_sense_lib``,
``ora_rebuilt_lib``, ``ora_corrected_mapping``), the breast and Xenium refits and gene rankings, seeds (``default_rng(20261014 + i)``), top-100 genes, two gene
folds, the null and w25/w50/w75 conditions, the call variants, the tolerances (G2_TOL, G3), the Xenium reproduction against the stored Test B rows.
``breast``, ``xenium``, ``summary``, ``changed_flags`` and ``render_md`` are copies with ONLY the changes listed in ``tests/test_ora_stage_b.py::ALLOWED_DIFF``.

What changed (inputs and I/O only; no design change):
  * inputs come from the CLI (``--module3-dir`` = B1 dir with its ``p0/``, ``--nuclei-dir`` = B2 dir); nothing is read from ``results_*`` as an input;
  * the gene universe is the committed list ``ora_stage_gene_universe.json`` (the column names of the earlier SACE layers; carries no result), hash-pinned in
    ``ora_stage_pins.json``;
  * the legacy gate (task 17: the OLD inputs against the stored sense calls) is NOT part of the stage: its inputs include the v3 stand-in, which the chain
    refuses. The reproduction proof is the fidelity test (stamped inputs against ``results_rebuilt_20261103_1001``);
  * G4 still compares p0 with the plain marker-share table of the earlier arm, as a fail-closed precondition on 12 files pinned by sha256 in
    ``ora_stage_pins.json`` (the files are recorded in the artifact, never in the stamp inputs);
  * the comparison with the OLD frozen Test B (``old_frozen_B``, change flags A and C_rebuild_vs_old_signs) is dropped: it prints a v41-era result;
  * gates raise ``GateFailure`` (exit 2, ``gate_failure_<tag>.json``, no artifact, no stamp) instead of ``sys.exit(3)``;
  * the stamp is written by ``summary`` after every task of the array left its marker.

Usage (the wrapper ``.commands/submit-v2-ora.sh`` fills every path from ``chain_preflight.py``):
  ora_stage_b.py task N --run-id ID --out-root DIR --module3-dir DIR --nuclei-dir DIR      N = 0..11 breast, 12..16 Xenium regions 0..4
  ora_stage_b.py summary          (same arguments)   checks the 17 markers, writes the artifact ``ora_test_b_results.json`` and its stamp
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "midkine_remediation"))
import chain_guard as PF  # noqa: E402
from ora_check import ora_corrected_mapping as M  # noqa: E402
from ora_check import ora_frozen_manifest as F  # noqa: E402
from ora_check import ora_lib as L  # noqa: E402
from ora_check import ora_rebuilt_lib as RL  # noqa: E402
from ora_check import ora_sense_lib as S  # noqa: E402
from ora_check import (  # noqa: E402  (exploratory runner: imported for the unchanged helpers Rec, write_new, write_tables, read_tables, check_vs_stored)
    run_ora_rebuilt as RR,
)
from ora_check import run_ora_sense as RS  # noqa: E402

STAGE = "B19"
ARTIFACT = "ora_test_b_results.json"
PRODUCED_BY = "Benchmarking/ora_check/ora_stage_b.py"
N_BREAST, N_XEN = RR.N_BREAST, RR.N_XEN
N_TASKS = N_BREAST + N_XEN
G2_TOL = RR.G2_TOL
PRIMARY = RS.PRIMARY
MAPPINGS, VARIANTS = RR.MAPPINGS, RR.VARIANTS
PINS_FILE = HERE / "ora_stage_pins.json"
UNIVERSE_FILE = HERE / "ora_stage_gene_universe.json"
RULE_FILE = HERE / "ora_rebuilt_DECISION_RULE.md"
sha = RR.sha
write_new = RR.write_new
Rec = RR.Rec
write_tables = RR.write_tables
read_tables = RR.read_tables
check_vs_stored = RR.check_vs_stored
S_sign = RR.S_sign


class GateFailure(RuntimeError):
    """A fail-closed gate failed; ``tag`` names it, ``info`` is what it saw."""

    def __init__(self, tag, info):
        self.tag, self.info = tag, info
        super().__init__(f"GATE FAILED {tag}: {json.dumps(info, default=str)[:600]}")


def fail(run, tag, info):
    raise GateFailure(tag, info)


# ------------------------------------------------------------------ configuration (inputs from the CLI, not from module constants)
# repo-relative homes of untracked data (constants, so that configure() can be called twice and with any namespace root)
GT_REL = RS.GT_PROP.relative_to(RS.REPO)
H5AD_REL = Path("Benchmarking/xenium_pseudovisium/data_protein_gt/h5ad_objects")


class Cfg:
    m3: Path = None  # B1 dir (v2b1), holds <s>_cell_prop_finetuned_results.csv and p0/
    nuc: Path = None  # B2 dir (v2b2)
    ns: Path = None  # namespace root: where the (untracked) frozen result files and the earlier marker-share table live
    frozen: Path = None  # <ns>/Benchmarking/ora_check


def configure(module3_dir, nuclei_dir, namespace_root=None):
    ns = Path(namespace_root or os.environ.get("V2_NAMESPACE_ROOT") or REPO)
    Cfg.m3, Cfg.nuc, Cfg.ns = Path(module3_dir), Path(nuclei_dir), ns
    Cfg.frozen = ns / "Benchmarking" / "ora_check"
    # UNTRACKED data read through modules that locate it from their own file (RS.GT_PROP, data_io.H5AD_DIR, sace_props_ablation.DATA) must come from the
    # namespace root, not from the code worktree (which holds committed files only). Committed code and libraries stay in the worktree.
    RS.GT_PROP = ns / GT_REL
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io

    data_io.H5AD_DIR = str(ns / H5AD_REL)


def data_paths() -> dict:
    """Every untracked-data path the task and summary steps resolve after ``configure``; the test requires all of them under the namespace root."""
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io

    return {
        "frozen_results_root": Cfg.frozen,
        "marker_share_dir": Cfg.ns / pins()["g4_reference"]["dir"],
        "xenium_gt": RS.GT_PROP,
        "xenium_h5ad_dir": Path(data_io.H5AD_DIR),
        "xenium_h5ad_cite_0": Path(data_io.h5ad_path(0, "CITE")),
    }


def pins() -> dict:
    return json.loads(PINS_FILE.read_text())


def gene_universe(sample: str) -> list:
    """Committed gene list of the SACE model of ``sample`` (hash-pinned in ora_stage_pins.json)."""
    want = pins()["gene_universe"]["sha256"]
    got = sha(UNIVERSE_FILE)
    if got != want:
        raise GateFailure("universe_pin", {"file": str(UNIVERSE_FILE), "sha256": got, "pinned": want})
    return json.loads(UNIVERSE_FILE.read_text())["genes"][sample]


def g4_reference(sample: str) -> Path:
    """Plain marker-share table of the earlier arm for ``sample``; refused unless it hashes to the pin (precondition only)."""
    p = pins()["g4_reference"]
    f = Cfg.ns / p["dir"] / p["file_pattern"].format(sample=sample)
    if not f.is_file():
        raise GateFailure("g4_reference_missing", {"file": str(f)})
    if sha(f) != p["sha256"][sample]:
        raise GateFailure("g4_reference_pin", {"file": str(f), "sha256": sha(f), "pinned": p["sha256"][sample]})
    return f


def frozen_check():
    problems = F.check(Cfg.frozen)
    if problems:
        raise GateFailure("frozen_files", {"problems": problems[:20]})


# ------------------------------------------------------------------ breast (copy of run_ora_rebuilt.breast, mode "rebuilt" only)
def breast_inputs(sample):
    return dict(
        p0=Cfg.m3 / "p0" / sample / f"{sample}_cell_prop_global_results.csv",
        v2=Cfg.m3 / f"{sample}_cell_prop_finetuned_results.csv",
        nuc=Cfg.nuc / sample / f"{sample}_nuclei_counts_fullres_patchwise.csv",
    )


def breast(idx, run):
    sys.path.insert(0, str(REPO / "Benchmarking" / "brightness_estimator"))
    import breast_predict as bp
    from ora_check.run_sace_layers import load_v41_profile

    from CITEgeist.model.deconvolution.rna_refinement_adata import run_rna_refinement

    run = Path(run)
    sample = RS.SPECIMENS[idx]
    out = run
    seed = RS.BASE_SEED + idx
    mapping = RS.load_mapping("term_mapping_7type.yaml")
    markers = L.marker_gene_set(load_v41_profile())
    inp = breast_inputs(sample)
    for k, p in inp.items():
        assert p.is_file(), f"missing input {k}: {p}"
    rec = Rec()
    p0_full = pd.read_csv(inp["p0"], index_col=0)
    nuc = pd.read_csv(inp["nuc"], index_col=0).iloc[:, 0]
    stored_v2 = pd.read_csv(inp["v2"], index_col=0)
    keep = p0_full.index[p0_full.index.isin(nuc.index)]
    ref = pd.read_csv(g4_reference(sample), index_col=0)  # G4
    common = p0_full.index.intersection(ref.index)
    g4 = {
        "same_spots": bool(set(stored_v2.index) == set(p0_full.index) == set(nuc.index) == set(ref.index)),
        "p0_vs_marker_share_max_abs": float((p0_full.loc[common] - ref.loc[common]).abs().max().max()),
        "n_common_with_marker_share": int(len(common)),
        "rows_sum_to_1_v2": bool(np.allclose(stored_v2.sum(axis=1), 1.0, atol=1e-6)),
        "rows_sum_to_1_p0": bool(np.allclose(p0_full.sum(axis=1), 1.0, atol=1e-6)),
        "columns_equal": bool(list(stored_v2.columns) == list(p0_full.columns)),
    }
    if not (
        g4["same_spots"]
        and g4["p0_vs_marker_share_max_abs"] == 0.0
        and g4["rows_sum_to_1_v2"]
        and g4["rows_sum_to_1_p0"]
        and g4["columns_equal"]
    ):
        fail(run, f"breast_{idx:02d}_G4", g4)
    spots = list(stored_v2.index)
    layer_genes = gene_universe(sample)
    gex_full = bp.load_raw_gex(sample, list(p0_full.index))
    gex_full.layers["raw_counts"] = gex_full.X.copy()
    gex = gex_full[list(keep)].copy()
    p0 = p0_full.loc[list(keep)]
    types = list(p0.columns)
    gn = list(gex.var_names)
    raw = gex.layers["raw_counts"]
    Y = np.asarray(raw.toarray() if hasattr(raw, "toarray") else raw, float)
    pos = {g: i for i, g in enumerate(gn)}
    obs_pos = {b: i for i, b in enumerate(gex.obs_names)}
    Yc = Y[np.array([obs_pos[b] for b in spots])]
    lib = Yc.sum(axis=1)
    ok = lib > 0
    uni_genes = [g for g in layer_genes if g in pos]
    mean = Yc[ok].mean(axis=0)
    expressed = [g for g in uni_genes if mean[pos[g]] >= S.MIN_MEAN]
    cand = [g for g in expressed if g not in markers]
    bg = [g.upper() for g in expressed]
    gi = np.array([pos[g] for g in uni_genes])
    Xl = np.log1p(1e4 * Yc[ok][:, gi] / lib[ok][:, None])
    ex, rl = S.prep_gene_ranks(Xl, lib[ok])
    folds = S.make_folds(len(gn), seed)
    gene_pos_in_uni = np.full(len(gn), -1)
    for i, g in enumerate(uni_genes):
        gene_pos_in_uni[pos[g]] = i
    meta = {
        "sample": sample,
        "mode": "rebuilt",
        "n_spots": len(spots),
        "n_spots_p0": int(len(p0_full)),
        "n_genes_universe": len(uni_genes),
        "n_expressed": len(expressed),
        "n_candidates": len(cand),
        "n_zero_nuclei_spots": int((nuc.loc[list(keep)] == 0).sum()),
        "G4": g4,
    }

    def refit(keepmask):
        g2 = gex[:, np.flatnonzero(keepmask)].copy()
        p, _ = run_rna_refinement(g2, p0, nuc.loc[list(keep)], n_hvg=2000, c=0.1, n_passes=3)
        return p.loc[spots].to_numpy(float)

    full = refit(np.ones(len(gn), dtype=bool))
    meta["repro_max_abs_v2"] = float(np.abs(full - stored_v2.loc[spots][types].to_numpy(float)).max())
    if meta["repro_max_abs_v2"] > G2_TOL:
        fail(run, f"breast_{idx:02d}_G2", meta)
    shares = {"v41": [p0.loc[spots].to_numpy(float)] * len(folds), "v2": []}
    for f in folds:
        km = np.ones(len(gn), dtype=bool)
        km[f] = False
        shares["v2"].append(refit(km))
    print("refits done", meta, flush=True)

    def rho_for(arm, transform):
        o = np.full((len(uni_genes), len(types)), np.nan)
        for k, f in enumerate(folds):
            sh = transform(shares[arm][k])[ok]
            fu = gene_pos_in_uni[f]
            fu = fu[fu >= 0]
            o[fu] = S.partial_rho(ex[:, fu], rl, sh)
        return pd.DataFrame(o, index=uni_genes, columns=types)

    rows, tables = [], {}
    for arm in ("v41", "v2"):
        evaluable = {t: all(S.evaluable_share(sh[:, j]) for sh in shares[arm]) for j, t in enumerate(types)}
        conds = {"real": lambda x: x, "null": lambda x: S.permute_fraction(x, 100, seed, null_seed=seed)}
        for w in (25, 50, 75):
            conds[f"w{w}"] = (lambda ww: (lambda x: S.permute_fraction(x, ww, seed + 1000 * ww)))(w)
        for cond, tf in conds.items():
            rho = rho_for(arm, tf)
            genes = RS.gene_lists_from_scores(rho, cand)
            for t in types:
                row, tb = rec.row(genes[t], bg, t, mapping)
                tk = f"breast|{arm}|{cond}|{sample}|{t}"
                tables[tk] = tb
                rows.append(
                    dict(
                        test="B",
                        dataset="breast",
                        arm=arm,
                        condition=cond,
                        unit=sample,
                        type=t,
                        evaluable=evaluable[t],
                        table_key=tk,
                        **row,
                    )
                )
                if cond == "null":
                    rows.append({**rows[-1], "condition": "w100"})  # w = 100 is the full permutation (rule); same table
        print("B done", sample, arm, flush=True)
    df = pd.DataFrame(rows)
    write_new(out / "calls" / f"B_breast_{idx:02d}.csv", df.to_csv(index=False))
    write_new(out / "calls" / f"B_breast_{idx:02d}_meta.json", json.dumps(meta, indent=2))
    write_tables(out / "tables" / f"B_breast_{idx:02d}.parquet", tables)
    inputs = {k: {"path": str(v), "sha256": sha(v)} for k, v in inp.items()}
    inputs["gene_universe"] = {"path": str(UNIVERSE_FILE), "sha256": sha(UNIVERSE_FILE)}
    inputs["g4_reference"] = {"path": str(g4_reference(sample)), "sha256": sha(g4_reference(sample))}
    raw_h5 = Path(bp.DATA_ROOT) / sample / "outs" / "filtered_feature_bc_matrix.h5"
    inputs["raw_h5_unstamped"] = {"path": str(raw_h5), "sha256": sha(raw_h5)}
    write_new(out / "inputs" / f"breast_{idx:02d}.json", json.dumps(inputs, indent=2))


# ------------------------------------------------------------------ Xenium (copy of run_ora_rebuilt.xenium)
def xenium(r, run):
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    run = Path(run)
    rec = Rec()
    mapping = RS.load_mapping("term_mapping_xenium7.yaml")
    markers = RS._xen_markers()
    seed = RS.XEN_SEED + r
    _, data = data_io.load_xenium(7)
    d = data[r]
    gt = RS._gt(r)
    obs = list(d["obs"])
    genes = [g.upper() for g in d["genes"]]
    p0 = np.asarray(rs.ms_share(d["X"], d["names"], rs.PROF[7]))
    Y = d["Y"]
    lib_all = Y.sum(axis=1)
    sel = np.array([i for i, b in enumerate(obs) if b in gt.index and lib_all[i] > 0])
    mean = Y[sel].mean(axis=0)
    expressed = [g for g, m in zip(genes, mean) if m >= S.MIN_MEAN]
    cand = [g for g in expressed if g not in markers]
    tb_frozen = RS._testable(genes, mapping)
    ex, rl = S.prep_gene_ranks(np.log1p(1e4 * Y[sel] / lib_all[sel][:, None]), lib_all[sel])
    folds = S.make_folds(len(genes), seed)
    shares = {"v41": [p0] * 2, "oracle": [gt.reindex(obs).to_numpy(float)] * 2, "v2": []}
    for f in folds:
        mask = np.ones(len(genes), dtype=bool)
        mask[f] = False
        p, _ = refine_proportions(Y, d["nc"], p0, c=0.1, n_passes=3, model="multinomial", gene_mask=mask)
        shares["v2"].append(p)

    def rho_for(arm, tf):
        o = np.full((len(genes), len(RS.XTYPES)), np.nan)
        for k, f in enumerate(folds):
            o[f] = S.partial_rho(ex[:, f], rl, tf(shares[arm][k][sel]))
        return pd.DataFrame(o, index=genes, columns=RS.XTYPES)

    rows, tables = [], {}
    for arm in ("oracle", "v2", "v41"):
        for cond, tf in {
            "real": lambda x: x,
            "null": lambda x: S.permute_fraction(x, 100, seed, null_seed=seed),
        }.items():
            gl = RS.gene_lists_from_scores(rho_for(arm, tf), cand)
            for j, t in enumerate(RS.XTYPES):
                row, tb = rec.row(
                    gl[t], expressed, t, mapping
                )  # ORA for every type; testability is applied per mapping later
                tk = f"xenium|{arm}|{cond}|region_{r}|{t}"
                tables[tk] = tb
                ev = bool(all(S.evaluable_share(sh[sel][:, j]) for sh in shares[arm]))
                rows.append(
                    dict(
                        test="B",
                        dataset="xenium",
                        arm=arm,
                        condition=cond,
                        unit=f"region_{r}",
                        type=t,
                        evaluable=ev,
                        testable=bool(tb_frozen[t] >= S.MIN_PANEL_GENES),
                        n_panel_genes=tb_frozen[t],
                        table_key=tk,
                        **row,
                    )
                )
    df = pd.DataFrame(rows)
    # reproduction against the stored Test B rows (frozen mapping; every type is testable there, so the stored rows have ORA for all)
    stored_csv = Cfg.frozen / "results_sense_20261101_1550" / "calls" / f"B_xenium_{r}.csv"
    stored = pd.read_csv(stored_csv, keep_default_na=False, na_values=[""])
    chk = check_vs_stored(df, stored_csv, f"xenium region {r} vs stored Test B")
    chk["stored_all_testable"] = bool(stored["testable"].astype(bool).all())
    chk["stored_file_sha256"] = sha(stored_csv)
    write_new(run / "calls" / f"B_xenium_{r}.csv", df.to_csv(index=False))
    write_new(
        run / "calls" / f"B_xenium_{r}_meta.json",
        json.dumps(
            {
                "region": r,
                "n_spots_selected": int(len(sel)),
                "n_expressed": len(expressed),
                "n_candidates": len(cand),
                "panel_genes": genes,
                "repro": chk,
            },
            indent=2,
        ),
    )
    write_tables(run / "tables" / f"B_xenium_{r}.parquet", tables)
    if chk["n_mismatch"] or not chk["compared"] or not chk["stored_all_testable"]:
        fail(run, f"xenium_{r}_repro", chk)


# ------------------------------------------------------------------ summary (copy of run_ora_rebuilt.summary / changed_flags / render_md)
def summary(run):
    run = Path(run)
    ov = M.load_overrides()
    mp_b, mp_x = RS.load_mapping("term_mapping_7type.yaml"), RS.load_mapping("term_mapping_xenium7.yaml")
    bcalls = pd.concat(
        [
            pd.read_csv(f, keep_default_na=False, na_values=[""])
            for f in sorted((run / "calls").glob("B_breast_[0-9][0-9].csv"))
        ],
        ignore_index=True,
    )
    xcalls = pd.concat(
        [
            pd.read_csv(f, keep_default_na=False, na_values=[""])
            for f in sorted((run / "calls").glob("B_xenium_[0-9].csv"))
        ],
        ignore_index=True,
    )
    for df in (bcalls, xcalls):
        df["condition"] = df["condition"].astype(str)
        for c in ("evaluable", "passed"):
            df[c] = df[c].astype(bool)
    if bcalls["unit"].nunique() != S.N_SPECIMENS or xcalls["unit"].nunique() != S.N_REGIONS:
        fail(
            run,
            "summary_units",
            {"breast_units": int(bcalls["unit"].nunique()), "xenium_units": int(xcalls["unit"].nunique())},
        )
    tables = {}
    for f in sorted((run / "tables").glob("B_*.parquet")):
        tables.update(read_tables(f))
    panel = {
        f"region_{r}": json.loads((run / "calls" / f"B_xenium_{r}_meta.json").read_text())["panel_genes"]
        for r in range(S.N_REGIONS)
    }
    gmt = RS.libs()[PRIMARY]
    repro = {"G2_max_abs_v2": {}, "G3": {}, "G4": {}}
    for f in sorted((run / "calls").glob("B_breast_[0-9][0-9]_meta.json")):
        m = json.loads(f.read_text())
        repro["G2_max_abs_v2"][m["sample"]] = m["repro_max_abs_v2"]
        repro["G4"][m["sample"]] = m["G4"]
    repro["xenium_repro"] = {
        r: json.loads((run / "calls" / f"B_xenium_{r}_meta.json").read_text())["repro"] for r in range(S.N_REGIONS)
    }

    res, allv = {}, []
    for mname, lv in MAPPINGS.items():
        fn = M.corrected_lineage_fn(ov, lv)
        lb, lx = RL.LineageCache(fn, mp_b), RL.LineageCache(fn, mp_x)
        tcount = {u: RL.testable_counts([g.upper() for g in panel[u]], gmt, mp_x, lx) for u in panel}
        for var in VARIANTS:
            b = RL.variant_table(bcalls, tables, mp_b, lb, var)
            x = RL.variant_table(xcalls, tables, mp_x, lx, var, testable=tcount)
            if mname == "frozen" and var == "current":  # G3: the table path reproduces the in-job calls exactly
                dp = int(
                    (b["passed"].to_numpy() != bcalls["passed"].to_numpy()).sum()
                    + (x["passed"].to_numpy() != xcalls["passed"].to_numpy()).sum()
                )
                dn = float(
                    max(
                        np.abs(b["nlp"].to_numpy(float) - bcalls["nlp"].to_numpy(float)).max(),
                        np.abs(x["nlp"].to_numpy(float) - xcalls["nlp"].to_numpy(float)).max(),
                    )
                )
                repro["G3"] = {"n_pass_differences": dp, "max_abs_nlp_difference": dn}
                if dp or dn > 1e-9:
                    fail(run, "G3", repro["G3"])
            r = S.readout(b, x)
            res[f"{mname}|{var}"] = {
                "readout": r,
                "xenium_testable_counts": tcount,
                "xenium_pass_totals": x[x["condition"] == "real"].groupby("arm")["passed"].sum().to_dict(),
                "xenium_null_oracle": int(x[(x["arm"] == "oracle") & (x["condition"] == "null")]["passed"].sum()),
                "breast_null_totals": {
                    a: int(b[(b["arm"] == a) & (b["condition"] == "null")]["passed"].sum()) for a in ("v41", "v2")
                },
                "breast_pass_by_type": {
                    a: b[(b["arm"] == a) & (b["condition"] == "real")].groupby("type")["passed"].sum().to_dict()
                    for a in ("v41", "v2")
                },
            }
            keep = ["dataset", "arm", "condition", "unit", "type", "evaluable", "passed", "best_rank", "nlp"]
            allv.append(pd.concat([b[keep], x[keep]]).assign(mapping=mname, variant=var))
    variants = pd.concat(allv, ignore_index=True)
    variants.to_csv(run / "variants.csv", index=False)
    flags = changed_flags(res)
    comp_libs = {
        a: {
            c: int(bcalls[(bcalls["arm"] == a) & (bcalls["condition"] == "real")][c].sum())
            for c in bcalls.columns
            if c.startswith("pass_")
        }
        for a in ("v41", "v2")
    }
    out = {"results": res, "changed_flags": flags, "comparison_library_passes_frozen_mapping_real": comp_libs}
    (run / "readout_rebuilt.json").write_text(
        json.dumps(out, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    )
    (run / "repro.json").write_text(json.dumps(repro, indent=2, default=str))
    (run / "readout_rebuilt.md").write_text(render_md(out, repro))
    print(json.dumps({k: v["readout"]["verdict"] for k, v in res.items()}), flush=True)
    print(json.dumps(flags, indent=1), flush=True)
    return out, repro, variants, bcalls, xcalls


def changed_flags(res):
    """Rule section 'what counts as a changed conclusion': (B) mapping, (C) per-type sign, (D) calibration. Flags A and C_rebuild_vs_old_signs compared with the OLD frozen Test B and are dropped."""
    f = {}
    base = res["frozen|current"]["readout"]
    types = sorted(base["per_type"]["v2"])

    def signs(r):
        return {t: S_sign(r["per_type"]["v2"][t] - r["per_type"]["v41"][t]) for t in types}

    for var in VARIANTS:
        ref = res[f"frozen|{var}"]["readout"]
        for m in ("L1", "L2", "L3"):
            r = res[f"{m}|{var}"]["readout"]
            v3r, v3c = ref["V3"], r["V3"]
            f[f"{m}|{var}"] = {
                "B_verdict": {
                    "frozen": ref["verdict"],
                    "corrected": r["verdict"],
                    "flag": "SAME" if ref["verdict"] == r["verdict"] else "CHANGED",
                },
                "C_signs": {
                    "frozen": signs(ref),
                    "corrected": signs(r),
                    "flag": "SAME" if signs(ref) == signs(r) else "CHANGED",
                },
                "D_calibration": {
                    k: {"frozen": v3r.get(k), "corrected": v3c.get(k)}
                    for k in ("status", "oracle_ge_methods", "gap_ok")
                },
            }
            f[f"{m}|{var}"]["D_calibration"]["flag"] = (
                "SAME"
                if all(
                    v["frozen"] == v["corrected"]
                    for v in f[f"{m}|{var}"]["D_calibration"].values()
                    if isinstance(v, dict)
                )
                else "CHANGED"
            )
    return f


def render_md(out, repro):
    res = out["results"]
    L_ = [
        "# Corrected ORA Test B on rebuilt breast spots (primary library CellMarker_2024) -- post-hoc mapping sensitivity",
        "",
        "Rule: `ora_rebuilt_DECISION_RULE.md`. Plausibility check, not an accuracy measure. The author decides; nothing here chooses a mapping.",
        "",
        "## Totals of 84, side by side",
        "",
        "| mapping | call | v41 | v2 | delta | verdict | V1 gap v41/v2 | V2 ok v41/v2 | V3 | null v41/v2 | Xenium oracle/v2/v41 (of 35) oracle null |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for k, v in res.items():
        r = v["readout"]
        m, c = k.split("|")
        xp = v["xenium_pass_totals"]
        L_.append(
            f"| {m} | {c} | {r['total']['v41']} | {r['total']['v2']} | {r['delta']} | {r['verdict']} | {r['V1']['v41']['gap']}/{r['V1']['v2']['gap']} | "
            f"{r['V2']['v41']['ok']}/{r['V2']['v2']['ok']} | {r['V3'].get('status')} | {v['breast_null_totals']['v41']}/{v['breast_null_totals']['v2']} | "
            f"{xp.get('oracle')}/{xp.get('v2')}/{xp.get('v41')} null {v['xenium_null_oracle']} |"
        )
    L_ += ["", "## Passes of 12 by type (v41 / v2), real shares", ""]
    types = sorted(res["frozen|current"]["readout"]["per_type"]["v2"])
    L_ += ["| type | " + " | ".join(res.keys()) + " |", "|---|" + "---|" * len(res)]
    for t in types:
        L_.append(
            f"| {t} | "
            + " | ".join(
                f"{v['readout']['per_type']['v41'][t]} / {v['readout']['per_type']['v2'][t]}" for v in res.values()
            )
            + " |"
        )
    L_ += [
        "",
        "## Changed-conclusion flags (rule; the comparison with the old frozen Test B is not part of the stage)",
        "",
        "```",
        json.dumps(out["changed_flags"], indent=1, default=str),
        "```",
        "",
        "## Reproduction (G2-G4, Xenium)",
        "",
        "```",
        json.dumps(repro, indent=1, default=str)[:6000],
        "```",
        "",
    ]
    return "\n".join(L_) + "\n"


# ------------------------------------------------------------------ stage driver: markers, artifact, stamp
def _git_head() -> str:
    return subprocess.run(
        ["git", "--no-optional-locks", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True
    ).stdout.strip()


def _clean(o):
    """JSON-safe copy: numpy scalars to Python, NaN/inf to None, tuples to lists, keys to str."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.generic):
        return _clean(o.item())
    if isinstance(o, float):
        return o if math.isfinite(o) else None
    return o


def records(df: pd.DataFrame) -> list:
    return _clean(json.loads(df.to_json(orient="records")))


def run_task(n, run, module3_dir, nuclei_dir):
    """One array task: ``n`` = 0..11 breast, 12..16 Xenium. Writes its outputs and then its marker ``tasks/task_NN.json`` (never overwritten)."""
    run = Path(run)
    configure(module3_dir, nuclei_dir)
    run.mkdir(parents=True, exist_ok=True)
    frozen_check()
    if not 0 <= n < N_TASKS:
        raise SystemExit(f"task id {n} out of range 0..{N_TASKS - 1}")
    t0 = time.time()
    breast(n, run) if n < N_BREAST else xenium(n - N_BREAST, run)
    write_new(
        run / "tasks" / f"task_{n:02d}.json",
        json.dumps(
            {
                "task": n,
                "slurm_job_id": os.environ.get("SLURM_JOB_ID") or "no-slurm",
                "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
                "code_commit": _git_head(),
                "wall_s": time.time() - t0,
            },
            indent=1,
        ),
    )


def config_block() -> dict:
    return {
        "descriptive_only": True,
        "post_hoc": True,
        "rule_file": RULE_FILE.name,
        "rule_file_sha256": sha(RULE_FILE),
        "seeds": {
            "breast_base": RS.BASE_SEED,
            "breast_rule": "default_rng-based folds and permutations with seed = base + specimen index",
            "xenium_base": RS.XEN_SEED,
        },
        "top_n_genes": S.TOP_N,
        "n_gene_folds": 2,
        "min_mean": S.MIN_MEAN,
        "adj_p": S.ADJ_P,
        "top_k_terms": S.TOP_K,
        "n_hvg": 2000,
        "c": 0.1,
        "n_passes": 3,
        "w_grid_partial_permutation": [25, 50, 75],
        "g2_tol": G2_TOL,
        "g3_tol_nlp": 1e-9,
        "min_panel_genes": S.MIN_PANEL_GENES,
        "primary_library": PRIMARY,
        "mappings": list(MAPPINGS),
        "calls": list(VARIANTS),
        "labels_printed_as_values": "readout verdict (V2_WORSE, MIXED, ...) and V1/V2/V3 statuses are values of the readout, not gates",
        "specimens": list(RS.SPECIMENS),
        "xenium_types": list(RS.XTYPES),
        "removed_from_exploratory": [
            "legacy gate (task 17)",
            "old_frozen_B and change flags A, C_rebuild_vs_old_signs (print a v41-era result)",
        ],
        "gene_universe_pin": pins()["gene_universe"],
        "libraries": {
            k: v["sha256"]
            for k, v in RS.R.deep_merge(RS.R.PROTO, json.loads((HERE / "ora_protocol_v41_7type.json").read_text()))[
                "ora"
            ]["libraries"].items()
        },
    }


def build_artifact(run, run_id, out, repro, variants, bcalls, xcalls, markers) -> dict:
    sib = {
        name: sha(run / name) for name in ("variants.csv", "readout_rebuilt.json", "readout_rebuilt.md", "repro.json")
    }
    pre = {
        "frozen_files_check": F.check(Cfg.frozen),
        "n_frozen_files": len(json.loads(F.MANIFEST.read_text())),
        "g4_reference_sha256": pins()["g4_reference"]["sha256"],
        "unstamped_external_inputs_sha256": external_inputs(),
        "unstamped_external_inputs_note": "raw Visium h5 (outside the repository) and Xenium h5ad / ground-truth tables: read by the code, not pinnable by a namespace-relative stamp input",
    }
    return {
        "stage": STAGE,
        "run_id": run_id,
        "produced_by": PRODUCED_BY,
        "descriptive_only": True,
        "post_hoc": True,
        "git_commit": _git_head(),
        "config": config_block(),
        "preconditions": pre,
        **_clean(out),
        "repro": _clean(repro),
        "tables": {"variants": records(variants), "calls_breast": records(bcalls), "calls_xenium": records(xcalls)},
        "tasks": markers,
        "sibling_files_sha256": sib,
    }


def external_inputs() -> dict:
    """sha256 of the raw inputs no stamp can pin (the stamp lists namespace-relative files of earlier stamped stages): the Visium h5 of every breast specimen, and for
    Xenium the CITE and GEX h5ad and the ground-truth proportion table of every region. Recorded in the artifact, not in the stamp.
    """
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io

    out = {}
    for s in RS.SPECIMENS:
        f = Path(os.environ.get("CITEGEIST_DATA_ROOT", "data")) / RS.R.PROTO["data_root"] / s / "outs" / "filtered_feature_bc_matrix.h5"
        out[str(f)] = sha(f)
    for r in range(N_XEN):
        for f in (
            Path(data_io.h5ad_path(r, "CITE")),
            Path(data_io.h5ad_path(r, "GEX")),
            RS.GT_PROP / f"Xenium_region_{r}_prop.csv",
        ):
            out[str(f)] = sha(f)
    return out


def check_markers(run) -> list:
    run = Path(run)
    bad = sorted(p.name for p in run.glob("gate_failure_*.json"))
    if bad:
        raise GateFailure("earlier_gate_failure", {"files": bad})
    markers = []
    for n in range(N_TASKS):
        f = run / "tasks" / f"task_{n:02d}.json"
        if not f.is_file():
            raise GateFailure("task_missing", {"task": n, "expected": str(f)})
        markers.append(json.loads(f.read_text()))
    commits = {m["code_commit"] for m in markers}
    if commits != {_git_head()}:
        raise GateFailure("code_commit_mixed", {"task_commits": sorted(commits), "summary_commit": _git_head()})
    return markers


def run_summary(run, run_id, module3_dir, nuclei_dir, stamp=None):
    run = Path(run)
    configure(module3_dir, nuclei_dir)
    frozen_check()
    markers = check_markers(run)
    out, repro, variants, bcalls, xcalls = summary(run)
    art = build_artifact(run, run_id, out, repro, variants, bcalls, xcalls, markers)
    (run / ARTIFACT).write_text(json.dumps(art, indent=1, default=str))
    (run / "_COMPLETE.json").write_text(
        json.dumps(
            {
                "path": str(run),
                "jobs": sorted({m["slurm_job_id"] for m in markers} | {os.environ.get("SLURM_JOB_ID") or "no-slurm"}),
                "sha": _git_head(),
                "completed": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "produced_by": PRODUCED_BY,
            },
            indent=1,
        )
    )
    return art


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["task", "summary"])
    ap.add_argument("n", nargs="?", type=int, help="array task id for `task` (0..16)")
    ap.add_argument("--run-id", required=True, help="the B19 id (v42_run_ids.yaml)")
    ap.add_argument(
        "--out-root", type=Path, required=True, help="Benchmarking/ora_check/runs/<sha12>/<B19 id>; never overwritten"
    )
    ap.add_argument(
        "--module3-dir",
        type=Path,
        required=True,
        help="the B1 dir (output/module3_rna_refined_7type_<ts>/v2b1), with its p0/",
    )
    ap.add_argument("--nuclei-dir", type=Path, required=True, help="the B2 dir (output/nuclei_rebuilt_7type_<ts>/v2b2)")
    ap.add_argument("--protocol", type=Path, default=HERE.parent / "midkine_remediation" / "protocol.json")
    a = ap.parse_args(argv)
    PF.guard_entry(
        STAGE, a.protocol, run_id=a.run_id, out_root=a.out_root, module3_dir=a.module3_dir, nuclei_dir=a.nuclei_dir
    )
    stamp = PF.stage_stamp_begin(
        STAGE, a.protocol, a.out_root / ARTIFACT, Path(__file__).resolve(), a.run_id
    )  # job-start stamp gate, before any compute
    try:
        if a.cmd == "task":
            if a.n is None:
                raise SystemExit("task needs an id")
            run_task(a.n, a.out_root, a.module3_dir, a.nuclei_dir)
        else:
            if (a.out_root / ARTIFACT).exists():
                raise SystemExit(f"refusing: {a.out_root / ARTIFACT} exists (a run dir is never overwritten)")
            run_summary(a.out_root, a.run_id, a.module3_dir, a.nuclei_dir)
            if stamp:
                stamp.finish(module3_dir=a.module3_dir, nuclei_dir=a.nuclei_dir)
    except (GateFailure, AssertionError) as e:  # an input assertion of the copied code is a gate failure too
        if not isinstance(e, GateFailure):
            e = GateFailure("assertion", str(e))
        a.out_root.mkdir(parents=True, exist_ok=True)
        (a.out_root / f"gate_failure_{e.tag}.json").write_text(
            json.dumps(
                {"stage": STAGE, "run_id": a.run_id, "cmd": a.cmd, "task": a.n, "tag": e.tag, "info": e.info},
                indent=1,
                default=str,
            )
        )
        print(str(e), file=sys.stderr, flush=True)
        return 2
    print("done", a.cmd, a.n, a.run_id, a.out_root, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
