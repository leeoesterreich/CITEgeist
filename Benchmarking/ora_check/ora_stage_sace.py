"""B20: ORA on the SACE layers of the rebuilt breast spots, with a Xenium calibration (stamped stage; post-hoc, descriptive). Promoted from ``run_ora_sace_rebuilt.py``.

Rule file: ``ora_sace_rebuilt_DECISION_RULE.md`` (unchanged). Part 1 (protocol: specificity genes of the SACE layers, real and three nulls), Part 2 (Test A: the
delta of the layer against the uniform-share layer, real, null, w25/50/75/100), Part 3 (Xenium calibration, XP and XA), each under 4 mappings x 2 calls. A label
(PROPORTION_SENSITIVE / DICTIONARY_DRIVEN, the readout verdicts) is a VALUE, never a gate. ``run_sace_layers.py`` and ``run_ora.py`` run UNCHANGED (subprocess).

Unchanged from ``run_ora_sace_rebuilt.py`` (exploratory run ``results_sace_rebuilt_20261103_1115``): the pure helpers (``proportion_label``, ``gate_format``,
``compare_protocol``, ``discordance``, ``top_terms``, ``lineage_lean``, ``protocol_rows``, ``testA_rows``, ``x_summary``, ...) are IMPORTED from it, so the
analysis is the same code. ``breast``, ``xenium``, ``second_pass``, ``check_inputs_g4``, ``summary``, ``render_md`` are copies with only the changes listed in
``tests/test_ora_stage_sace.py::ALLOWED_DIFF``.

What changed (inputs and I/O only, plus one NAMED post-hoc block; no design change):
  * inputs come from the CLI (``--module3-dir`` = B1 dir with ``p0/``): v41 arm = the stamped p0, v2 arm = the stamped v2b1 table;
  * Xenium: the exploratory run read SACE layers of five arms written by ``run_ora_sense xenium-sace`` from ``gexgate_20261031_1832/props/v2_full`` and
    ``output_marker_share_p138_7type`` (v41/exploratory outputs, forbidden tokens, not hash-pinned). The stage REGENERATES the five proportion tables from the
    Xenium data (oracle, oracle_perm, uniform from the ground truth; v41 = ``ms_share``; v2 = ``refine_proportions`` full fit, c 0.1, n_passes 3, multinomial) and
    runs the same ``sace_props_ablation.run_region(..., "ms_v41", cell_mode="no_round")`` on them. The exploratory gate ``gex_gate.json`` measured these two
    tables against the files to 1.1e-16 (max abs); the fidelity test compares the resulting calls with the exploratory ones;
  * the legacy gate (task 17) is not part of the stage (its inputs include the v3 stand-in, which the chain refuses; the fidelity test replaces it);
  * G4 compares p0 with the plain marker-share table of the earlier arm as a fail-closed precondition on 12 files pinned in ``ora_stage_pins.json``;
  * the comparison with the OLD spots (``old_spots_results_v2`` from ``results_v2_20261101_1501/readout.json``, a v41-era result) is dropped;
  * NEW named block ``post_hoc_ungated_no_abundance_rule`` (``ungated_layer_ora``): Part 1 real_marker without the 2 percent abundance rule, see its docstring;
  * gates raise ``GateFailure`` (exit 2, ``gate_failure_<tag>.json``, no artifact, no stamp); the stamp is written by ``summary`` after the 17 task markers.

Usage (the wrapper ``.commands/submit-v2-ora.sh`` fills every path from ``chain_preflight.py``):
  ora_stage_sace.py task N --run-id ID --out-root DIR --module3-dir DIR --nuclei-dir DIR   N = 0..11 breast, 12..16 Xenium regions 0..4
  ora_stage_sace.py summary                                                                 checks the 17 markers, writes ``ora_sace_results.json`` and its stamp
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
from ora_check import ora_stage_b as SB  # noqa: E402
from ora_check import ora_v2_lib as V  # noqa: E402
from ora_check import (  # noqa: E402  (exploratory runner: the unchanged pure helpers are imported from it)
    run_ora_sace_rebuilt as SR,
)
from ora_check import run_ora_sense as RS  # noqa: E402

STAGE = "B20"
ARTIFACT = "ora_sace_results.json"
PRODUCED_BY = "Benchmarking/ora_check/ora_stage_sace.py"
RULE_FILE = HERE / "ora_sace_rebuilt_DECISION_RULE.md"
N_BREAST, N_XEN = SR.N_BREAST, SR.N_XEN
N_TASKS = N_BREAST + N_XEN
SPECIMENS, P_VARIANTS, BLENDS, ARMS = SR.SPECIMENS, SR.P_VARIANTS, SR.BLENDS, SR.ARMS
MIN_ABUNDANCE, MAPPINGS, VARIANTS, PRIMARY = SR.MIN_ABUNDANCE, SR.MAPPINGS, SR.VARIANTS, SR.PRIMARY
A_CONDS_FULL, LABEL_GAP, SP = SR.A_CONDS_FULL, SR.LABEL_GAP, SR.SP
proportion_label, gate_format, compare_protocol, discordance = (
    SR.proportion_label,
    SR.gate_format,
    SR.compare_protocol,
    SR.discordance,
)
top_terms, lineage_lean, tkey, protocol_rows, testA_rows = (
    SR.top_terms,
    SR.lineage_lean,
    SR.tkey,
    SR.protocol_rows,
    SR.testA_rows,
)
load_calls, x_summary, tab, lean_lines, to_bool, rec = (
    SR.load_calls,
    SR.x_summary,
    SR.tab,
    SR.lean_lines,
    SR.to_bool,
    SR.rec,
)
run_sace, run_ora_unchanged, write_override, libs_all = (
    SR.run_sace,
    SR.run_ora_unchanged,
    SR.write_override,
    SR.libs_all,
)
sha, write_new, write_tables, read_tables, check_vs_stored = (
    SB.sha,
    SB.write_new,
    SB.write_tables,
    SB.read_tables,
    SB.check_vs_stored,
)
GateFailure, fail, Cfg, configure, frozen_check, g4_reference = (
    SB.GateFailure,
    SB.fail,
    SB.Cfg,
    SB.configure,
    SB.frozen_check,
    SB.g4_reference,
)
_clean, records, _git_head = SB._clean, SB.records, SB._git_head


# ============================================================ inputs
def breast_inputs(sample: str) -> dict:
    return dict(
        v41=Cfg.m3 / "p0" / sample / f"{sample}_cell_prop_global_results.csv",
        v2=Cfg.m3 / f"{sample}_cell_prop_finetuned_results.csv",
    )


def check_inputs_g4(sample: str, inp: dict, types: list) -> dict:
    """G4: column order, equal barcode sets, row sums 1, p0 equals the official marker share on shared barcodes."""
    p0, v2 = (pd.read_csv(inp[k], index_col=0) for k in ARMS)
    out = {}
    for k, d in (("v41", p0), ("v2", v2)):
        assert list(d.columns) == types, f"{k}: columns {list(d.columns)} != {types}"
        assert float(np.abs(d.sum(axis=1) - 1).max()) < 1e-6, f"{k}: row sums"
    assert set(p0.index) == set(v2.index) and not p0.index.has_duplicates, "barcode sets differ"
    ref = pd.read_csv(g4_reference(sample), index_col=0)
    sh = p0.index.intersection(ref.index)
    diff = float((p0.loc[sh] - ref.loc[sh][types]).abs().max().max())
    assert diff == 0.0, f"p0 differs from module3_marker_share_7type: {diff}"
    out.update(n_spots=int(len(p0)), n_shared_with_official=int(len(sh)), p0_vs_official_max_abs=diff)
    return out


# ============================================================ NEW named post-hoc block
def ungated_layer_ora(Pt: pd.DataFrame) -> dict:
    """POST HOC (named output; not in the rule file, computed afterwards in the exploratory run from ``variants.csv`` of ``results_sace_rebuilt_20261103_1115``).

    Layer ORA (Part 1, condition ``real_marker``, the call of this readout's mapping and call type) WITHOUT the 2 percent abundance rule: the pass count of each arm
    is the sum of ``passed`` over the 84 (specimen, type) pairs, ignoring ``evaluable``; the discordance is ``discordance()`` on the same ungated calls
    (v41_only = v41 passes and v2 fails, v2_only = the reverse, exact two-sided binomial at 0.5). For ``frozen|current`` the exploratory run gives v41 37, v2 34, 6 v41-only,
    3 v2-only, P 0.5078 (the same numbers for ``L1|current``). ``Pt`` is ``RL.variant_table`` of the Part 1 calls (columns arm, condition, unit, type, passed).
    """
    s = {
        a: Pt[(Pt["arm"] == a) & (Pt["condition"] == "real_marker")].set_index(["unit", "type"])["passed"].astype(bool)
        for a in ARMS
    }
    return {
        "label": "POST_HOC_NO_ABUNDANCE_RULE",
        "condition": "real_marker",
        "rule": "passed ignoring the 2 percent abundance rule (evaluable not applied)",
        "passes_of_84": {a: int(s[a].sum()) for a in ARMS},
        "discordance": discordance(s["v41"], s["v2"]),
    }


# ============================================================ SACE + ORA second pass (copy; gates raise GateFailure)


def second_pass(root: Path, sample: str, idx: int, conds_a: dict, variants: list, run: Path, tag: str):
    """Both arms: protocol rows + Test A rows + tables; asserts the second pass equals run_ora.py (G3). Returns (P rows, A rows, tables, G3 info)."""
    mapping = RS.load_mapping("term_mapping_7type.yaml")
    from ora_check.run_sace_layers import load_v41_profile

    markers = L.marker_gene_set(load_v41_profile())
    prow, arow, tabs, g3 = [], [], {}, {}
    libs = libs_all()
    for arm in ARMS:
        d = root / arm / "sace_layers" / sample
        r, t = protocol_rows(d, arm, sample, mapping, markers, variants)
        prow += r
        tabs.update(t)
        r, t = testA_rows(d, arm, sample, mapping, markers, conds_a)
        arow += r
        tabs.update(t)
        o = pd.read_csv(root / arm / f"ora_calls{tag}.csv", keep_default_na=False, na_values=[""])
        mine = gate_format(pd.DataFrame([x for x in prow if x["arm"] == arm]), libs)
        bad = compare_protocol(mine, o, libs)
        g3[arm] = {"compared": int(len(o)), "n_mismatch": len(bad), "mismatches": bad[:20]}
        if bad or not len(o):
            fail(run, f"G3_{tag}_{arm}", g3[arm])
    return pd.DataFrame(prow), pd.DataFrame(arow), tabs, g3


# ============================================================ breast task


def breast(idx: int, run: Path):
    sample = SPECIMENS[idx]
    frozen_check()
    inp = breast_inputs(sample)
    for k, p in inp.items():
        assert p.is_file(), f"missing input {k}: {p}"
    types = list(pd.read_csv(inp["v41"], index_col=0).columns)
    g4 = check_inputs_g4(sample, inp, types)
    write_new(
        run / "inputs" / f"breast_{idx:02d}.json",
        json.dumps(
            {k: {"path": str(p), "sha256": sha(p)} for k, p in inp.items()}
            | {"g4": g4, "g4_reference_sha256": sha(g4_reference(sample))},
            indent=1,
        ),
    )
    root = run / "rebuilt"
    run_sace(idx, root, inp, P_VARIANTS + BLENDS)
    meta = {}
    for arm in ARMS:
        m = json.loads((root / arm / "sace_layers" / sample / "run_meta.json").read_text())
        meta[arm] = {k: m[k] for k in ("n_spots", "n_arm_table_spots", "n_model_spots", "n_arm_spots_dropped")}
        if m["n_arm_spots_dropped"] != 0 or m["n_spots"] != m["n_model_spots"] or m["n_spots"] != g4["n_spots"]:
            fail(run, f"breast_{idx:02d}_spots_{arm}", {"meta": m, "g4": g4})
    ov = run / "overrides" / f"spec_{idx:02d}.json"
    write_override(ov, sample)
    for arm in ARMS:
        run_ora_unchanged(root / arm, ov, f"_{idx:02d}")
    prow, arow, tabs, g3 = second_pass(root, sample, idx, A_CONDS_FULL, P_VARIANTS, run, f"_{idx:02d}")
    write_new(run / "calls" / f"P_breast_{idx:02d}.csv", prow.to_csv(index=False))
    write_new(run / "calls" / f"A_breast_{idx:02d}.csv", arow.to_csv(index=False))
    write_tables(run / "tables" / f"breast_{idx:02d}.parquet", tabs)
    write_new(
        run / "calls" / f"breast_{idx:02d}_meta.json",
        json.dumps({"sample": sample, "spots": meta, "g4": g4, "G3": g3}, indent=2),
    )


# ============================================================ Xenium task: the five arms' proportions are regenerated, then the same SACE and ORA code as before
def xenium_arm_props(r: int, props_root: Path):
    """Write the five proportion tables of region ``r`` (``<props_root>/<arm>/Xenium_region_<r>/predictions/baseline.csv.gz``) as ``run_ora_sense.xenium_sace`` and
    ``run_gex_gate.region_job`` did, but from the Xenium data instead of reading earlier outputs: oracle / oracle_perm / uniform from the ground truth, v41 =
    ``ms_share`` (the marker share, == output_marker_share_p138 to 1.1e-16), v2 = the full-gene ``refine_proportions`` fit (== gexgate v2_full to 1.1e-16).
    """
    from Benchmarking.xenium_benchmarking.rna_refinement import data_io
    from Benchmarking.xenium_benchmarking.rna_refinement import reproduce_spike as rs
    from Benchmarking.xenium_benchmarking.rna_refinement import run_gex_gate as gg
    from CITEgeist.model.deconvolution.rna_refinement import refine_proportions

    types, data = data_io.load_xenium(7)
    d = data[r]
    gt = RS._gt(r)
    p0 = np.asarray(rs.ms_share(d["X"], d["names"], rs.PROF[7]))
    p_full, _ = refine_proportions(d["Y"], d["nc"], p0, c=0.1, n_passes=3, model="multinomial")
    perm, _ = L.spot_permuted(gt, RS.XEN_SEED + r)
    uni = pd.DataFrame(1.0 / len(RS.XTYPES), index=d["obs"], columns=RS.XTYPES)
    arms = {
        "oracle": gt,
        "oracle_perm": perm,
        "uniform": uni,
        "v2": pd.DataFrame(p_full, index=d["obs"], columns=types),
        "v41": pd.DataFrame(p0, index=d["obs"], columns=types),
    }
    for arm, df in arms.items():
        gg._write_props(df, props_root, arm, r)
    return {arm: props_root / arm for arm in arms}


def _sace_props_ablation():
    """The committed module (from the code tree) with its untracked-data root moved to the namespace root."""
    sys.path.insert(0, str(RS.XB / "spots_baseline"))
    import sace_props_ablation as sp

    sp.DATA = Cfg.ns / SB.H5AD_REL
    return sp


def data_paths() -> dict:
    return {**SB.data_paths(), "sace_props_ablation_DATA": Path(_sace_props_ablation().DATA)}


def xenium_sace(r: int, run: Path):
    sp = _sace_props_ablation()

    dirs = xenium_arm_props(r, run / "xenium_props")
    for arm, pdir in dirs.items():
        sp.run_region("marker_share", r, run / "xenium_sace" / "sace" / arm, str(pdir), "ms_v41", cell_mode="no_round")
        print("sace done", arm, r, flush=True)
    return dirs


def xenium(r: int, run: Path):
    frozen_check()
    dirs = xenium_sace(r, run)
    mapping = RS.load_mapping("term_mapping_xenium7.yaml")
    markers = RS._xen_markers()
    shas = {}
    for arm in ("oracle", "oracle_perm", "uniform", "v2", "v41"):
        f = dirs[arm] / f"Xenium_region_{r}" / "predictions" / "baseline.csv.gz"
        shas[str(f.relative_to(run))] = sha(f)
    write_new(run / "inputs" / f"xenium_{r}.json", json.dumps(shas, indent=1))
    lay = {a: RS._xen_layer_means(run / "xenium_sace", a, r) for a in ("oracle", "oracle_perm", "uniform", "v2", "v41")}
    panel = list(lay["oracle"][0].columns)
    tb_frozen = RS._testable(panel, mapping)
    tabs, xp, xa = {}, [], []

    def add(rows, test, arm, cond, genes, bg, props):
        for t in RS.XTYPES:
            row, tb = rec().row(genes[t], bg, t, mapping)
            k = tkey(test, arm, cond, f"region_{r}", t)
            tabs[k] = tb
            rows.append(
                dict(
                    test=test,
                    dataset="xenium",
                    arm=arm,
                    condition=cond,
                    unit=f"region_{r}",
                    type=t,
                    evaluable=bool(S.evaluable_share(props[t])),
                    testable=bool(tb_frozen[t] >= S.MIN_PANEL_GENES),
                    n_panel_genes=tb_frozen[t],
                    table_key=k,
                    **row,
                )
            )

    for arm, src, cond in (
        ("oracle", "oracle", "real"),
        ("oracle", "oracle_perm", "null"),
        ("uniform", "uniform", "real"),
        ("v2", "v2", "real"),
        ("v41", "v41", "real"),
    ):
        lm, props = lay[src]
        layers = {t: lm.loc[[t]] for t in lm.index}
        genes = L.specificity_genes(
            layers, markers, n=SP["top_n"], min_mean=SP["min_mean"], pseudocount=SP["pseudocount"]
        )
        bg = L.expressed_background(layers, SP["min_mean"])
        add(xp, "XP", arm, cond, genes, bg, props)
    uni = lay["uniform"][0]
    oracle_lay = lay["oracle"][0]
    for arm, src, cond in (
        ("oracle", "oracle", "real"),
        ("oracle", "oracle_perm", "null"),
        ("v2", "v2", "real"),
        ("v41", "v41", "real"),
    ):
        lm, props = lay[src]
        real_lay = lm if cond == "real" else oracle_lay
        expressed = real_lay.sum(axis=0) >= S.MIN_MEAN
        cand = [g for g in real_lay.columns[expressed.to_numpy()] if g not in markers]
        bg = list(real_lay.columns[expressed.to_numpy()])
        genes = RS.gene_lists_from_scores(S.delta_table(lm, uni).T, cand)
        add(xa, "XA", arm, cond, genes, bg, props)
    xa_df = pd.DataFrame(xa)
    stored_csv = Cfg.frozen / "results_sense_20261101_1550" / "calls" / f"A_xenium_{r}.csv"
    chk = check_vs_stored(xa_df, stored_csv, f"XA region {r} vs stored A_xenium")
    chk["stored_file_sha256"] = sha(stored_csv)
    all_testable = all(v >= S.MIN_PANEL_GENES for v in tb_frozen.values())
    write_new(run / "calls" / f"XP_{r}.csv", pd.DataFrame(xp).to_csv(index=False))
    write_new(run / "calls" / f"XA_{r}.csv", xa_df.to_csv(index=False))
    write_tables(run / "tables" / f"xenium_{r}.parquet", tabs)
    write_new(
        run / "calls" / f"xenium_{r}_meta.json",
        json.dumps({"region": r, "panel_genes": panel, "repro": chk, "all_testable_frozen": all_testable}, indent=2),
    )
    if chk["n_mismatch"] or chk["compared"] != 28 or not all_testable:
        fail(run, f"xenium_{r}_repro", {"repro": chk, "all_testable_frozen": all_testable})


def summary(run: Path):
    frozen_check()
    ov = M.load_overrides()
    mp_b, mp_x = RS.load_mapping("term_mapping_7type.yaml"), RS.load_mapping("term_mapping_xenium7.yaml")
    P = load_calls(run, "P_breast", "[0-9][0-9]")
    Ab = load_calls(run, "A_breast", "[0-9][0-9]")
    XP = load_calls(run, "XP", "[0-9]")
    XA = load_calls(run, "XA", "[0-9]")
    if (
        P["unit"].nunique() != S.N_SPECIMENS
        or Ab["unit"].nunique() != S.N_SPECIMENS
        or XP["unit"].nunique() != S.N_REGIONS
        or XA["unit"].nunique() != S.N_REGIONS
    ):
        fail(
            run,
            "summary_units",
            {
                "P": int(P["unit"].nunique()),
                "A": int(Ab["unit"].nunique()),
                "XP": int(XP["unit"].nunique()),
                "XA": int(XA["unit"].nunique()),
            },
        )
    tables = {}
    for f in sorted((run / "tables").glob("*.parquet")):
        tables.update(read_tables(f))
    panel = {
        f"region_{r}": json.loads((run / "calls" / f"xenium_{r}_meta.json").read_text())["panel_genes"]
        for r in range(S.N_REGIONS)
    }
    gmt = RS.libs()[PRIMARY]
    res, allv = {}, []
    g3 = {}
    for mname, lv in MAPPINGS.items():
        fn = M.corrected_lineage_fn(ov, lv)
        lb, lx = RL.LineageCache(fn, mp_b), RL.LineageCache(fn, mp_x)
        tcount = {u: RL.testable_counts([g.upper() for g in panel[u]], gmt, mp_x, lx) for u in panel}
        for var in VARIANTS:
            collapsed = var == "collapsed"
            Pt = RL.variant_table(P, tables, mp_b, lb, var)
            At = RL.variant_table(Ab, tables, mp_b, lb, var)
            XPt = RL.variant_table(XP, tables, mp_x, lx, var, testable=tcount)
            XAt = RL.variant_table(XA, tables, mp_x, lx, var, testable=tcount)
            if mname == "frozen" and var == "current":  # G3: table path equals the in-job calls
                tst = np.array([tcount[u][t] >= S.MIN_PANEL_GENES for u, t in zip(XP["unit"], XP["type"])])
                tsa = np.array([tcount[u][t] >= S.MIN_PANEL_GENES for u, t in zip(XA["unit"], XA["type"])])
                g3 = {
                    "P": int((Pt["passed"].to_numpy() != P["passed"].to_numpy()).sum()),
                    "A": int((At["passed"].to_numpy() != Ab["passed"].to_numpy()).sum()),
                    "XP": int(((XPt["passed"].to_numpy() != XP["passed"].to_numpy()) & tst).sum()),
                    "XA": int(((XAt["passed"].to_numpy() != XA["passed"].to_numpy()) & tsa).sum()),
                    "nlp_max_abs": float(
                        max(
                            np.abs(Pt["nlp"].to_numpy(float) - P["nlp"].to_numpy(float)).max(),
                            np.abs(At["nlp"].to_numpy(float) - Ab["nlp"].to_numpy(float)).max(),
                        )
                    ),
                }
                if any(g3[k] for k in ("P", "A", "XP", "XA")) or g3["nlp_max_abs"] > 1e-9:
                    fail(run, "G3_summary", g3)
            Pt["eff"] = Pt["passed"] & Pt["evaluable"]  # protocol call: a non-evaluable pair does not pass
            # ---- Part 1 (protocol)
            calls_v2 = {
                a: pd.DataFrame(
                    {"specimen": g["unit"], "variant": g["condition"], "type": g["type"], "passed": g["eff"]}
                )
                for a, g in Pt.groupby("arm")
            }
            ro = V.readout(calls_v2)
            pv = {a: {v: int(g[g["condition"] == v]["eff"].sum()) for v in P_VARIANTS} for a, g in Pt.groupby("arm")}
            pv_ungated = {
                a: {v: int(g[g["condition"] == v]["passed"].sum()) for v in P_VARIANTS} for a, g in Pt.groupby("arm")
            }
            pt = {
                a: {v: {t: int(h["eff"].sum()) for t, h in g[g["condition"] == v].groupby("type")} for v in P_VARIANTS}
                for a, g in Pt.groupby("arm")
            }
            ev_by_type = {
                a: {t: int(h["evaluable"].sum()) for t, h in g[g["condition"] == "real_marker"].groupby("type")}
                for a, g in Pt.groupby("arm")
            }
            labelP = {
                a: {
                    "gap_vs_spot_permuted_marker": pv[a]["real_marker"] - pv[a]["spot_permuted_marker"],
                    "label": proportion_label(pv[a]["real_marker"], pv[a]["spot_permuted_marker"]),
                }
                for a in ARMS
            }
            series = {
                a: Pt[(Pt["arm"] == a) & (Pt["condition"] == "real_marker")].set_index(["unit", "type"])["eff"]
                for a in ARMS
            }
            disc_p = discordance(series["v41"], series["v2"])
            pairs_b = [k for k in series["v41"].index if series["v41"][k] and not series["v2"][k]]
            pairs_c = [k for k in series["v41"].index if not series["v41"][k] and series["v2"][k]]
            kfn = lambda arm, u, t: tkey("P", arm, "real_marker", u, t)  # noqa: E731
            lean_p = {
                "v41_only": lineage_lean(tables, pairs_b, kfn, lb, collapsed),
                "v2_only": lineage_lean(tables, pairs_c, kfn, lb, collapsed),
                "all_pairs": lineage_lean(tables, list(series["v41"].index), kfn, lb, collapsed),
            }
            # ---- Part 2 (Test A)
            ra = S.readout(At, XAt)
            pa = {
                a: {c: int(g[g["condition"] == c]["passed"].sum()) for c in A_CONDS_FULL} for a, g in At.groupby("arm")
            }
            labelA = {
                a: {
                    "gap_vs_null": pa[a]["real"] - pa[a]["null"],
                    "label": proportion_label(pa[a]["real"], pa[a]["null"]),
                }
                for a in ARMS
            }
            sa = {
                a: At[(At["arm"] == a) & (At["condition"] == "real")].set_index(["unit", "type"])["passed"]
                for a in ARMS
            }
            disc_a = discordance(sa["v41"], sa["v2"])
            pb_a = [k for k in sa["v41"].index if sa["v41"][k] and not sa["v2"][k]]
            pc_a = [k for k in sa["v41"].index if not sa["v41"][k] and sa["v2"][k]]
            kfa = lambda arm, u, t: tkey("A", arm, "real", u, t)  # noqa: E731
            lean_a = {
                "v41_only": lineage_lean(tables, pb_a, kfa, lb, collapsed),
                "v2_only": lineage_lean(tables, pc_a, kfa, lb, collapsed),
                "all_pairs": lineage_lean(tables, list(sa["v41"].index), kfa, lb, collapsed),
            }
            v3p, v3a = S.v3_status(XPt), S.v3_status(XAt)
            res[f"{mname}|{var}"] = {
                "P": {
                    "post_hoc_ungated_no_abundance_rule": ungated_layer_ora(Pt),
                    "passes_of_84": pv,
                    "passes_ungated_of_84": pv_ungated,
                    "by_type_by_variant": pt,
                    "evaluable_real_marker_by_type": ev_by_type,
                    "label": labelP,
                    "ora_v2_readout": ro,
                    "discordance_real_marker": disc_p,
                    "lineage_lean": lean_p,
                },
                "A": {
                    "passes_of_84": pa,
                    "label": labelA,
                    "sense_readout": ra,
                    "discordance_real": disc_a,
                    "lineage_lean": lean_a,
                },
                "XP": {"v3": v3p, "summary": x_summary(XPt), "testable_counts": tcount},
                "XA": {"v3": v3a, "summary": x_summary(XAt)},
            }
            keep = ["dataset", "test", "arm", "condition", "unit", "type", "evaluable", "passed", "best_rank", "nlp"]
            allv.append(
                pd.concat(
                    [Pt.assign(passed_protocol=Pt["eff"])[keep + ["passed_protocol"]], At[keep], XPt[keep], XAt[keep]]
                ).assign(mapping=mname, variant=var)
            )
    variants = pd.concat(allv, ignore_index=True)
    variants.to_csv(run / "variants.csv", index=False)
    meta = {f: json.loads(f.read_text()) for f in sorted((run / "calls").glob("breast_[0-9][0-9]_meta.json"))}
    repro = {
        "G3_summary": g3,
        "G3_runora_vs_second_pass": {m["sample"]: m["G3"] for m in meta.values()},
        "G4": {m["sample"]: m["g4"] for m in meta.values()},
        "spots": {m["sample"]: m["spots"] for m in meta.values()},
        "xenium_XA_repro": {
            r: json.loads((run / "calls" / f"xenium_{r}_meta.json").read_text())["repro"] for r in range(S.N_REGIONS)
        },
    }
    out = {"results": res, "frozen_files_check": F.check(Cfg.frozen)}
    (run / "readout_sace.json").write_text(
        json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    )
    (run / "repro.json").write_text(json.dumps(repro, indent=1, default=str))
    (run / "readout_sace.md").write_text(render_md(out))
    print(
        json.dumps(
            {
                k: {
                    "P_labels": {a: v["P"]["label"][a]["label"] for a in ARMS},
                    "A_labels": {a: v["A"]["label"][a]["label"] for a in ARMS},
                }
                for k, v in res.items()
            }
        ),
        flush=True,
    )
    return out, repro, variants, P, Ab, XP, XA


def render_md(out):
    res = out["results"]
    L_ = [
        "# ORA on the SACE layers, rebuilt breast spots, plus Xenium calibration (primary library CellMarker_2024)",
        "",
        "Exploratory, post-hoc. Rule: `ora_sace_rebuilt_DECISION_RULE.md`. Plausibility of gene-to-type links, not accuracy. Nothing here chooses a number.",
        "",
        "The comparison with the old spots (a v41-era result) is not part of the stage.",
        "",
    ]
    for key, v in res.items():
        L_ += [f"## Readout `{key}` (mapping | call)", ""]
        p = v["P"]
        L_ += ["### Part 1: protocol (specificity genes), passes of 84 (evaluable = share >= 0.02)", ""]
        L_ += tab(
            [[var, p["passes_of_84"]["v41"][var], p["passes_of_84"]["v2"][var]] for var in P_VARIANTS],
            ["variant", "v41", "v2"],
        ) + [""]
        L_ += tab(
            [[a, p["label"][a]["gap_vs_spot_permuted_marker"], p["label"][a]["label"]] for a in ARMS],
            ["arm", "real_marker - spot_permuted_marker", "label"],
        ) + [""]
        r = p["ora_v2_readout"]
        L_ += [
            f"`ora_v2_lib.readout`: verdict {r['verdict']}, totals {json.dumps(r.get('total'))}, delta {r.get('delta')}, loss v41-v2 {json.dumps(r.get('loss'))}, not sensitive to proportions {json.dumps(r.get('not_sensitive_to_proportions'))}.",
            "",
        ]
        types = sorted(p["by_type_by_variant"]["v41"]["real_marker"])
        L_ += ["Passes of 12 by type, v41 / v2:", ""]
        L_ += tab(
            [
                [t]
                + [
                    f"{p['by_type_by_variant']['v41'][var].get(t, 0)} / {p['by_type_by_variant']['v2'][var].get(t, 0)}"
                    for var in P_VARIANTS
                ]
                for t in types
            ],
            ["type"] + P_VARIANTS,
        ) + [""]
        L_ += [
            f"Evaluable (of 12) real_marker by type: {json.dumps(p['evaluable_real_marker_by_type'])}. Ungated passes of 84 (real_marker): {json.dumps({a: p['passes_ungated_of_84'][a]['real_marker'] for a in ARMS})}.",
            "",
        ]
        d = p["discordance_real_marker"]
        L_ += [f"Discordance real_marker: overall {json.dumps(d['overall'])}.", ""]
        L_ += tab(
            [
                [t, x["passes_v41"], x["passes_v2"], x["v41_only"], x["v2_only"], f"{x['p_exact']:.4f}"]
                for t, x in d["by_type"].items()
            ],
            ["type", "v41 passes", "v2 passes", "v41 only", "v2 only", "exact P"],
        ) + [""]
        L_ += ["Lineage of the top-5 terms (discordant pairs):"] + lean_lines(p["lineage_lean"]) + [""]
        a = v["A"]
        L_ += ["### Part 2: Test A (delta of layers), passes of 84", ""]
        L_ += tab(
            [[c, a["passes_of_84"]["v41"][c], a["passes_of_84"]["v2"][c]] for c in A_CONDS_FULL],
            ["condition", "v41", "v2"],
        ) + [""]
        L_ += tab(
            [[x, a["label"][x]["gap_vs_null"], a["label"][x]["label"]] for x in ARMS], ["arm", "real - null", "label"]
        ) + [""]
        s = a["sense_readout"]
        L_ += [
            f"`ora_sense_lib.readout`: verdict {s['verdict']}, totals {json.dumps(s.get('total'))}, delta {s.get('delta')}, loss {json.dumps(s.get('loss'))}, V1 {json.dumps(s.get('V1'))}, V2 {json.dumps(s.get('V2'))}.",
            "",
        ]
        da = a["discordance_real"]
        L_ += [f"Discordance real: overall {json.dumps(da['overall'])}.", ""]
        L_ += tab(
            [
                [t, x["passes_v41"], x["passes_v2"], x["v41_only"], x["v2_only"], f"{x['p_exact']:.4f}"]
                for t, x in da["by_type"].items()
            ],
            ["type", "v41 passes", "v2 passes", "v41 only", "v2 only", "exact P"],
        ) + [""]
        L_ += ["Lineage of the top-5 terms (discordant pairs):"] + lean_lines(a["lineage_lean"]) + [""]
        L_ += ["### Part 3: Xenium calibration (35 calls)", ""]
        for nm in ("XP", "XA"):
            x = v[nm]
            L_ += [
                f"{nm} ({'protocol-style specificity' if nm == 'XP' else 'Test A delta'}): V3 = {json.dumps(x['v3'])}",
                "",
                f"passes of 35: {json.dumps(x['summary']['passes'])}; mean nlp: {json.dumps(x['summary']['mean_nlp'])}",
                "",
            ]
            arms = list(x["summary"]["by_type"])
            xt = sorted(x["summary"]["by_type"][arms[0]])
            L_ += tab([[t] + [x["summary"]["by_type"][k].get(t, 0) for k in arms] for t in xt], ["type"] + arms) + [""]
    return "\n".join(L_) + "\n"


# ------------------------------------------------------------------ stage driver: markers, artifact, stamp
def run_task(n, run, module3_dir, nuclei_dir):
    """One array task: ``n`` = 0..11 breast (SACE + ORA), 12..16 Xenium. Writes its outputs and then its marker ``tasks/task_NN.json`` (never overwritten)."""
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
        "label_gap": LABEL_GAP,
        "label_rule": "PROPORTION_SENSITIVE iff passes(real) - passes(null) >= LABEL_GAP of 84, else DICTIONARY_DRIVEN (a value, not a gate)",
        "min_layer_abundance": MIN_ABUNDANCE,
        "p_variants": list(P_VARIANTS),
        "blends": list(BLENDS),
        "a_conditions": A_CONDS_FULL,
        "seeds": {
            "breast_base": RS.BASE_SEED,
            "xenium_base": RS.XEN_SEED,
            "sace": "ora_protocol.json seeds.base + specimen index (run_sace_layers.py, unchanged)",
        },
        "top_n": SR.SP["top_n"],
        "min_mean": SR.SP["min_mean"],
        "pseudocount": SR.SP["pseudocount"],
        "g3_tol_nlp": 1e-9,
        "primary_library": PRIMARY,
        "mappings": list(MAPPINGS),
        "calls": list(VARIANTS),
        "specimens": list(SPECIMENS),
        "xenium_types": list(RS.XTYPES),
        "post_hoc_named_block": "results[<mapping>|<call>].P.post_hoc_ungated_no_abundance_rule (see ungated_layer_ora)",
        "removed_from_exploratory": ["legacy gate (task 17)", "old_spots_results_v2 (a v41-era result)"],
        "xenium_inputs": "regenerated from the Xenium data and ground truth in the stage (xenium_arm_props), not read from earlier outputs",
    }


def build_artifact(run, run_id, out, repro, variants, P, Ab, XP, XA, markers) -> dict:
    sib = {name: sha(run / name) for name in ("variants.csv", "readout_sace.json", "readout_sace.md", "repro.json")}
    pre = {
        "frozen_files_check": F.check(Cfg.frozen),
        "n_frozen_files": len(json.loads(F.MANIFEST.read_text())),
        "g4_reference_sha256": SB.pins()["g4_reference"]["sha256"],
        "unstamped_external_inputs_sha256": SB.external_inputs(),
        "unstamped_external_inputs_note": "raw Visium h5 (outside the repository, read by run_sace_layers.load_model) and Xenium h5ad / ground-truth tables (data_io, sace_props_ablation): not pinnable by a namespace-relative stamp input",
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
        "tables": {
            "variants": records(variants),
            "calls_P_breast": records(P),
            "calls_A_breast": records(Ab),
            "calls_XP": records(XP),
            "calls_XA": records(XA),
        },
        "tasks": markers,
        "sibling_files_sha256": sib,
    }


def run_summary(run, run_id, module3_dir, nuclei_dir):
    run = Path(run)
    configure(module3_dir, nuclei_dir)
    frozen_check()
    markers = SB.check_markers(run)
    out, repro, variants, P, Ab, XP, XA = summary(run)
    art = build_artifact(run, run_id, out, repro, variants, P, Ab, XP, XA, markers)
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
    ap.add_argument("--run-id", required=True, help="the B20 id (v42_run_ids.yaml)")
    ap.add_argument(
        "--out-root", type=Path, required=True, help="Benchmarking/ora_check/runs/<sha12>/<B20 id>; never overwritten"
    )
    ap.add_argument(
        "--module3-dir",
        type=Path,
        required=True,
        help="the B1 dir (output/module3_rna_refined_7type_<ts>/v2b1), with its p0/",
    )
    ap.add_argument(
        "--nuclei-dir",
        type=Path,
        required=True,
        help="the B2 dir (output/nuclei_rebuilt_7type_<ts>/v2b2); checked by the guard, spot-set equality is gate G4's job",
    )
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
