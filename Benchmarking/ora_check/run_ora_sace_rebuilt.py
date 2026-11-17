"""ORA on the SACE layers, rebuilt breast spots, plus a Xenium calibration (rule: ora_sace_rebuilt_DECISION_RULE.md). CPU only (htc).

  task N RUN   N = 0..11 breast (SACE layers + ORA on rebuilt inputs), 12..16 Xenium regions 0..4, 17 = legacy gate (creates RUN)
  summary RUN  needs gate_legacy/G1_PASS.json; recomputes 8 mapping x call readouts, writes readout_sace.{json,md}, variants.csv, repro.json

Exploratory and post-hoc; not a v42 stage. run_sace_layers.py and run_ora.py run UNCHANGED (subprocess); this file adds a second ORA pass that keeps
the full primary-library tables (for collapsed calls, nlp, term lineage) and asserts it equals run_ora.py's calls.
"""
import argparse
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE.parent))
from ora_check import ora_audit_lib as A  # noqa: E402
from ora_check import ora_corrected_mapping as M  # noqa: E402
from ora_check import ora_frozen_manifest as F  # noqa: E402
from ora_check import ora_lib as L  # noqa: E402
from ora_check import ora_rebuilt_lib as RL  # noqa: E402
from ora_check import ora_sense_lib as S  # noqa: E402
from ora_check import ora_v2_lib as V  # noqa: E402
from ora_check import run_ora_rebuilt as RR  # noqa: E402
from ora_check import run_ora_sense as RS  # noqa: E402

V2_TS = "20261102T200000Z"
M3 = REPO / "output" / f"module3_rna_refined_7type_{V2_TS}" / "v2b1"
SENSE = HERE / "results_sense_20261101_1550"
OLDV2 = HERE / "results_v2_20261101_1501"
PROTO = RS.R.PROTO
SPECIMENS = PROTO["specimens"]
P_VARIANTS = list(PROTO["variants"])
BLENDS = ["blend25_marker", "blend50_marker", "blend75_marker"]
ARMS = ("v41", "v2")
N_BREAST, N_XEN, GATE_TASK = 12, 5, 17
LABEL_GAP = 5  # PROPORTION_SENSITIVE iff passes(real) - passes(null) >= 5 of 84
MIN_ABUNDANCE = PROTO["min_layer_abundance"]["threshold"]
MAPPINGS, VARIANTS = RR.MAPPINGS, RR.VARIANTS
PRIMARY = RS.PRIMARY
STROMAL = ("fibroblast", "endothelial")
OVERRIDE = json.loads((HERE / "ora_protocol_v41_7type.json").read_text())
SP = PROTO["specificity"]
_REC = None


def rec():
    global _REC
    if _REC is None:
        _REC = RR.Rec()
    return _REC


def to_bool(s):
    return s.map(lambda x: str(x).strip().lower() == "true") if s.dtype == object else s.astype(bool)


# ============================================================ pure helpers (unit tested)
def proportion_label(real: int, null: int, thr: int = LABEL_GAP) -> str:
    return "PROPORTION_SENSITIVE" if real - null >= thr else "DICTIONARY_DRIVEN"


def gate_format(mine: pd.DataFrame, libs) -> pd.DataFrame:
    """Second-pass protocol rows (ungated pass_*) -> run_ora.py's ora_calls format (a pass counts only for an evaluable pair)."""
    o = mine.rename(columns={"unit": "specimen", "condition": "variant"}).copy()
    for lib in libs:
        o[f"pass_{lib}"] = o[f"pass_{lib}"].astype(bool) & o["evaluable"].astype(bool)
    o["passed"] = o["passed"].astype(bool) & o["evaluable"].astype(bool)
    return o


def compare_protocol(a: pd.DataFrame, b: pd.DataFrame, libs, tol: float = 1e-9) -> list:
    """Mismatches between two run_ora-format tables on (specimen, variant, type): evaluable, passed, pass_<lib>, abundance."""
    key = ["specimen", "variant", "type"]
    cols = ["evaluable", "passed"] + [f"pass_{l}" for l in libs]
    m = a[key + ["abundance"] + cols].merge(b[key + ["abundance"] + cols], on=key, how="outer", suffixes=("_a", "_b"), indicator=True)
    bad = []
    for _, r in m[m["_merge"] != "both"].iterrows():
        bad.append(("unmatched", r["specimen"], r["variant"], r["type"], str(r["_merge"])))
    both = m[m["_merge"] == "both"]
    for c in cols:
        x, y = to_bool(both[f"{c}_a"]), to_bool(both[f"{c}_b"])
        for i in np.flatnonzero((x != y).to_numpy()):
            r = both.iloc[i]
            bad.append((c, r["specimen"], r["variant"], r["type"], bool(x.iloc[i]), bool(y.iloc[i])))
    d = (both["abundance_a"].astype(float) - both["abundance_b"].astype(float)).abs()
    for i in np.flatnonzero((d > tol).to_numpy()):
        r = both.iloc[i]
        bad.append(("abundance", r["specimen"], r["variant"], r["type"], float(d.iloc[i])))
    return bad


def discordance(v41: pd.Series, v2: pd.Series) -> dict:
    """v41, v2: bool Series indexed by (unit, type). b = v41 passes and v2 fails; c = v2 passes and v41 fails; exact two-sided binomial at 0.5."""
    d = pd.concat({"v41": v41.astype(bool), "v2": v2.astype(bool)}, axis=1)
    if d.isna().any().any():
        raise ValueError("pairs missing in one arm")

    def one(x):
        b = int((x["v41"] & ~x["v2"]).sum())
        c = int((~x["v41"] & x["v2"]).sum())
        p = 1.0 if b + c == 0 else float(binomtest(min(b, c), b + c, 0.5).pvalue)
        return {"n": int(len(x)), "both": int((x["v41"] & x["v2"]).sum()), "neither": int((~x["v41"] & ~x["v2"]).sum()),
                "v41_only": b, "v2_only": c, "p_exact": p, "passes_v41": int(x["v41"].sum()), "passes_v2": int(x["v2"].sum())}

    out = {"overall": one(d), "by_type": {}}
    for t, g in d.groupby(level=1):
        out["by_type"][t] = one(g)
    return out


def top_terms(tb: pd.DataFrame, collapsed: bool, k: int = 5) -> list:
    """Top-k significant terms (adjusted P, ties by term); collapsed: first k distinct parsed names (best-P term kept)."""
    terms = list(RL.significant(tb)["Term"])
    if not collapsed:
        return terms[:k]
    seen, kept = set(), []
    for t in terms:
        key = A.name_key(A.parse_term(t)[0])
        if key in seen:
            continue
        seen.add(key)
        kept.append(t)
        if len(kept) == k:
            break
    return kept


def lineage_lean(tables: dict, pairs: list, key_fn, lineage, collapsed: bool) -> dict:
    """For (unit, type) pairs: per arm, lineage counts of the top-5 terms of its `real` table, stromal share, top-10 terms by name."""
    out = {}
    for arm in ARMS:
        cnt, names, n = Counter(), Counter(), 0
        for unit, typ in pairs:
            for t in top_terms(tables[key_fn(arm, unit, typ)], collapsed):
                lin = lineage(t) or "none"
                cnt[lin] += 1
                names[(A.parse_term(t)[0], lin)] += 1
                n += 1
        st = sum(cnt[l] for l in STROMAL)
        out[arm] = {"n_pairs": len(pairs), "n_terms": n, "by_lineage": dict(cnt), "stromal_terms": st,
                    "stromal_share": (st / n) if n else None,
                    "top_terms": [[nm, lin, c] for (nm, lin), c in names.most_common(10)]}
    return out


# ============================================================ inputs
def breast_inputs(sample: str, mode: str) -> dict:
    if mode == "legacy":
        return dict(v41=RS.PLAIN_DIR / f"{sample}_cell_prop_global_results.csv", v2=RS.V2_PROPS / sample / "cell_prop.csv")
    return dict(v41=M3 / "p0" / sample / f"{sample}_cell_prop_global_results.csv", v2=M3 / f"{sample}_cell_prop_finetuned_results.csv")


def check_inputs_g4(sample: str, inp: dict, types: list) -> dict:
    """G4: column order, equal barcode sets, row sums 1, p0 equals the official marker share on shared barcodes."""
    p0, v2 = (pd.read_csv(inp[k], index_col=0) for k in ARMS)
    out = {}
    for k, d in (("v41", p0), ("v2", v2)):
        assert list(d.columns) == types, f"{k}: columns {list(d.columns)} != {types}"
        assert float(np.abs(d.sum(axis=1) - 1).max()) < 1e-6, f"{k}: row sums"
    assert set(p0.index) == set(v2.index) and not p0.index.has_duplicates, "barcode sets differ"
    ref = pd.read_csv(RS.PLAIN_DIR / f"{sample}_cell_prop_global_results.csv", index_col=0)
    sh = p0.index.intersection(ref.index)
    diff = float((p0.loc[sh] - ref.loc[sh][types]).abs().max().max())
    assert diff == 0.0, f"p0 differs from module3_marker_share_7type: {diff}"
    out.update(n_spots=int(len(p0)), n_shared_with_official=int(len(sh)), p0_vs_official_max_abs=diff)
    return out


# ============================================================ SACE + ORA drivers (unchanged scripts, subprocess)
def run_sace(idx: int, root: Path, inp: dict, variants: list, write_arms: str = ""):
    cmd = [sys.executable, "-u", str(HERE / "run_sace_layers.py"), "--sample-index", str(idx), "--out-root", str(root), "--profile", "v41",
           "--subset-arm", "v2", "--arm", f"v41={inp['v41']}", "--arm", f"v2={inp['v2']}", "--variants", ",".join(variants)]
    if write_arms:
        cmd += ["--write-arms", write_arms]
    subprocess.run(cmd, check=True)


def run_ora_unchanged(res_dir: Path, override_path: Path, suffix: str):
    subprocess.run([sys.executable, "-u", str(HERE / "run_ora.py"), "--results-dir", str(res_dir), "--protocol-override", str(override_path),
                    "--suffix", suffix], check=True)


def write_override(path: Path, sample: str):
    RR.write_new(path, json.dumps({**OVERRIDE, "specimens": [sample]}, indent=1))


def libs_all():
    return list(RS.libs().keys())


# ============================================================ second pass (full primary tables)
def tkey(test, arm, cond, unit, typ):
    return f"{test}|{arm}|{cond}|{unit}|{typ}"


def protocol_rows(layer_dir: Path, arm: str, sample: str, mapping: dict, markers: set, variants: list):
    props = pd.read_csv(layer_dir / "props.csv", index_col=0)
    abundance = props.mean(axis=0)
    rows, tables = [], {}
    for variant in variants:
        m = pd.read_parquet(layer_dir / f"layers_{variant}.parquet")
        layers = {t: m.loc[[t]] for t in m.index}
        genes = L.specificity_genes(layers, markers, n=SP["top_n"], min_mean=SP["min_mean"], pseudocount=SP["pseudocount"])
        bg = [g.upper() for g in L.expressed_background(layers, SP["min_mean"])]
        for t in m.index:
            row, tb = rec().row(genes[t], bg, t, mapping)
            k = tkey("P", arm, variant, sample, t)
            tables[k] = tb
            rows.append(dict(test="P", dataset="breast", arm=arm, condition=variant, unit=sample, type=t, abundance=float(abundance[t]),
                             evaluable=bool(abundance[t] >= MIN_ABUNDANCE), table_key=k, **row))
    return rows, tables


def testA_rows(layer_dir: Path, arm: str, sample: str, mapping: dict, markers: set, conds: dict):
    """ora_sense Test A on one arm (copy of run_ora_sense.breast_A with the layer directory and conditions as arguments)."""
    real = pd.read_parquet(layer_dir / "layers_real_marker.parquet")
    uni = pd.read_parquet(layer_dir / "layers_uniform_marker.parquet")
    props = pd.read_csv(layer_dir / "props.csv", index_col=0)
    expressed = real.sum(axis=0) >= S.MIN_MEAN
    cand = [g for g in real.columns[expressed.to_numpy()] if g not in markers]
    bg = [g.upper() for g in real.columns[expressed.to_numpy()]]
    rows, tables = [], {}
    for cond, variant in conds.items():
        lay = pd.read_parquet(layer_dir / f"layers_{variant}.parquet").reindex(columns=real.columns)
        genes = RS.gene_lists_from_scores(S.delta_table(lay, uni).T, cand)
        for t in real.index:
            row, tb = rec().row(genes[t], bg, t, mapping)
            k = tkey("A", arm, cond, sample, t)
            tables[k] = tb
            rows.append(dict(test="A", dataset="breast", arm=arm, condition=cond, unit=sample, type=t, evaluable=S.evaluable_share(props[t]),
                             abundance=float(props[t].mean()), table_key=k, **row))
    return rows, tables


A_CONDS_FULL = {"real": "real_marker", "null": "spot_permuted_marker", "w100": "uniform_marker", "w25": "blend25_marker", "w50": "blend50_marker",
                "w75": "blend75_marker"}
A_CONDS_GATE = {"real": "real_marker", "null": "spot_permuted_marker", "w100": "uniform_marker"}


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
            RR.fail(run, f"G3_{tag}_{arm}", g3[arm])
    return pd.DataFrame(prow), pd.DataFrame(arow), tabs, g3


# ============================================================ legacy gate
def gate(run: Path):
    run.mkdir(parents=True, exist_ok=False)
    g = run / "gate_legacy"
    g.mkdir()
    problems = F.check()
    if problems:
        RR.fail(run, "frozen_files", {"problems": problems[:20]})
    idx, sample = 0, SPECIMENS[0]
    inp = breast_inputs(sample, "legacy")
    run_sace(idx, g / "sace", inp, P_VARIANTS, write_arms="v41")
    ov = g / "override_gate.json"
    write_override(ov, sample)
    v41dir = g / "sace" / "v41"
    run_ora_unchanged(v41dir, ov, "_gate")
    libs = libs_all()
    stored = pd.read_csv(OLDV2 / "v41" / "ora_calls.csv", keep_default_na=False, na_values=[""])
    stored = stored[stored["specimen"] == sample]
    new = pd.read_csv(v41dir / "ora_calls_gate.csv", keep_default_na=False, na_values=[""])
    info = {"specimen": sample, "n_stored_rows": int(len(stored)), "n_new_rows": int(len(new))}
    bad_ora = compare_protocol(new, stored, libs)
    info["run_ora_vs_stored"] = {"n_mismatch": len(bad_ora), "mismatches": bad_ora[:20]}
    # second pass on the gate layers
    mapping = RS.load_mapping("term_mapping_7type.yaml")
    from ora_check.run_sace_layers import load_v41_profile

    markers = L.marker_gene_set(load_v41_profile())
    d = v41dir / "sace_layers" / sample
    prow, tabs = protocol_rows(d, "v41", sample, mapping, markers, P_VARIANTS)
    bad_sp = compare_protocol(gate_format(pd.DataFrame(prow), libs), new, libs)
    info["second_pass_vs_run_ora"] = {"n_mismatch": len(bad_sp), "mismatches": bad_sp[:20]}
    arow, tabs_a = testA_rows(d, "v41", sample, mapping, markers, A_CONDS_GATE)
    chk = RR.check_vs_stored(pd.DataFrame(arow), SENSE / "calls" / "A_breast_00.csv", "Test A v41 real/null/w100 vs stored A_breast_00")
    info["testA_vs_stored"] = chk
    # layers diagnostic
    diffs = {}
    for v in P_VARIANTS:
        a = pd.read_parquet(d / f"layers_{v}.parquet")
        b = pd.read_parquet(OLDV2 / "v41" / "sace_layers" / sample / f"layers_{v}.parquet")
        diffs[v] = float(np.abs(a.to_numpy() - b.reindex(index=a.index, columns=a.columns).to_numpy()).max())
    info["layers_max_abs_diff_vs_stored_diagnostic"] = diffs
    ok = (not bad_ora and not bad_sp and chk["n_mismatch"] == 0 and chk["compared"] == 21 and len(stored) == 49 and len(new) == 49)
    (g / "gate_report.json").write_text(json.dumps(info, indent=2, default=str))
    if not ok:
        RR.fail(run, "G1_legacy_gate", info)
    RR.write_new(g / "G1_PASS.json", json.dumps(info, indent=2, default=str))
    print("G1 PASS", flush=True)


# ============================================================ breast task
def breast(idx: int, run: Path):
    sample = SPECIMENS[idx]
    if not (run / "gate_legacy" / "G1_PASS.json").exists():
        RR.fail(run, f"breast_{idx:02d}_no_G1", {"error": "gate_legacy/G1_PASS.json missing"})
    problems = F.check()
    if problems:
        RR.fail(run, f"breast_{idx:02d}_frozen", {"problems": problems[:20]})
    inp = breast_inputs(sample, "rebuilt")
    for k, p in inp.items():
        assert p.is_file(), f"missing input {k}: {p}"
    types = list(pd.read_csv(inp["v41"], index_col=0).columns)
    g4 = check_inputs_g4(sample, inp, types)
    RR.write_new(run / "inputs" / f"breast_{idx:02d}.json", json.dumps({k: {"path": str(p), "sha256": RR.sha(p)} for k, p in inp.items()} | {"g4": g4}, indent=1))
    root = run / "rebuilt"
    run_sace(idx, root, inp, P_VARIANTS + BLENDS)
    meta = {}
    for arm in ARMS:
        m = json.loads((root / arm / "sace_layers" / sample / "run_meta.json").read_text())
        meta[arm] = {k: m[k] for k in ("n_spots", "n_arm_table_spots", "n_model_spots", "n_arm_spots_dropped")}
        if m["n_arm_spots_dropped"] != 0 or m["n_spots"] != m["n_model_spots"] or m["n_spots"] != g4["n_spots"]:
            RR.fail(run, f"breast_{idx:02d}_spots_{arm}", {"meta": m, "g4": g4})
    ov = run / "overrides" / f"spec_{idx:02d}.json"
    write_override(ov, sample)
    for arm in ARMS:
        run_ora_unchanged(root / arm, ov, f"_{idx:02d}")
    prow, arow, tabs, g3 = second_pass(root, sample, idx, A_CONDS_FULL, P_VARIANTS, run, f"_{idx:02d}")
    RR.write_new(run / "calls" / f"P_breast_{idx:02d}.csv", prow.to_csv(index=False))
    RR.write_new(run / "calls" / f"A_breast_{idx:02d}.csv", arow.to_csv(index=False))
    RR.write_tables(run / "tables" / f"breast_{idx:02d}.parquet", tabs)
    RR.write_new(run / "calls" / f"breast_{idx:02d}_meta.json", json.dumps({"sample": sample, "spots": meta, "g4": g4, "G3": g3}, indent=2))


# ============================================================ Xenium task (existing layers of ora_sense)
def xenium(r: int, run: Path):
    if not (run / "gate_legacy" / "G1_PASS.json").exists():
        RR.fail(run, f"xenium_{r}_no_G1", {"error": "gate_legacy/G1_PASS.json missing"})
    mapping = RS.load_mapping("term_mapping_xenium7.yaml")
    markers = RS._xen_markers()
    shas = {}
    for arm in ("oracle", "oracle_perm", "uniform", "v2", "v41"):
        base = SENSE / "sace" / arm / f"Xenium_region_{r}"
        for f in [base / "props.csv"] + [base / "layers" / f"{RS.XFILE[t]}_layer.csv" for t in RS.XTYPES]:
            shas[str(f.relative_to(HERE))] = RR.sha(f)
    RR.write_new(run / "inputs" / f"xenium_{r}.json", json.dumps(shas, indent=1))
    lay = {a: RS._xen_layer_means(SENSE, a, r) for a in ("oracle", "oracle_perm", "uniform", "v2", "v41")}
    panel = list(lay["oracle"][0].columns)
    tb_frozen = RS._testable(panel, mapping)
    tabs, xp, xa = {}, [], []

    def add(rows, test, arm, cond, genes, bg, props):
        for t in RS.XTYPES:
            row, tb = rec().row(genes[t], bg, t, mapping)
            k = tkey(test, arm, cond, f"region_{r}", t)
            tabs[k] = tb
            rows.append(dict(test=test, dataset="xenium", arm=arm, condition=cond, unit=f"region_{r}", type=t,
                             evaluable=bool(S.evaluable_share(props[t])), testable=bool(tb_frozen[t] >= S.MIN_PANEL_GENES), n_panel_genes=tb_frozen[t],
                             table_key=k, **row))

    for arm, src, cond in (("oracle", "oracle", "real"), ("oracle", "oracle_perm", "null"), ("uniform", "uniform", "real"), ("v2", "v2", "real"), ("v41", "v41", "real")):
        lm, props = lay[src]
        layers = {t: lm.loc[[t]] for t in lm.index}
        genes = L.specificity_genes(layers, markers, n=SP["top_n"], min_mean=SP["min_mean"], pseudocount=SP["pseudocount"])
        bg = L.expressed_background(layers, SP["min_mean"])
        add(xp, "XP", arm, cond, genes, bg, props)
    uni = lay["uniform"][0]
    oracle_lay = lay["oracle"][0]
    for arm, src, cond in (("oracle", "oracle", "real"), ("oracle", "oracle_perm", "null"), ("v2", "v2", "real"), ("v41", "v41", "real")):
        lm, props = lay[src]
        real_lay = lm if cond == "real" else oracle_lay
        expressed = real_lay.sum(axis=0) >= S.MIN_MEAN
        cand = [g for g in real_lay.columns[expressed.to_numpy()] if g not in markers]
        bg = list(real_lay.columns[expressed.to_numpy()])
        genes = RS.gene_lists_from_scores(S.delta_table(lm, uni).T, cand)
        add(xa, "XA", arm, cond, genes, bg, props)
    xa_df = pd.DataFrame(xa)
    chk = RR.check_vs_stored(xa_df, SENSE / "calls" / f"A_xenium_{r}.csv", f"XA region {r} vs stored A_xenium")
    all_testable = all(v >= S.MIN_PANEL_GENES for v in tb_frozen.values())
    RR.write_new(run / "calls" / f"XP_{r}.csv", pd.DataFrame(xp).to_csv(index=False))
    RR.write_new(run / "calls" / f"XA_{r}.csv", xa_df.to_csv(index=False))
    RR.write_tables(run / "tables" / f"xenium_{r}.parquet", tabs)
    RR.write_new(run / "calls" / f"xenium_{r}_meta.json", json.dumps({"region": r, "panel_genes": panel, "repro": chk, "all_testable_frozen": all_testable}, indent=2))
    if chk["n_mismatch"] or chk["compared"] != 28 or not all_testable:
        RR.fail(run, f"xenium_{r}_repro", {"repro": chk, "all_testable_frozen": all_testable})


# ============================================================ summary
def load_calls(run: Path, prefix: str, pattern: str) -> pd.DataFrame:
    fs = sorted((run / "calls").glob(f"{prefix}_{pattern}.csv"))
    df = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in fs], ignore_index=True)
    df["condition"] = df["condition"].astype(str)
    for c in [c for c in df.columns if c == "evaluable" or c == "passed" or c.startswith("pass_") or c == "testable"]:
        df[c] = to_bool(df[c])
    return df


def x_summary(x: pd.DataFrame) -> dict:
    r = x[x["evaluable"] & x["testable"]]
    out = {"passes": {f"{a}|{c}": int(g["passed"].sum()) for (a, c), g in x.groupby(["arm", "condition"])},
           "mean_nlp": {f"{a}|{c}": float(g["nlp"].mean()) for (a, c), g in r.groupby(["arm", "condition"])},
           "by_type": {f"{a}|{c}": {t: int(h["passed"].sum()) for t, h in g.groupby("type")} for (a, c), g in x.groupby(["arm", "condition"])},
           "n_evaluable": {f"{a}|{c}": int(len(g)) for (a, c), g in r.groupby(["arm", "condition"])}}
    return out


def summary(run: Path):
    if not (run / "gate_legacy" / "G1_PASS.json").exists():
        RR.fail(run, "summary_no_G1", {"error": "gate_legacy/G1_PASS.json missing"})
    problems = F.check()
    if problems:
        RR.fail(run, "summary_frozen", {"problems": problems[:20]})
    ov = M.load_overrides()
    mp_b, mp_x = RS.load_mapping("term_mapping_7type.yaml"), RS.load_mapping("term_mapping_xenium7.yaml")
    P = load_calls(run, "P_breast", "[0-9][0-9]")
    Ab = load_calls(run, "A_breast", "[0-9][0-9]")
    XP = load_calls(run, "XP", "[0-9]")
    XA = load_calls(run, "XA", "[0-9]")
    if P["unit"].nunique() != S.N_SPECIMENS or Ab["unit"].nunique() != S.N_SPECIMENS or XP["unit"].nunique() != S.N_REGIONS or XA["unit"].nunique() != S.N_REGIONS:
        RR.fail(run, "summary_units", {"P": int(P["unit"].nunique()), "A": int(Ab["unit"].nunique()), "XP": int(XP["unit"].nunique()), "XA": int(XA["unit"].nunique())})
    tables = {}
    for f in sorted((run / "tables").glob("*.parquet")):
        tables.update(RR.read_tables(f))
    panel = {f"region_{r}": json.loads((run / "calls" / f"xenium_{r}_meta.json").read_text())["panel_genes"] for r in range(S.N_REGIONS)}
    gmt = RS.libs()[PRIMARY]
    old = json.loads((OLDV2 / "readout.json").read_text())
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
                g3 = {"P": int((Pt["passed"].to_numpy() != P["passed"].to_numpy()).sum()), "A": int((At["passed"].to_numpy() != Ab["passed"].to_numpy()).sum()),
                      "XP": int(((XPt["passed"].to_numpy() != XP["passed"].to_numpy()) & tst).sum()), "XA": int(((XAt["passed"].to_numpy() != XA["passed"].to_numpy()) & tsa).sum()),
                      "nlp_max_abs": float(max(np.abs(Pt["nlp"].to_numpy(float) - P["nlp"].to_numpy(float)).max(), np.abs(At["nlp"].to_numpy(float) - Ab["nlp"].to_numpy(float)).max()))}
                if any(g3[k] for k in ("P", "A", "XP", "XA")) or g3["nlp_max_abs"] > 1e-9:
                    RR.fail(run, "G3_summary", g3)
            Pt["eff"] = Pt["passed"] & Pt["evaluable"]  # protocol call: a non-evaluable pair does not pass
            # ---- Part 1 (protocol)
            calls_v2 = {a: pd.DataFrame({"specimen": g["unit"], "variant": g["condition"], "type": g["type"], "passed": g["eff"]}) for a, g in Pt.groupby("arm")}
            ro = V.readout(calls_v2)
            pv = {a: {v: int(g[g["condition"] == v]["eff"].sum()) for v in P_VARIANTS} for a, g in Pt.groupby("arm")}
            pv_ungated = {a: {v: int(g[g["condition"] == v]["passed"].sum()) for v in P_VARIANTS} for a, g in Pt.groupby("arm")}
            pt = {a: {v: {t: int(h["eff"].sum()) for t, h in g[g["condition"] == v].groupby("type")} for v in P_VARIANTS} for a, g in Pt.groupby("arm")}
            ev_by_type = {a: {t: int(h["evaluable"].sum()) for t, h in g[g["condition"] == "real_marker"].groupby("type")} for a, g in Pt.groupby("arm")}
            labelP = {a: {"gap_vs_spot_permuted_marker": pv[a]["real_marker"] - pv[a]["spot_permuted_marker"],
                          "label": proportion_label(pv[a]["real_marker"], pv[a]["spot_permuted_marker"])} for a in ARMS}
            series = {a: Pt[(Pt["arm"] == a) & (Pt["condition"] == "real_marker")].set_index(["unit", "type"])["eff"] for a in ARMS}
            disc_p = discordance(series["v41"], series["v2"])
            pairs_b = [k for k in series["v41"].index if series["v41"][k] and not series["v2"][k]]
            pairs_c = [k for k in series["v41"].index if not series["v41"][k] and series["v2"][k]]
            kfn = lambda arm, u, t: tkey("P", arm, "real_marker", u, t)  # noqa: E731
            lean_p = {"v41_only": lineage_lean(tables, pairs_b, kfn, lb, collapsed), "v2_only": lineage_lean(tables, pairs_c, kfn, lb, collapsed),
                      "all_pairs": lineage_lean(tables, list(series["v41"].index), kfn, lb, collapsed)}
            # ---- Part 2 (Test A)
            ra = S.readout(At, XAt)
            pa = {a: {c: int(g[g["condition"] == c]["passed"].sum()) for c in A_CONDS_FULL} for a, g in At.groupby("arm")}
            labelA = {a: {"gap_vs_null": pa[a]["real"] - pa[a]["null"], "label": proportion_label(pa[a]["real"], pa[a]["null"])} for a in ARMS}
            sa = {a: At[(At["arm"] == a) & (At["condition"] == "real")].set_index(["unit", "type"])["passed"] for a in ARMS}
            disc_a = discordance(sa["v41"], sa["v2"])
            pb_a = [k for k in sa["v41"].index if sa["v41"][k] and not sa["v2"][k]]
            pc_a = [k for k in sa["v41"].index if not sa["v41"][k] and sa["v2"][k]]
            kfa = lambda arm, u, t: tkey("A", arm, "real", u, t)  # noqa: E731
            lean_a = {"v41_only": lineage_lean(tables, pb_a, kfa, lb, collapsed), "v2_only": lineage_lean(tables, pc_a, kfa, lb, collapsed),
                      "all_pairs": lineage_lean(tables, list(sa["v41"].index), kfa, lb, collapsed)}
            v3p, v3a = S.v3_status(XPt), S.v3_status(XAt)
            res[f"{mname}|{var}"] = {
                "P": {"passes_of_84": pv, "passes_ungated_of_84": pv_ungated, "by_type_by_variant": pt, "evaluable_real_marker_by_type": ev_by_type, "label": labelP,
                      "ora_v2_readout": ro, "discordance_real_marker": disc_p, "lineage_lean": lean_p},
                "A": {"passes_of_84": pa, "label": labelA, "sense_readout": ra, "discordance_real": disc_a, "lineage_lean": lean_a},
                "XP": {"v3": v3p, "summary": x_summary(XPt), "testable_counts": tcount},
                "XA": {"v3": v3a, "summary": x_summary(XAt)},
            }
            keep = ["dataset", "test", "arm", "condition", "unit", "type", "evaluable", "passed", "best_rank", "nlp"]
            allv.append(pd.concat([Pt.assign(passed_protocol=Pt["eff"])[keep + ["passed_protocol"]], At[keep], XPt[keep], XAt[keep]]).assign(mapping=mname, variant=var))
    pd.concat(allv, ignore_index=True).to_csv(run / "variants.csv", index=False)
    meta = {f: json.loads(f.read_text()) for f in sorted((run / "calls").glob("breast_[0-9][0-9]_meta.json"))}
    repro = {"G1": json.loads((run / "gate_legacy" / "G1_PASS.json").read_text()), "G3_summary": g3,
             "G3_runora_vs_second_pass": {m["sample"]: m["G3"] for m in meta.values()}, "G4": {m["sample"]: m["g4"] for m in meta.values()},
             "spots": {m["sample"]: m["spots"] for m in meta.values()},
             "xenium_XA_repro": {r: json.loads((run / "calls" / f"xenium_{r}_meta.json").read_text())["repro"] for r in range(S.N_REGIONS)}}
    old_keep = {k: old[k] for k in ("total", "per_type", "delta", "loss", "verdict", "null_gaps", "pass_rate_by_variant")}
    out = {"results": res, "old_spots_results_v2": old_keep, "frozen_files_check": F.check()}
    (run / "readout_sace.json").write_text(json.dumps(out, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    (run / "repro.json").write_text(json.dumps(repro, indent=1, default=str))
    (run / "readout_sace.md").write_text(render_md(out))
    print(json.dumps({k: {"P_labels": {a: v["P"]["label"][a]["label"] for a in ARMS}, "A_labels": {a: v["A"]["label"][a]["label"] for a in ARMS}} for k, v in res.items()}), flush=True)


def tab(rows, head):
    return ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]


def lean_lines(lean):
    out = []
    for d in ("v41_only", "v2_only", "all_pairs"):
        for arm in ARMS:
            x = lean[d][arm]
            sh = "n/a" if x["stromal_share"] is None else f"{x['stromal_share']:.3f}"
            out.append(f"- {d}, {arm}: pairs {x['n_pairs']}, top-5 terms {x['n_terms']}, by lineage {json.dumps(x['by_lineage'], sort_keys=True)}, stromal {x['stromal_terms']} (share {sh})")
    return out


def render_md(out):
    res = out["results"]
    o = out["old_spots_results_v2"]
    L_ = ["# ORA on the SACE layers, rebuilt breast spots, plus Xenium calibration (primary library CellMarker_2024)", "",
          "Exploratory, post-hoc. Rule: `ora_sace_rebuilt_DECISION_RULE.md`. Plausibility of gene-to-type links, not accuracy. Nothing here chooses a number.", "",
          f"Old spots (`results_v2_20261101_1501`, real_marker, passes of 84): {json.dumps(o['total'])}, verdict {o['verdict']}; "
          f"null totals v41 {json.dumps({k: v['null_total'] for k, v in o['null_gaps']['v41']['marker'].items()})}, v2 {json.dumps({k: v['null_total'] for k, v in o['null_gaps']['v2']['marker'].items()})}.", ""]
    for key, v in res.items():
        L_ += [f"## Readout `{key}` (mapping | call)", ""]
        p = v["P"]
        L_ += ["### Part 1: protocol (specificity genes), passes of 84 (evaluable = share >= 0.02)", ""]
        L_ += tab([[var, p["passes_of_84"]["v41"][var], p["passes_of_84"]["v2"][var]] for var in P_VARIANTS], ["variant", "v41", "v2"]) + [""]
        L_ += tab([[a, p["label"][a]["gap_vs_spot_permuted_marker"], p["label"][a]["label"]] for a in ARMS], ["arm", "real_marker - spot_permuted_marker", "label"]) + [""]
        r = p["ora_v2_readout"]
        L_ += [f"`ora_v2_lib.readout`: verdict {r['verdict']}, totals {json.dumps(r.get('total'))}, delta {r.get('delta')}, loss v41-v2 {json.dumps(r.get('loss'))}, not sensitive to proportions {json.dumps(r.get('not_sensitive_to_proportions'))}.", ""]
        types = sorted(p["by_type_by_variant"]["v41"]["real_marker"])
        L_ += ["Passes of 12 by type, v41 / v2:", ""]
        L_ += tab([[t] + [f"{p['by_type_by_variant']['v41'][var].get(t, 0)} / {p['by_type_by_variant']['v2'][var].get(t, 0)}" for var in P_VARIANTS] for t in types], ["type"] + P_VARIANTS) + [""]
        L_ += [f"Evaluable (of 12) real_marker by type: {json.dumps(p['evaluable_real_marker_by_type'])}. Ungated passes of 84 (real_marker): {json.dumps({a: p['passes_ungated_of_84'][a]['real_marker'] for a in ARMS})}.", ""]
        d = p["discordance_real_marker"]
        L_ += [f"Discordance real_marker: overall {json.dumps(d['overall'])}.", ""]
        L_ += tab([[t, x["passes_v41"], x["passes_v2"], x["v41_only"], x["v2_only"], f"{x['p_exact']:.4f}"] for t, x in d["by_type"].items()], ["type", "v41 passes", "v2 passes", "v41 only", "v2 only", "exact P"]) + [""]
        L_ += ["Lineage of the top-5 terms (discordant pairs):"] + lean_lines(p["lineage_lean"]) + [""]
        a = v["A"]
        L_ += ["### Part 2: Test A (delta of layers), passes of 84", ""]
        L_ += tab([[c, a["passes_of_84"]["v41"][c], a["passes_of_84"]["v2"][c]] for c in A_CONDS_FULL], ["condition", "v41", "v2"]) + [""]
        L_ += tab([[x, a["label"][x]["gap_vs_null"], a["label"][x]["label"]] for x in ARMS], ["arm", "real - null", "label"]) + [""]
        s = a["sense_readout"]
        L_ += [f"`ora_sense_lib.readout`: verdict {s['verdict']}, totals {json.dumps(s.get('total'))}, delta {s.get('delta')}, loss {json.dumps(s.get('loss'))}, V1 {json.dumps(s.get('V1'))}, V2 {json.dumps(s.get('V2'))}.", ""]
        da = a["discordance_real"]
        L_ += [f"Discordance real: overall {json.dumps(da['overall'])}.", ""]
        L_ += tab([[t, x["passes_v41"], x["passes_v2"], x["v41_only"], x["v2_only"], f"{x['p_exact']:.4f}"] for t, x in da["by_type"].items()], ["type", "v41 passes", "v2 passes", "v41 only", "v2 only", "exact P"]) + [""]
        L_ += ["Lineage of the top-5 terms (discordant pairs):"] + lean_lines(a["lineage_lean"]) + [""]
        L_ += ["### Part 3: Xenium calibration (35 calls)", ""]
        for nm in ("XP", "XA"):
            x = v[nm]
            L_ += [f"{nm} ({'protocol-style specificity' if nm == 'XP' else 'Test A delta'}): V3 = {json.dumps(x['v3'])}", "", f"passes of 35: {json.dumps(x['summary']['passes'])}; mean nlp: {json.dumps(x['summary']['mean_nlp'])}", ""]
            arms = list(x["summary"]["by_type"])
            xt = sorted(x["summary"]["by_type"][arms[0]])
            L_ += tab([[t] + [x["summary"]["by_type"][k].get(t, 0) for k in arms] for t in xt], ["type"] + arms) + [""]
    return "\n".join(L_) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["task", "summary"])
    ap.add_argument("a")
    ap.add_argument("b", nargs="?")
    a = ap.parse_args()
    t0 = time.time()
    if a.cmd == "task":
        n, run = int(a.a), Path(a.b)
        if n == GATE_TASK:
            gate(run)
        elif n < N_BREAST:
            breast(n, run)
        elif n < N_BREAST + N_XEN:
            xenium(n - N_BREAST, run)
        else:
            raise SystemExit("task id out of range")
    else:
        summary(Path(a.a))
    print(f"done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
