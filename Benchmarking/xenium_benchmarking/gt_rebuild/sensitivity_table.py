"""Annotation-sensitivity supplementary table (6-type only).

Primary column: strict 6-type result on the rnaprot GT. Sensitivity columns: the RNA-only (rnagt)
strict 6-type result, and the per-cell variant GTs scored by reference_audit/singler_variants
(2026-10-21 skip-NaN evaluator; per-type r there is NOT comparable to the strict per-type headline).

CITEgeist_v41p (the paper's "CITEgeist") has no row in singler_variants/scores_all.csv, so
``rescore_v41p`` scores only it on the same variant GT dirs with the same evaluator call and
writes ``reference_audit/singler_variants_v41p/scores_v41p.csv``. scores_all.csv is never touched.

v41 (``--with-p138``): adds the row ``CITEgeist (P138)`` = plain marker share with the ms_v41 dictionary, the
paper's v41 CITEgeist. Its primary column is the strict rnaprot P138 result (``results_rnaprot_p138``), its RNA-only
column the strict rnagt P138 result (``results_rnagt_p138``, spots_baseline/score_p138_strict.py --labelset rnagt),
and its variant columns come from ``rescore_p138`` (``reference_audit/singler_variants_p138/scores_p138.csv``, same
evaluator call as the other rows). The v41 table is written to ``results_rnaprot_v41/`` so the v40/v41p table in
``results_rnaprot/`` is never touched.

Entry point (from Benchmarking/xenium_benchmarking)::

    XEN_LABELSET=rnaprot python -m gt_rebuild.sensitivity_table [--rescore]
    XEN_LABELSET=rnaprot python -m gt_rebuild.sensitivity_table --with-p138 [--rescore-p138]
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

XB = Path(__file__).resolve().parents[1]
REPO = XB.parents[1]
for _p in (str(XB), str(XB / "evaluation" / "src"), str(REPO / "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import benchmark_paths as bp  # noqa: E402
from gt_rebuild import manifest  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/gt_rebuild/sensitivity_table.py"
SV = XB / "reference_audit" / "singler_variants"
SV_V41P = XB / "reference_audit" / "singler_variants_v41p"
SCORES_ALL = SV / "scores_all.csv"
SCORES_V41P = SV_V41P / "scores_v41p.csv"
SV_P138 = XB / "reference_audit" / "singler_variants_p138"
SCORES_P138 = SV_P138 / "scores_p138.csv"
P138_ID = "CITEgeist_p138"
P138_DISPLAY = "CITEgeist (P138)"  # display name of the P138 arm
V2_ID = "CITEgeist_v2"  # v42 frozen-v2 arm; the display name below is the one the v42 builder expects (SENS_TABLE_V2_DISPLAY)
V2_DISPLAY = "CITEgeist (v2)"
P138_RESULT = "SPOTS_P138_6type_results.json"  # score_p138_strict.py output name (both labelsets)
CT6 = ["B cells", "T cells", "Macrophages", "Fibroblasts", "Epithelial", "Endothelial"]
VARIANT_GTS = ["canonical", "control", "A_NKmerged", "A_NKremoved", "B"]  # recomputed_canonical dropped
VARIANT_DISPLAY = {
    "canonical": "March per-cell labels",
    "control": "May per-cell labels",
    "A_NKmerged": "A_NKmerged",
    "A_NKremoved": "A_NKremoved",
    "B": "B",
}
# method id -> (display label, displayed, strict result json name)
METHODS = {
    "CITEgeist_v41p": ("CITEgeist", True, "CITEgeist_v41p_6type_results.json"),
    "MarkerShare_plain": ("CITEgeist (plain marker share)", True, "CITEgeist_6type_results.json"),
    "RCTD": ("RCTD", True, "RCTD_6type_results.json"),
    "Cell2Location": ("Cell2Location", True, "Cell2Location_6type_results.json"),
    "CARD": ("CARD", True, "CARD_results.json"),
    "Seurat": ("Seurat", True, "Seurat_6type_results.json"),
    "Tangram": ("Tangram", True, "Tangram_6type_results.json"),
    "CITEgeist": ("CITEgeist (QP)", False, None),
    "MarkerShare_corrected": ("CITEgeist (corrected marker share)", False, None),
}
NOTE = (
    "6-type only (singler_variants scored 6-type predictions only). Primary = rnaprot strict 6-type "
    "(per-type r = mean over types of the region-mean per-type r; flat r = overall_mean_pearson_r). "
    "Annotation-sensitivity columns (RNA-only GT, per-cell variant GTs) are NOT primary. Variant per-type r "
    "comes from the 2026-10-21 skip-NaN evaluator and is NOT comparable to the strict per-type headline "
    "(Seurat B cells are NaN there)."
)


def _variant_row(gt_dir, pred_dir, gtn, method):
    """Same evaluator call and row layout as reference_audit/singler_variants/score_variants.py."""
    from evaluate_benchmark import evaluate_all_regions

    s = evaluate_all_regions(
        gt_dir=str(gt_dir),
        pred_dir=str(pred_dir),
        n_regions=5,
        output_path=None,
        prefix="Xenium",
        use_achievable_7=False,
        combine_t_cells=True,
    )["summary"]
    return dict(
        gt=gtn,
        method=method,
        pooled_r=s["overall_mean_pearson_r"],
        std_r=s["overall_std_pearson_r"],
        rmse=s["overall_mean_rmse"],
        mae=s["overall_mean_mae"],
        jsd=s.get("overall_mean_jsd", np.nan),
        T_r=s["T cells_mean_r"],
        total_spots=s["total_spots"],
        **{f"{c}_r": s[f"{c}_mean_r"] for c in CT6 if c != "T cells"},
    )


def variant_gt_dirs(sv=SV):
    d = {g: Path(sv) / "gt" / g for g in VARIANT_GTS if g != "canonical"}
    d["canonical"] = XB / "ground_truth_singler" / "ground_truth_6type"
    return d


def rescore_v41p(pred_dir=None, out_csv=SCORES_V41P, gt_dirs=None):
    """Score CITEgeist_v41p 6-type predictions (read-only) on the variant GTs; write a NEW csv."""
    pred_dir = bp.citegeist_v41p_pred_dir(6, "rnaprot") if pred_dir is None else pred_dir
    gt_dirs = variant_gt_dirs() if gt_dirs is None else gt_dirs
    rows = [_variant_row(gt_dirs[g], pred_dir, g, "CITEgeist_v41p") for g in VARIANT_GTS]
    df = pd.DataFrame(rows)
    out_csv = Path(out_csv)
    if out_csv.exists():
        raise FileExistsError(f"refusing to overwrite {out_csv}")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def rescore_p138(pred_dir=None, out_csv=SCORES_P138, gt_dirs=None, method_id=P138_ID):
    """Score CITEgeist P138 6-type predictions (read-only) on the variant GTs; write a NEW csv.

    ``method_id`` names the scored arm in the csv (default P138; the v2 arm passes CITEgeist_v2 with its own pred_dir and out_csv)."""
    pred_dir = bp.citegeist_p138_pred_dir(6) if pred_dir is None else pred_dir
    gt_dirs = variant_gt_dirs() if gt_dirs is None else gt_dirs
    rows = [_variant_row(gt_dirs[g], pred_dir, g, method_id) for g in VARIANT_GTS]
    df = pd.DataFrame(rows)
    out_csv = Path(out_csv)
    if out_csv.exists():
        raise FileExistsError(f"refusing to overwrite {out_csv}")
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    return df


def _strict(json_dir, name):
    if name is None or not (Path(json_dir) / name).is_file():
        return np.nan, np.nan
    s = json.loads((Path(json_dir) / name).read_text())["summary"]
    return s["per_type_mean_r"], s["overall_mean_pearson_r"]


def build_table(
    scores_csv,
    v41p_scores_csv,
    primary_json_dir_rnaprot,
    primary_json_dir_rnagt,
    rnagt_v41p_json_dir=None,
    p138=None,
):
    """One row per method. Raises if a variant GT is missing from the scores.

    p138: optional dict(scores_csv=, rnaprot_dir=, rnagt_dir=) adding the ``CITEgeist (P138)`` row; both strict result
    files must exist (a missing one raises rather than leaving a NaN in the v41 CITEgeist row). Optional keys method_id /
    display / result / hide (method ids to keep in the csv but not display) (defaults: the P138 constants) let the same slot carry the frozen-v2 arm (v42 stage B15)."""
    parts = [pd.read_csv(scores_csv), pd.read_csv(v41p_scores_csv)]
    methods = dict(METHODS)
    pid = (p138 or {}).get("method_id", P138_ID)
    if p138 is not None:
        pdisp, pres = p138.get("display", P138_DISPLAY), p138.get("result", P138_RESULT)
        for k in ("rnaprot_dir", "rnagt_dir"):
            f = Path(p138[k]) / pres
            if not f.is_file():
                raise FileNotFoundError(f"{pid} strict result missing: {f}")
        parts.append(pd.read_csv(p138["scores_csv"]))
        for hidden in p138.get("hide", ()):  # v2 table: v41-lineage CITEgeist rows are kept in the csv but not displayed
            label, _disp, jname = methods[hidden]
            methods[hidden] = (label, False, jname)
        methods = {pid: (pdisp, True, pres), **methods}
    scores = pd.concat(parts, ignore_index=True)
    missing = [g for g in VARIANT_GTS if g not in set(scores["gt"])]
    if missing:
        raise KeyError(f"variant GT(s) missing from scores: {missing}")
    rnagt_v41p_json_dir = primary_json_dir_rnagt if rnagt_v41p_json_dir is None else rnagt_v41p_json_dir
    rows = []
    for mid, (label, displayed, jname) in methods.items():
        r = {"method_id": mid, "display_name": label, "displayed": displayed}
        r["primary_gt"] = "rnaprot (strict 6-type)"
        prim_dir = p138["rnaprot_dir"] if mid == pid and p138 is not None else primary_json_dir_rnaprot
        r["primary_per_type_r"], r["primary_flat_r"] = _strict(prim_dir, jname)
        sens_dir = {"CITEgeist_v41p": rnagt_v41p_json_dir, pid: (p138 or {}).get("rnagt_dir")}.get(
            mid, primary_json_dir_rnagt
        )
        r["sensitivity_rnagt_per_type_r"], r["sensitivity_rnagt_flat_r"] = _strict(sens_dir, jname)
        for g in VARIANT_GTS:
            sub = scores[(scores["gt"] == g) & (scores["method"] == mid)]
            if len(sub) != 1:
                raise KeyError(f"expected one score row for {mid} on {g}, found {len(sub)}")
            x = sub.iloc[0]
            r[f"sens_{g}_pooled_r"] = x["pooled_r"]
            for c in CT6:
                r[f"sens_{g}_{c}_r"] = x["T_r"] if c == "T cells" else x[f"{c}_r"]
        r["variant_role"] = "annotation sensitivity (not primary)"
        r["note"] = NOTE
        rows.append(r)
    df = pd.DataFrame(rows)
    df.attrs["variant_display"] = VARIANT_DISPLAY
    return df


def displayed_table(df):
    """Paper-facing view: displayed rows only, raw variant ids replaced by paper wording."""
    d = df[df["displayed"]].drop(columns=["displayed", "method_id"])
    ren = {}
    for g, disp in VARIANT_DISPLAY.items():
        for c in d.columns:
            if c.startswith(f"sens_{g}_"):
                ren[c] = c.replace(f"sens_{g}_", f"{disp} | ", 1)
    return d.rename(columns=ren)


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rescore", action="store_true", help="first create scores_v41p.csv (fails if present)")
    ap.add_argument("--with-p138", action="store_true", help="v41: add the CITEgeist (P138) row; write to results_rnaprot_v41/")
    ap.add_argument("--rescore-p138", action="store_true", help="first create scores_p138.csv (fails if present)")
    ap.add_argument("--with-v2", action="store_true", help="v42 stage B15: the frozen-v2 row in the P138 slot; write to results_rnaprot_v2/")
    ap.add_argument("--rescore-v2", action="store_true", help="first create singler_variants_v2/scores_v2.csv (fails if present)")
    ap.add_argument("--v2-display", default=V2_DISPLAY, help="display name of the v2 row")
    args = ap.parse_args(argv)
    if args.rescore_p138 and not args.with_p138:
        ap.error("--rescore-p138 requires --with-p138")
    if args.rescore_v2 and not args.with_v2:
        ap.error("--rescore-v2 requires --with-v2")
    if args.with_v2 and args.with_p138:
        ap.error("--with-v2 and --with-p138 are exclusive (one arm in the slot)")
    manifest.require_rnagt_signoff(labelset="rnaprot")
    manifest.require_rnaprot_spec()
    import provenance

    if args.rescore:
        rescore_v41p()
    if args.rescore_p138:
        rescore_p138()
    p138 = None
    if args.with_p138:
        p138 = {
            "scores_csv": SCORES_P138,
            "rnaprot_dir": XB / "evaluation" / "results_rnaprot_p138" / "method_comparison_singler_6type",
            "rnagt_dir": XB / "evaluation" / "results_rnagt_p138" / "method_comparison_singler_6type",
        }
    if args.with_v2:
        sys.path.insert(0, str(REPO))
        from Benchmarking.xenium_benchmarking.rna_refinement import v2_official as vo

        v2_scores = vo.variants_scores_dir() / "scores_v2.csv"
        if args.rescore_v2:
            rescore_p138(pred_dir=vo.pred_root(6), out_csv=v2_scores, method_id=V2_ID)
        p138 = {
            "scores_csv": v2_scores,
            "rnaprot_dir": vo.results_dir("rnaprot") / "method_comparison_singler_6type",
            "rnagt_dir": vo.results_dir("rnagt") / "method_comparison_singler_6type",
            "method_id": V2_ID,
            "display": args.v2_display,
            "result": "CITEgeist_6type_results.json",
            "hide": ("CITEgeist_v41p", "MarkerShare_plain"),
        }
    out_dir = vo.results_dir("rnaprot") if args.with_v2 else bp.eval_results_dir("rnaprot_v41" if args.with_p138 else "rnaprot")
    csv_p, json_p = out_dir / "sensitivity_supplement_6type.csv", out_dir / "sensitivity_supplement_6type.json"
    for p in (csv_p, json_p):
        if p.exists():
            raise FileExistsError(f"refusing to overwrite {p}")
    jg = bp.eval_results_dir("rnagt") / "method_comparison_singler_6type"
    jv = bp.eval_results_dir("rnagt_v41p") / "method_comparison_singler_6type"
    jr = bp.eval_results_dir("rnaprot") / "method_comparison_singler_6type"  # primary competitor results (v40/v41p)
    df = build_table(SCORES_ALL, SCORES_V41P, jr, jg, jv, p138=p138)
    score_files = [SCORES_ALL, SCORES_V41P] + ([Path(p138["scores_csv"])] if p138 else [])
    inputs = {str(p.relative_to(REPO)): _sha(p) for p in score_files}
    for d in (jr, jg, jv) + ((Path(p138["rnaprot_dir"]), Path(p138["rnagt_dir"])) if p138 else ()):
        for p in sorted(d.glob("*_results.json")):
            inputs[str(p.relative_to(REPO))] = _sha(p)
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_p, index=False)
    provenance.write_stamped_json(
        json_p,
        {
            "note": NOTE,
            "n_types": 6,
            "inputs_sha256": inputs,
            "displayed": json.loads(displayed_table(df).to_json(orient="records")),
        },
        PRODUCER,
    )
    return df


if __name__ == "__main__":
    main()
