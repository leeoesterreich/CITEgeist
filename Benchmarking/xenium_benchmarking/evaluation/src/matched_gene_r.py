"""Matched-gene GEX r (manuscript Fig 3F metric): CITEgeist SACE pooled r restricted to Cell2Location's gene set.

Wraps spots_baseline/eval_sace_arms.matched_pooled_r. Layer dirs come from benchmark_paths (XEN_LABELSET), so
v40 reproduces 0.727 and march gives the v41 value. Writes a stamped gex_eval_dir()/matched_gene_r.json:
  {"citegeist_matched_r", "c2l_flat_r", "n_c2l_genes"}   (c2l_flat_r = Cell2Location ALL row flattened_r_mean
  of xenium_gex_comparison_7type_summary.csv in the same dir)
Under rnagt / rnaprot (strict labelsets) the rnagt sign-off (+ the approved rnaprot spec) is required first, and the
output is a stamped JSON (code commit, input sha256s) with one entry per SACE variant (CITEgeist_SACE, and
CITEgeist_v41p_SACE under rnaprot) and Cell2Location, each listing the (region, type) blocks scored and excluded
(gt_all_zero / insufficient_overlap / ...): the variants can score different block sets.
Usage: matched_gene_r.py [--sace-dir DIR] [--c2l-dir DIR] [--out JSON]
"""
import argparse
import hashlib
import os
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[4]
XEN = REPO / "Benchmarking" / "xenium_benchmarking"
sys.path.insert(0, str(XEN))
sys.path.insert(0, str(XEN / "spots_baseline"))
sys.path.insert(0, str(XEN / "evaluation" / "src"))
sys.path.insert(0, str(REPO / "scripts"))
import benchmark_paths as bp  # noqa: E402
import eval_sace_arms as E  # noqa: E402
from gt_rebuild import manifest  # noqa: E402

PRODUCER = "Benchmarking/xenium_benchmarking/evaluation/src/matched_gene_r.py"


def c2l_union_genes(c2l_dir):
    """Union of upper-cased gene names over every Cell2Location layer CSV (all regions/types)."""
    genes = set()
    for reg in range(5):
        for f in (Path(c2l_dir) / f"Xenium_region_{reg}" / "layers").glob("*.csv"):
            genes |= set(pd.read_csv(f, index_col=0, nrows=1).columns.str.upper())
    return genes


def compute(sace_dir, c2l_dir, summary_csv):
    genes = c2l_union_genes(c2l_dir)
    assert genes, f"no Cell2Location layer CSVs under {c2l_dir}"
    summ = pd.read_csv(summary_csv)
    row = summ[(summ["method"] == "Cell2Location") & (summ["cell_type"] == "ALL")]
    assert len(row) == 1, f"expected one Cell2Location ALL row in {summary_csv}"
    return {
        "citegeist_matched_r": E.matched_pooled_r(sace_dir, genes),
        "c2l_flat_r": float(row["flattened_r_mean"].iloc[0]),
        "n_c2l_genes": len(genes),
    }


def sace_roster(labelset):
    """SACE layer dirs to score: the original CITEgeist SACE, plus v41p under the v41p labelsets and P138 under
    the P138 labelsets."""
    r = {"CITEgeist_SACE": bp.citegeist_sace_dir(labelset)}
    if labelset in bp.V41P_LABELSETS:
        r["CITEgeist_v41p_SACE"] = bp.citegeist_v41p_sace_dir(labelset)
    if labelset in bp.P138_LABELSETS:
        r["CITEgeist_p138_SACE"] = bp.citegeist_p138_sace_dir(labelset)
    return r


def _rel(p):
    if p is None:
        return None
    try:
        return str(Path(p).relative_to(REPO))
    except ValueError:
        return str(p)


def compute_strict(sace_dirs, c2l_dir, summary_csv, primary="CITEgeist_SACE", all_spots=False):
    """Matched-gene pooled r for each SACE variant and for Cell2Location itself, on C2L's gene set.

    ``primary`` names the roster entry reported as ``citegeist_matched_r`` (default unchanged; the v2 arm passes CITEgeist_v2_SACE)."""
    genes = c2l_union_genes(c2l_dir)
    assert genes, f"no Cell2Location layer CSVs under {c2l_dir}"
    summ = pd.read_csv(summary_csv)
    row = summ[(summ["method"] == "Cell2Location") & (summ["cell_type"] == "ALL")]
    assert len(row) == 1, f"expected one Cell2Location ALL row in {summary_csv}"
    methods, inputs = {}, {}
    for name, d in {**sace_dirs, "Cell2Location": c2l_dir}.items():
        r, blocks = E.matched_pooled_r(d, genes, return_blocks=True, all_spots=all_spots)
        scored = [b for b in blocks if b["scored"]]
        excl = [b for b in blocks if not b["scored"]]
        for b in blocks:
            for k in ("gt_path", "pred_path"):
                b[k] = _rel(b[k])
        for b in scored:
            for k in ("gt_path", "pred_path"):
                f = REPO / b[k]
                inputs[b[k]] = hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else None
        methods[name] = {
            "matched_pooled_r": r,
            "n_blocks_scored": len(scored),
            "n_blocks_excluded": len(excl),
            "n_blocks_gt_all_zero": sum(b["reason"] == "gt_all_zero" for b in excl),
            "n_spots_total": int(sum(b["n_spots"] for b in scored)),
            "blocks_scored": [{k: b[k] for k in ("region", "cell_type", "n_gt_spots", "n_spots", "n_genes")} for b in scored],
            "blocks_excluded": [{k: b[k] for k in ("region", "cell_type", "reason", "n_gt_spots", "n_spots") if k in b} for b in excl],
        }
    return {
        "citegeist_matched_r": methods[primary]["matched_pooled_r"],
        **(
            {"citegeist_v41p_matched_r": methods["CITEgeist_v41p_SACE"]["matched_pooled_r"]}
            if "CITEgeist_v41p_SACE" in methods
            else {}
        ),
        **(
            {"citegeist_p138_matched_r": methods["CITEgeist_p138_SACE"]["matched_pooled_r"]}
            if "CITEgeist_p138_SACE" in methods
            else {}
        ),
        "c2l_flat_r": float(row["flattened_r_mean"].iloc[0]),
        "c2l_matched_pooled_r": methods["Cell2Location"]["matched_pooled_r"],
        "like_for_like_partner": {"citegeist_matched_r": "c2l_matched_pooled_r"},
        "n_c2l_genes": len(genes),
        "note": (
            "matched_pooled_r = one Pearson r over all scored blocks (log1p, pred>0 or GT>0) on Cell2Location's gene set; "
            "c2l_flat_r = Cell2Location ALL row of the GEX summary (per-block canonical metric). "
            "Only c2l_matched_pooled_r is the like-for-like partner of citegeist_matched_r (same pooled metric, same genes, same "
            "blocks); c2l_flat_r is a different, per-block metric. C2L is dense over all spots while SACE layers cover only "
            "detected spots, so the spots entering the pooled r differ (see methods.*.n_spots_total) even though genes and "
            "blocks match. Variants can score different block sets: compare blocks_scored."
        ),
        "methods": methods,
        "input_sha256": inputs,
    }


def main(argv=None):
    manifest.require_rnagt_signoff()  # first: before reading or writing anything; no-op outside rnagt / rnaprot
    ap = argparse.ArgumentParser()
    ap.add_argument("--sace-dir", default=str(bp.citegeist_sace_dir()))
    ap.add_argument("--c2l-dir", default=str(bp.pred_dir("Cell2Location", 7)))
    ap.add_argument("--out", default=str(bp.gex_eval_dir() / "matched_gene_r.json"))
    ap.add_argument("--v2-sace-dir", default=None, help="v42 stage X4: score the frozen-v2 SACE layers (as CITEgeist_v2_SACE) instead of the labelset roster")
    ap.add_argument("--all-spots", action="store_true", help="v42 fair standard: score all GT spots, an absent spot is zero (default off)")
    ap.add_argument("--summary", default=None, help="GEX summary csv with the Cell2Location ALL row (default: next to --out)")
    a = ap.parse_args(argv)
    from provenance import write_stamped_json

    summary = Path(a.summary) if a.summary else Path(a.out).with_name("xenium_gex_comparison_7type_summary.csv")
    ls = os.environ.get("XEN_LABELSET", "")
    if ls in bp.STRICT_LABELSETS:
        if a.v2_sace_dir:
            roster, primary = {"CITEgeist_v2_SACE": Path(a.v2_sace_dir)}, "CITEgeist_v2_SACE"
        else:
            roster, primary = sace_roster(ls), "CITEgeist_SACE"
            roster["CITEgeist_SACE"] = Path(a.sace_dir)
        out = compute_strict(roster, a.c2l_dir, summary, primary, all_spots=a.all_spots)
        out["labelset"] = ls
        out["gt_gex_dir"] = _rel(E.C.GT_ROOT)
    else:
        out = compute(a.sace_dir, a.c2l_dir, summary)
    print({k: v for k, v in out.items() if k not in ("methods", "input_sha256")}, flush=True)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    write_stamped_json(Path(a.out), out, PRODUCER)


if __name__ == "__main__":
    main()
