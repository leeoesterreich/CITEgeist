"""Analyze unassigned rate and doublet rate in gating-based cell classification."""

import json, logging, sys
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).parent.parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "CITEgeist"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
logger = logging.getLogger(__name__)

ASSIGN_PATH = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_gating_validation" / "region_2" / "assignments_no_neg.csv"
GT_PATH = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt" / "cell_type_assignments.csv"
THRESH_PATH = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_gating_validation" / "region_2" / "thresholds.json"
OUTPUT_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_gating_validation" / "region_2"

RCC_PROFILE_DICT = {
    "B cells": {"Major": ["CD20"], "Minor": ["CD45RA"], "Negative": ["CD3E", "CD68"]},
    "CD4+ T cells": {"Major": ["CD3E", "CD4"], "Minor": ["CD45RO"], "Negative": ["CD8A"]},
    "CD8+ T cells": {"Major": ["CD3E", "CD8A"], "Minor": ["GranzymeB"], "Negative": ["CD4"]},
    "Macrophages": {"Major": ["CD68", "CD163"], "Minor": ["CD16", "HLA-DR"], "Negative": ["CD3E"]},
    "Endothelial": {"Major": ["CD31"], "Minor": [], "Negative": ["PanCK", "alphaSMA"]},
    "Epithelial": {"Major": ["PanCK"], "Minor": [], "Negative": ["CD45", "CD31"]},
    "Fibroblasts": {"Major": ["alphaSMA"], "Minor": [], "Negative": ["PanCK", "CD31"]},
}
PROFILE_NAMES = list(RCC_PROFILE_DICT.keys())



def pct_str(num, denom):
    if denom == 0:
        return "0.0%"
    return f"{100 * num / denom:.1f}%"


def main():
    sep = "=" * 70
    logger.info(sep)
    logger.info("UNASSIGNED & DOUBLET ANALYSIS - Region 2 Gating Classification")
    logger.info(sep)
    logger.info("Loading data...")

    assign_df = pd.read_csv(ASSIGN_PATH, index_col=0)
    logger.info(f"Assignments: {assign_df.shape[0]} cells, columns={list(assign_df.columns)}")

    gt_df = pd.read_csv(GT_PATH, index_col=0)
    logger.info(f"Ground truth: {gt_df.shape[0]} cells")

    with open(THRESH_PATH) as f:
        thresholds = json.load(f)
    logger.info(f"Thresholds for {len(thresholds)} markers")

    from load_xenium_singlecell import load_xenium_singlecell
    logger.info("Loading Xenium single-cell protein data (region 2, 10000 cells)...")
    adata_gex, adata_protein = load_xenium_singlecell(region_id=2, max_cells=10000, seed=42)

    X_protein = adata_protein.X
    if hasattr(X_protein, "toarray"):
        X_protein = X_protein.toarray()
    X_protein = np.asarray(X_protein, dtype=np.float64)
    marker_names = list(adata_protein.var_names)
    cell_ids = list(adata_protein.obs_names)
    logger.info(f"Protein data: {X_protein.shape[0]} cells x {X_protein.shape[1]} markers")

    protein_df = pd.DataFrame(X_protein, index=cell_ids, columns=marker_names)
    common_cells = sorted(set(assign_df.index) & set(gt_df.index) & set(protein_df.index))
    logger.info(f"Common cells: {len(common_cells)}")

    assign_df = assign_df.loc[common_cells]
    gt_labels = gt_df.loc[common_cells, "cell_type"]
    protein_df = protein_df.loc[common_cells]
    pred_labels = assign_df["assigned_type"]
    doublet_flags = assign_df["doublet_flag"]
    n_total = len(common_cells)
    gt_in_profiles = set(PROFILE_NAMES)

    # ======== 1. UNASSIGNED ANALYSIS ========
    logger.info("\n" + sep)
    logger.info("1. UNASSIGNED ANALYSIS")

    unassigned_mask = pred_labels == "Unassigned"
    n_unassigned = unassigned_mask.sum()
    logger.info(f"Total cells: {n_total}, Unassigned: {n_unassigned} ({pct_str(n_unassigned, n_total)})")

    gt_of_unassigned = gt_labels[unassigned_mask]
    gt_counts = gt_of_unassigned.value_counts()
    logger.info("GT types of Unassigned cells:")
    n_capturable = 0
    n_uncapturable = 0
    for ct, count in gt_counts.items():
        tag = "IN profiles" if ct in gt_in_profiles else "NOT in profiles"
        logger.info(f"  {ct}: {count} ({pct_str(count, n_unassigned)}) [{tag}]")
        if ct in gt_in_profiles:
            n_capturable += count
        else:
            n_uncapturable += count
    logger.info(f"Capturable: {n_capturable} ({pct_str(n_capturable, n_unassigned)}), Uncapturable: {n_uncapturable} ({pct_str(n_uncapturable, n_unassigned)})")

    # 1b. Marker gate failures
    logger.info("Marker gate failures for capturable unassigned cells:")
    capturable_unassigned = unassigned_mask & gt_labels.isin(gt_in_profiles)
    capturable_idx = [c for c in common_cells if capturable_unassigned.loc[c]]
    marker_failure_by_type = {}
    for ct in PROFILE_NAMES:
        major_markers = RCC_PROFILE_DICT[ct]["Major"]
        ct_unassigned = [c for c in capturable_idx if gt_labels.loc[c] == ct]
        if not ct_unassigned:
            continue
        failures = {}
        for marker in major_markers:
            if marker not in marker_names:
                continue
            thresh = thresholds[marker]["threshold"]
            values = protein_df.loc[ct_unassigned, marker].values
            n_below = int((values < thresh).sum())
            n_ct = len(ct_unassigned)
            failures[marker] = {"n_below_threshold": n_below, "n_total": n_ct, "pct_failing": round(100 * n_below / max(n_ct, 1), 1), "mean_expr": round(float(values.mean()), 2), "median_expr": round(float(np.median(values)), 2), "threshold": round(thresh, 2)}
        marker_failure_by_type[ct] = failures
        logger.info(f"  {ct} ({len(ct_unassigned)} cells):")
        for marker, info in failures.items():
            logger.info(f"    {marker}: {info['n_below_threshold']}/{info['n_total']} ({info['pct_failing']}%%) below thresh={info['threshold']}, mean={info['mean_expr']}, median={info['median_expr']}")

    # 1c. Multi-marker failure patterns
    logger.info("Multi-marker failure patterns:")
    for ct in PROFILE_NAMES:
        major_markers = RCC_PROFILE_DICT[ct]["Major"]
        if len(major_markers) < 2:
            continue
        ct_unassigned = [c for c in capturable_idx if gt_labels.loc[c] == ct]
        if not ct_unassigned:
            continue
        n_ct = len(ct_unassigned)
        fail_patterns = Counter()
        for c in ct_unassigned:
            failed = tuple(sorted(m for m in major_markers if m in marker_names and protein_df.loc[c, m] < thresholds[m]["threshold"]))
            fail_patterns[failed] += 1
        logger.info(f"  {ct} ({n_ct} cells):")
        for pattern, count in fail_patterns.most_common():
            if pattern:
                logger.info(f"    Failed {list(pattern)}: {count} ({pct_str(count, n_ct)})")
            else:
                logger.info(f"    All PASSED but still unassigned: {count} ({pct_str(count, n_ct)})")

    # ======== 2. DOUBLET ANALYSIS ========
    logger.info("\n" + sep)
    logger.info("2. DOUBLET ANALYSIS")

    doublet_mask = doublet_flags.astype(str).str.lower() == "true"
    n_doublets = doublet_mask.sum()
    logger.info(f"Doublet cells: {n_doublets} ({pct_str(n_doublets, n_total)})")

    doublet_cells_df = assign_df.loc[doublet_mask, PROFILE_NAMES]
    doublet_pairs = Counter()
    doublet_type_counts = Counter()
    for idx, row in doublet_cells_df.iterrows():
        pos = [ct for ct in PROFILE_NAMES if row[ct] == 1.0]
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                doublet_pairs[tuple(sorted([pos[i], pos[j]]))] += 1
        for ct in pos:
            doublet_type_counts[ct] += 1

    logger.info("Most common doublet type pairs:")
    for pair, count in doublet_pairs.most_common(15):
        logger.info(f"  {pair[0]} + {pair[1]}: {count} ({pct_str(count, n_doublets)})")
    logger.info("Types involved in doublets:")
    for ct, count in doublet_type_counts.most_common():
        logger.info(f"  {ct}: {count} ({pct_str(count, n_doublets)})")

    evaluable = gt_labels.isin(gt_in_profiles) & (pred_labels != "Unassigned")
    doublet_eval = evaluable & doublet_mask
    nondoublet_eval = evaluable & ~doublet_mask
    d_acc = float((pred_labels[doublet_eval] == gt_labels[doublet_eval]).mean()) if doublet_eval.sum() > 0 else float("nan")
    nd_acc = float((pred_labels[nondoublet_eval] == gt_labels[nondoublet_eval]).mean()) if nondoublet_eval.sum() > 0 else float("nan")
    logger.info(f"Doublet accuracy: {d_acc:.4f} (n={doublet_eval.sum()}), Non-doublet accuracy: {nd_acc:.4f} (n={nondoublet_eval.sum()})")

    # Misclassification patterns
    d_misclass = pred_labels[doublet_eval] != gt_labels[doublet_eval]
    if d_misclass.sum() > 0:
        logger.info("Doublet misclassification patterns:")
        cp = Counter()
        for g, p in zip(gt_labels[doublet_eval][d_misclass], pred_labels[doublet_eval][d_misclass]):
            cp[(g, p)] += 1
        for (g, p), count in cp.most_common(10):
            logger.info(f"  GT={g} -> Pred={p}: {count}")

    # ======== 3. GT UNASSIGNED RATE COMPARISON ========
    logger.info("\n" + sep)
    logger.info("3. GT UNASSIGNED RATE COMPARISON")

    all_gt_types = gt_labels.value_counts()
    logger.info(f"GT label distribution ({n_total} cells):")
    for ct, count in all_gt_types.items():
        tag = "IN profiles" if ct in gt_in_profiles else "NOT in profiles"
        logger.info(f"  {ct}: {count} ({pct_str(count, n_total)}) [{tag}]")

    gt_unknown = ~gt_labels.isin(gt_in_profiles)
    n_gt_unknown = gt_unknown.sum()
    logger.info(f"GT cells NOT in 7 profiles: {n_gt_unknown} ({pct_str(n_gt_unknown, n_total)})")
    logger.info(f"Pred unassigned: {n_unassigned} ({pct_str(n_unassigned, n_total)})")

    both = gt_unknown & unassigned_mask
    n_both = both.sum()
    logger.info(f"Both GT-unknown AND pred-unassigned: {n_both}")
    if n_gt_unknown > 0:
        logger.info(f"  fraction of GT-unknown: {pct_str(n_both, n_gt_unknown)}")
    if n_unassigned > 0:
        logger.info(f"  fraction of pred-unassigned: {pct_str(n_both, n_unassigned)}")

    n_true_miss = int((~gt_unknown & unassigned_mask).sum())
    n_false_capture = int((gt_unknown & ~unassigned_mask).sum())
    logger.info(f"True misses (GT in profiles, pred unassigned): {n_true_miss} ({pct_str(n_true_miss, n_total)})")
    logger.info(f"False captures (GT NOT in profiles, pred assigned): {n_false_capture} ({pct_str(n_false_capture, n_total)})")

    # ======== 4. THRESHOLD ANALYSIS ========
    logger.info("\n" + sep)
    logger.info("4. THRESHOLD ANALYSIS")
    logger.info("Per-marker threshold summary:")

    threshold_summary = []
    for marker in marker_names:
        if marker not in thresholds:
            continue
        t = thresholds[marker]
        tv = t["threshold"]
        method = t["method"]
        snr = t.get("snr", 0)
        quality = t.get("quality", "unknown")
        gm = t.get("gmm_means", [])
        vals = protein_df[marker].values
        pct_pos = 100 * (vals >= tv).sum() / len(vals)
        gl = gm[0] if gm else float("nan")
        gh = gm[-1] if gm else float("nan")
        logger.info(f"  {marker}: thresh={tv:.2f}, method={method}, SNR={snr:.2f}, quality={quality}, pos={pct_pos:.1f}%%, gmm=[{gl:.1f}, {gh:.1f}]")
        threshold_summary.append({"marker": marker, "threshold": tv, "method": method, "snr": snr, "quality": quality, "pct_positive": round(pct_pos, 1), "gmm_mean_low": round(gl, 1), "gmm_mean_high": round(gh, 1)})

    logger.info("Major marker expression in GT-positive cells:")
    for ct, profile in RCC_PROFILE_DICT.items():
        for marker in profile["Major"]:
            if marker not in marker_names:
                continue
            thresh = thresholds[marker]["threshold"]
            gpm = gt_labels == ct
            ngp = gpm.sum()
            if ngp == 0:
                continue
            vals = protein_df.loc[gpm, marker].values
            pa = 100 * (vals >= thresh).sum() / ngp
            logger.info(f"  {marker} ({ct}): thresh={thresh:.1f}, n={ngp}, above={pa:.1f}%%, mean={vals.mean():.1f}, P25={np.percentile(vals, 25):.1f}, P50={np.percentile(vals, 50):.1f}, P75={np.percentile(vals, 75):.1f}")

    logger.info("Markers where >20%% of GT+ cells fall below threshold:")
    found = False
    for ct, profile in RCC_PROFILE_DICT.items():
        for marker in profile["Major"]:
            if marker not in marker_names:
                continue
            thresh = thresholds[marker]["threshold"]
            gpm = gt_labels == ct
            ngp = gpm.sum()
            if ngp == 0:
                continue
            vals = protein_df.loc[gpm, marker].values
            pb = 100 * (vals < thresh).sum() / ngp
            if pb > 20:
                found = True
                logger.info(f"  WARNING: {marker} ({ct}): {pb:.1f}%% below thresh {thresh:.1f} (mean={vals.mean():.1f}, median={np.median(vals):.1f})")
    if not found:
        logger.info("  None found")

    # Save results
    results = {
        "n_total": n_total, "n_unassigned": int(n_unassigned),
        "pct_unassigned": round(100 * n_unassigned / n_total, 2),
        "n_doublets": int(n_doublets), "pct_doublets": round(100 * n_doublets / n_total, 2),
        "unassigned_gt_breakdown": {str(k): int(v) for k, v in gt_counts.items()},
        "unassigned_capturable": int(n_capturable), "unassigned_uncapturable": int(n_uncapturable),
        "marker_failure_by_type": marker_failure_by_type,
        "doublet_pairs_top10": [{"pair": list(p), "count": int(c)} for p, c in doublet_pairs.most_common(10)],
        "doublet_accuracy": round(d_acc, 4) if not np.isnan(d_acc) else None,
        "nondoublet_accuracy": round(nd_acc, 4) if not np.isnan(nd_acc) else None,
        "gt_unknown_count": int(n_gt_unknown), "gt_unknown_pct": round(100 * n_gt_unknown / n_total, 2),
        "overlap_both_unknown": int(n_both), "true_misses": n_true_miss,
        "false_captures": n_false_capture, "threshold_summary": threshold_summary,
    }
    out_path = OUTPUT_DIR / "unassigned_doublet_analysis.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    logger.info(f"Results saved to {out_path}")
    logger.info(sep)
    logger.info("ANALYSIS COMPLETE")


if __name__ == "__main__":
    main()
