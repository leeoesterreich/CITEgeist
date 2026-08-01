"""Aggregate competitor GEX biology validation results and classify failure modes.

Collects Moran's I, COMMOT, and Module 4 NMF program results across methods
and samples, compares against CITEgeist baselines, and produces failure-mode
classification tables.

Outputs (saved to output/competitor_gex/comparison/):
    morans_comparison.csv   - Candidate genes side-by-side across methods
    commot_comparison.csv   - Raw COMMOT pathway scores across methods
    failure_mode_table.csv  - Per method x sample x analysis classification
    failure_mode_summary.csv - Pivot counts (method x analysis x category)
"""
from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from constants import (
    CANDIDATE_GENES,
    COMMOT_PATHWAYS,
    OUTPUT_ROOT,
    PROJECT_ROOT,
    SAMPLES,
    SECRETORY_GENES,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

METHODS = ["cell2location", "tangram"]
COMPARISON_DIR = OUTPUT_ROOT / "comparison"

# CITEgeist baselines — canonical path (includes outputs/ subdir)
CITEGEIST_MORANS_PATH = PROJECT_ROOT / "mdk_analysis" / "outputs" / "if_analysis" / "bivariate_morans_results.csv"
CITEGEIST_COMMOT_DIR = PROJECT_ROOT / "output" / "commot_scores"

# Minimum secretory gene overlap to call a program "matched" (ceil(11/2) = 6)
SECRETORY_OVERLAP_THRESHOLD = 6


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------


def load_citegeist_morans() -> pd.DataFrame:
    """Load CITEgeist baseline Moran's I results.

    Raises FileNotFoundError if the canonical CSV is missing — never falls
    back to a hardcoded number.
    """
    if not CITEGEIST_MORANS_PATH.exists():
        raise FileNotFoundError(
            f"CITEgeist bivariate Moran's baseline CSV not found at "
            f"{CITEGEIST_MORANS_PATH}. Run mdk_analysis/scripts/spatial_morans.py first."
        )
    df = pd.read_csv(CITEGEIST_MORANS_PATH)
    df["method"] = "CITEgeist"
    logger.info("Loaded CITEgeist Moran's baseline: %d rows", len(df))
    return df


def load_citegeist_commot() -> pd.DataFrame:
    """Load CITEgeist COMMOT results from output/commot_scores/.

    Returns a concatenated DataFrame with columns: spot_barcode, pathway
    scores, method, sample.  Returns empty DataFrame if no files found.
    """
    frames = []
    for sample in SAMPLES:
        path = CITEGEIST_COMMOT_DIR / f"{sample}_commot_mdk.csv"
        if not path.exists():
            logger.warning("Missing CITEgeist COMMOT for %s", sample)
            continue
        df = pd.read_csv(path)
        # Rename 'barcode' to 'spot_barcode' if needed for consistency
        if "barcode" in df.columns and "spot_barcode" not in df.columns:
            df = df.rename(columns={"barcode": "spot_barcode"})
        df["method"] = "CITEgeist"
        df["sample"] = sample
        frames.append(df)

    if not frames:
        logger.warning("No CITEgeist COMMOT files found")
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    logger.info("Loaded CITEgeist COMMOT: %d rows across %d samples", len(combined), len(frames))
    return combined


# ---------------------------------------------------------------------------
# Competitor result aggregators
# ---------------------------------------------------------------------------


def aggregate_morans() -> pd.DataFrame:
    """Collect all Moran's I CSVs across methods and samples."""
    frames = []
    for method in METHODS:
        morans_dir = OUTPUT_ROOT / method / "morans"
        for sample in SAMPLES:
            path = morans_dir / f"{sample}_morans.csv"
            if not path.exists():
                logger.debug("Missing Moran's: %s / %s", method, sample)
                continue
            df = pd.read_csv(path)
            # Ensure method/sample columns present (may already be in CSV)
            df["method"] = method
            df["sample"] = sample
            frames.append(df)

    if not frames:
        logger.warning("No competitor Moran's CSVs found")
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    logger.info("Aggregated competitor Moran's: %d rows", len(combined))
    return combined


def aggregate_commot() -> pd.DataFrame:
    """Collect all COMMOT CSVs across methods and samples."""
    frames = []
    for method in METHODS:
        commot_dir = OUTPUT_ROOT / method / "commot"
        for sample in SAMPLES:
            path = commot_dir / f"{sample}_commot_mdk.csv"
            if not path.exists():
                logger.debug("Missing COMMOT: %s / %s", method, sample)
                continue
            df = pd.read_csv(path)
            df["method"] = method
            df["sample"] = sample
            frames.append(df)

    if not frames:
        logger.warning("No competitor COMMOT CSVs found")
        return pd.DataFrame()
    combined = pd.concat(frames, ignore_index=True)
    logger.info("Aggregated competitor COMMOT: %d rows", len(combined))
    return combined


def aggregate_programs() -> list[dict]:
    """Collect all program JSONs across methods and samples."""
    results = []
    for method in METHODS:
        prog_dir = OUTPUT_ROOT / method / "programs"
        for sample in SAMPLES:
            path = prog_dir / f"{sample}_programs.json"
            if not path.exists():
                logger.debug("Missing programs: %s / %s", method, sample)
                continue
            with open(path) as f:
                data = json.load(f)
            data["method"] = method
            data["sample"] = sample
            results.append(data)

    logger.info("Aggregated %d program JSONs", len(results))
    return results


# ---------------------------------------------------------------------------
# Failure-mode classifiers
# ---------------------------------------------------------------------------


def classify_morans(
    competitor_df: pd.DataFrame,
    citegeist_df: pd.DataFrame | None,
) -> list[dict]:
    """Classify Moran's I recovery for MDK per method x sample.

    Categories:
        full     -- I within +/-0.1 of CITEgeist AND p <= 0.001
        partial  -- p <= 0.05 (but not full)
        absent   -- p > 0.05 or NaN
    """
    records = []

    # Build CITEgeist per-sample MDK lookup
    cg_lookup: dict[str, float] = {}
    if citegeist_df is not None and not citegeist_df.empty:
        # CITEgeist CSV uses 'morans_I' column name
        i_col = "morans_I" if "morans_I" in citegeist_df.columns else "I"
        sample_col = "sample"
        if sample_col in citegeist_df.columns and i_col in citegeist_df.columns:
            for _, row in citegeist_df.iterrows():
                sample_key = row.get(sample_col, "")
                if sample_key and not pd.isna(row[i_col]):
                    cg_lookup[sample_key] = row[i_col]

    if competitor_df.empty:
        logger.warning("No competitor Moran's data to classify")
        return records

    # Filter to MDK rows
    mdk_mask = competitor_df["gene"].str.upper() == "MDK"
    mdk_df = competitor_df[mdk_mask]

    for method in METHODS:
        for sample in SAMPLES:
            mask = (mdk_df["method"] == method) & (mdk_df["sample"] == sample)
            subset = mdk_df[mask]

            if subset.empty:
                records.append(
                    {
                        "method": method,
                        "sample": sample,
                        "analysis": "morans_I",
                        "category": "absent",
                        "details": "no result file",
                    }
                )
                continue

            row = subset.iloc[0]
            I_val = row.get("I", np.nan)
            p_val = row.get("p_value", np.nan)

            # Baseline I for this sample (no fallback constant)
            baseline_I = cg_lookup.get(sample)

            if pd.isna(I_val) or pd.isna(p_val) or p_val > 0.05:
                category = "absent"
                details = f"I={I_val:.4f}, p={p_val:.4g}" if not pd.isna(I_val) else "NaN result"
            elif baseline_I is not None and p_val <= 0.001 and abs(I_val - baseline_I) <= 0.1:
                category = "full"
                details = f"I={I_val:.4f} (baseline={baseline_I:.4f}), p={p_val:.4g}"
            elif p_val <= 0.05:
                category = "partial"
                details = f"I={I_val:.4f}"
                if baseline_I is not None:
                    details += f" (baseline={baseline_I:.4f})"
                details += f", p={p_val:.4g}"
            else:
                category = "absent"
                details = f"I={I_val:.4f}, p={p_val:.4g}"

            records.append(
                {
                    "method": method,
                    "sample": sample,
                    "analysis": "morans_I",
                    "category": category,
                    "details": details,
                }
            )

    return records


def classify_commot(
    competitor_df: pd.DataFrame,
    citegeist_df: pd.DataFrame,
) -> list[dict]:
    """Classify COMMOT MDK-SDC4 pathway recovery per method x sample.

    Categories:
        full    -- Cancer/Epithelial is top sender AND score >= 50% of CITEgeist
        partial -- MDK-SDC4 detected but cancer not top or score < 50%
        absent  -- not detected or indistinguishable from zero
    """
    records = []
    pathway_col = "s-MDK-SDC4"
    cancer_names = {"cancer epithelial", "cancer_epithelial", "epithelial"}

    # Build CITEgeist per-sample mean MDK-SDC4 score
    cg_score_lookup: dict[str, float] = {}
    if not citegeist_df.empty and pathway_col in citegeist_df.columns:
        for sample in SAMPLES:
            s_mask = citegeist_df["sample"] == sample
            s_df = citegeist_df[s_mask]
            if not s_df.empty:
                cg_score_lookup[sample] = s_df[pathway_col].mean()

    if competitor_df.empty:
        logger.warning("No competitor COMMOT data to classify")
        for method in METHODS:
            for sample in SAMPLES:
                records.append(
                    {
                        "method": method,
                        "sample": sample,
                        "analysis": "commot_mdk_sdc4",
                        "category": "absent",
                        "details": "no result file",
                    }
                )
        return records

    for method in METHODS:
        for sample in SAMPLES:
            mask = (competitor_df["method"] == method) & (competitor_df["sample"] == sample)
            subset = competitor_df[mask]

            if subset.empty or pathway_col not in subset.columns:
                records.append(
                    {
                        "method": method,
                        "sample": sample,
                        "analysis": "commot_mdk_sdc4",
                        "category": "absent",
                        "details": "no result file",
                    }
                )
                continue

            total_score = subset[pathway_col].sum()

            # Check if indistinguishable from zero
            if total_score < 1e-6:
                records.append(
                    {
                        "method": method,
                        "sample": sample,
                        "analysis": "commot_mdk_sdc4",
                        "category": "absent",
                        "details": f"total_score={total_score:.6f} (near zero)",
                    }
                )
                continue

            # Determine top sender cell type (if cell_type column exists)
            if "cell_type" in subset.columns:
                ct_scores = subset.groupby("cell_type")[pathway_col].sum()
                top_sender = ct_scores.idxmax()
                cancer_is_top = top_sender.lower().strip() in cancer_names
            else:
                # No cell_type column -- cannot determine top sender
                cancer_is_top = False
                top_sender = "unknown"

            # Compare to CITEgeist baseline
            cg_baseline = cg_score_lookup.get(sample, None)
            mean_score = subset[pathway_col].mean()

            if cg_baseline is not None and cg_baseline > 0:
                ratio = mean_score / cg_baseline
            else:
                # No CITEgeist baseline -- classify based on detection only
                ratio = None

            if cancer_is_top and (ratio is None or ratio >= 0.5):
                category = "full"
                details = f"top_sender={top_sender}, mean={mean_score:.4f}"
                if ratio is not None:
                    details += f", ratio_to_cg={ratio:.2f}"
            elif total_score > 1e-6:
                category = "partial"
                details = f"top_sender={top_sender}, mean={mean_score:.4f}"
                if ratio is not None:
                    details += f", ratio_to_cg={ratio:.2f}"
                if not cancer_is_top:
                    details += " (cancer not top sender)"
                if ratio is not None and ratio < 0.5:
                    details += " (score <50% of CITEgeist)"
            else:
                category = "absent"
                details = "not detected"

            records.append(
                {
                    "method": method,
                    "sample": sample,
                    "analysis": "commot_mdk_sdc4",
                    "category": category,
                    "details": details,
                }
            )

    return records


def classify_programs(programs_list: list[dict]) -> list[dict]:
    """Classify NMF program recovery of MDK secretory program in cancer cells.

    A matched secretory program requires BOTH:
      (a) MDK in the program's top-50 genes, AND
      (b) >= SECRETORY_OVERLAP_THRESHOLD (6) of the 11 SECRETORY_GENES
          in the program's top-50 genes.

    Categories:
        full    -- MDK in top-50 AND >=6 secretory genes in top-50
                   (matched secretory program)
        partial -- MDK in top-50 but <6 secretory genes
                   (MDK detected but not a secretory program)
        absent  -- MDK not in any cancer program top-50
    """
    records = []
    cancer_names = {"cancer epithelial", "cancer_epithelial", "epithelial"}
    sec_set = set(SECRETORY_GENES)

    # Index by (method, sample) for lookup
    prog_lookup: dict[tuple[str, str], dict] = {}
    for entry in programs_list:
        key = (entry["method"], entry["sample"])
        prog_lookup[key] = entry

    for method in METHODS:
        for sample in SAMPLES:
            key = (method, sample)
            if key not in prog_lookup:
                records.append(
                    {
                        "method": method,
                        "sample": sample,
                        "analysis": "programs_mdk",
                        "category": "absent",
                        "details": "no result file",
                    }
                )
                continue

            entry = prog_lookup[key]
            programs_dict = entry.get("programs", {})

            # Find cancer cell type programs
            cancer_programs = []
            for ct, progs in programs_dict.items():
                if ct.lower().strip() in cancer_names:
                    cancer_programs.extend(progs)

            if not cancer_programs:
                records.append(
                    {
                        "method": method,
                        "sample": sample,
                        "analysis": "programs_mdk",
                        "category": "absent",
                        "details": "no cancer cell type programs found",
                    }
                )
                continue

            # Check each cancer program for matched secretory criterion
            best_category = "absent"
            best_details = "MDK not in any cancer program top-50"

            for p in cancer_programs:
                mdk_in_top50 = p.get("mdk_in_top50", False)
                if not mdk_in_top50:
                    continue

                # MDK is in top-50 -- check secretory overlap
                top_50_genes = p.get("top_50_genes", [])
                if not top_50_genes:
                    # Fallback: top_50_genes not saved (old format)
                    logger.warning(
                        "Program %s/%s prog %s: top_50_genes not saved, " "cannot check secretory overlap",
                        method,
                        sample,
                        p.get("program_idx", "?"),
                    )
                    if best_category == "absent":
                        best_category = "partial"
                        best_details = (
                            f"MDK in top-50 of program {p.get('program_idx', '?')} "
                            f"but top_50_genes not available for secretory check"
                        )
                    continue

                top_50_set = set(top_50_genes)
                sec_overlap = sec_set & top_50_set
                n_sec = len(sec_overlap)

                if n_sec >= SECRETORY_OVERLAP_THRESHOLD:
                    # Full match: MDK + secretory program
                    best_category = "full"
                    best_details = (
                        f"MDK + {n_sec}/{len(SECRETORY_GENES)} secretory genes "
                        f"in top-50 of program {p.get('program_idx', '?')} "
                        f"(genes: {sorted(sec_overlap)})"
                    )
                    break  # No need to check further programs
                else:
                    # Partial: MDK found but not enough secretory genes
                    if best_category != "full":
                        best_category = "partial"
                        best_details = (
                            f"MDK in top-50 of program {p.get('program_idx', '?')} "
                            f"but only {n_sec}/{len(SECRETORY_GENES)} secretory genes "
                            f"(need >={SECRETORY_OVERLAP_THRESHOLD})"
                        )

            records.append(
                {
                    "method": method,
                    "sample": sample,
                    "analysis": "programs_mdk",
                    "category": best_category,
                    "details": best_details,
                }
            )

    return records


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Aggregate all results, classify failure modes, and save comparison tables."""
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load CITEgeist baselines ---
    citegeist_morans = load_citegeist_morans()
    citegeist_commot = load_citegeist_commot()

    # --- Aggregate competitor results ---
    competitor_morans = aggregate_morans()
    competitor_commot = aggregate_commot()
    competitor_programs = aggregate_programs()

    # --- Save comparison CSVs ---
    # Moran's comparison: merge CITEgeist + competitors for candidate genes
    morans_frames = []
    if citegeist_morans is not None and not citegeist_morans.empty:
        morans_frames.append(citegeist_morans)
    if not competitor_morans.empty:
        morans_frames.append(competitor_morans)

    if morans_frames:
        all_morans = pd.concat(morans_frames, ignore_index=True)
        # Filter to candidate genes if column exists
        if "gene" in all_morans.columns:
            gene_mask = all_morans["gene"].str.upper().isin([g.upper() for g in CANDIDATE_GENES])
            all_morans = all_morans[gene_mask]
        out_path = COMPARISON_DIR / "morans_comparison.csv"
        all_morans.to_csv(out_path, index=False)
        logger.info("Saved %s (%d rows)", out_path, len(all_morans))
    else:
        logger.warning("No Moran's data to save")

    # COMMOT comparison
    commot_frames = []
    if not citegeist_commot.empty:
        commot_frames.append(citegeist_commot)
    if not competitor_commot.empty:
        commot_frames.append(competitor_commot)

    if commot_frames:
        all_commot = pd.concat(commot_frames, ignore_index=True)
        out_path = COMPARISON_DIR / "commot_comparison.csv"
        all_commot.to_csv(out_path, index=False)
        logger.info("Saved %s (%d rows)", out_path, len(all_commot))
    else:
        logger.warning("No COMMOT data to save")

    # --- Classify failure modes ---
    failure_records = []
    failure_records.extend(classify_morans(competitor_morans, citegeist_morans))
    failure_records.extend(classify_commot(competitor_commot, citegeist_commot))
    failure_records.extend(classify_programs(competitor_programs))

    if failure_records:
        failure_df = pd.DataFrame(failure_records)
        out_path = COMPARISON_DIR / "failure_mode_table.csv"
        failure_df.to_csv(out_path, index=False)
        logger.info("Saved %s (%d rows)", out_path, len(failure_df))

        # Pivot summary: method x analysis x category -> count
        summary = failure_df.groupby(["method", "analysis", "category"]).size().reset_index(name="count")
        # Pivot wider for readability
        summary_wide = summary.pivot_table(
            index=["method", "analysis"],
            columns="category",
            values="count",
            fill_value=0,
        ).reset_index()
        out_path = COMPARISON_DIR / "failure_mode_summary.csv"
        summary_wide.to_csv(out_path, index=False)
        logger.info("Saved %s", out_path)

        # Print summary to stdout
        print("\n" + "=" * 70)
        print("FAILURE MODE SUMMARY")
        print("=" * 70)
        for method in METHODS:
            print(f"\n--- {method} ---")
            m_df = failure_df[failure_df["method"] == method]
            for analysis in sorted(m_df["analysis"].unique()):
                a_df = m_df[m_df["analysis"] == analysis]
                counts = a_df["category"].value_counts()
                parts = [f"{cat}: {n}/{len(a_df)}" for cat, n in counts.items()]
                print(f"  {analysis}: {', '.join(parts)}")
        print("=" * 70)
    else:
        logger.warning("No failure mode classifications generated (no competitor results found)")


if __name__ == "__main__":
    main()
