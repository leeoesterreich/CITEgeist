"""
Diagnostic: GT calling validity and interaction with Module 2 / Module 3.5.

Tests whether the current quantile(0.8) GT approach produces valid binary labels
for sparse protein markers, and compares with GMM-based calling.

Examines:
  1. Per-marker protein distribution in each target cell type
  2. Zero-inflation rates (cause of single_class_gt collapse)
  3. Quantile(0.8) threshold behavior on zero-inflated data
  4. 2-component GMM alternative on cell-level protein data
  5. Module 3.5 spot-level ratio signal for the same pairs
  6. Module 2 marker interest signals (kurtosis, GMM SNR, Moran's I)

Run on login node (CPU-only, read-only).
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))
import _paths as _cg  # noqa: E402

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture

REPO_ROOT = Path(str(_cg.REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking"))
sys.path.insert(0, str(REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "src"))

from load_xenium_singlecell import load_xenium_singlecell

BENCHMARK_DIR = REPO_ROOT / "Benchmarking" / "xenium_benchmarking" / "CITEgeist" / "output_module3_5_benchmark_array"
PSEUDOVISIUM_DIR = REPO_ROOT / "Benchmarking" / "xenium_pseudovisium" / "data_protein_gt"

FUNCTIONAL_MARKERS = ["PD-1", "PD-L1", "PCNA", "Vimentin"]
ACTIVE_TYPES = {
    "PD-1": ["T cells"],
    "PD-L1": ["Macrophages", "Epithelial"],
    "PCNA": ["B cells", "T cells", "Macrophages", "Endothelial", "Epithelial", "Fibroblasts"],
    "Vimentin": ["Fibroblasts", "Endothelial", "Epithelial"],
}
GT_TO_6TYPE = {
    "B cells": "B cells",
    "CD4+ T cells": "T cells",
    "CD8+ T cells": "T cells",
    "Macrophages": "Macrophages",
    "Endothelial": "Endothelial",
    "Epithelial": "Epithelial",
    "Fibroblasts": "Fibroblasts",
}

REGION_ID = 0  # Examine region 0; key findings generalize


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def gmm_binary_calls(values: np.ndarray, random_state: int = 42) -> tuple[np.ndarray, dict]:
    """Fit 2-component GMM on log1p values; return binary calls and diagnostics."""
    log_vals = np.log1p(values).reshape(-1, 1)
    info = {}

    # Fraction zero
    frac_zero = float((values == 0).mean())
    info["frac_zero"] = frac_zero

    # Fraction positive above log1p(1) = 0.693
    frac_nonzero = float((values > 0).mean())
    info["frac_nonzero"] = frac_nonzero

    try:
        gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=random_state, n_init=5)
        gmm.fit(log_vals)

        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        pooled_std = np.sqrt(gmm.weights_[0] * stds[0] ** 2 + gmm.weights_[1] * stds[1] ** 2)
        separation = abs(means[0] - means[1])

        high_comp = int(np.argmax(means))
        posteriors = gmm.predict_proba(log_vals)
        gmm_calls = (posteriors[:, high_comp] > 0.5).astype(int)

        info["gmm_means"] = means.tolist()
        info["gmm_weights"] = gmm.weights_.tolist()
        info["gmm_separation"] = float(separation)
        info["gmm_pooled_std"] = float(pooled_std)
        info["gmm_separation_ratio"] = float(separation / (pooled_std + 1e-8))
        info["gmm_high_comp_weight"] = float(gmm.weights_[high_comp])
        info["gmm_bimodal"] = separation > 0.5 * pooled_std
        info["gmm_frac_positive"] = float(gmm_calls.mean())
        info["status"] = "ok"
        return gmm_calls, info

    except Exception as exc:  # pylint: disable=broad-except
        info["status"] = f"failed: {exc}"
        return np.zeros(len(values), dtype=int), info


def quantile_binary_calls(values: np.ndarray, q: float = 0.8) -> tuple[np.ndarray, dict]:
    """Current benchmark approach: quantile threshold on log1p values."""
    log_vals = np.log1p(values.astype(float))
    threshold = float(np.quantile(log_vals, q))
    calls = (log_vals >= threshold).astype(int)
    return calls, {
        "threshold": threshold,
        "frac_positive": float(calls.mean()),
        "q_is_zero": threshold == 0.0,
    }


# ---------------------------------------------------------------------------
# Section 1: Cell-level protein distribution analysis
# ---------------------------------------------------------------------------


def section1_protein_distributions(protein_df: pd.DataFrame, cell_types: pd.Series) -> pd.DataFrame:
    """Characterize per-(cell_type, marker) protein distributions."""
    rows = []
    for marker in FUNCTIONAL_MARKERS:
        if marker not in protein_df.columns:
            print(f"  WARNING: marker {marker} not in protein data, skipping")
            continue
        for ct in ACTIVE_TYPES[marker]:
            mask = cell_types == ct
            n_cells = int(mask.sum())
            if n_cells < 10:
                continue
            values = protein_df.loc[mask, marker].values.astype(float)

            frac_zero = float((values == 0).mean())
            mean_expr = float(values.mean())
            median_expr = float(np.median(values))
            q80 = float(np.quantile(np.log1p(values), 0.8))
            max_expr = float(values.max())

            q_calls, q_info = quantile_binary_calls(values)
            gmm_calls, gmm_info = gmm_binary_calls(values)

            rows.append(
                {
                    "pair_id": f"{ct}__{marker}",
                    "cell_type": ct,
                    "marker": marker,
                    "n_cells": n_cells,
                    "frac_zero": frac_zero,
                    "mean_raw": mean_expr,
                    "median_raw": median_expr,
                    "max_raw": max_expr,
                    "log1p_q80": q80,
                    "q80_is_zero": q_info["q_is_zero"],
                    "quantile_frac_positive": q_info["frac_positive"],
                    "gmm_bimodal": gmm_info.get("gmm_bimodal", False),
                    "gmm_separation_ratio": gmm_info.get("gmm_separation_ratio", np.nan),
                    "gmm_frac_positive": gmm_info.get("gmm_frac_positive", np.nan),
                    "gmm_status": gmm_info.get("status", "unknown"),
                    "gmm_high_comp_weight": gmm_info.get("gmm_high_comp_weight", np.nan),
                }
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Section 2: Spot-level ratio signal from Module 3.5 output
# ---------------------------------------------------------------------------


def section2_spot_ratio_signal(region_id: int) -> pd.DataFrame:
    """Load Module 3.5 projected functional calls and examine ratio distributions."""
    parquet_path = BENCHMARK_DIR / f"Xenium_region_{region_id}" / "projected_functional_calls.parquet"
    if not parquet_path.exists():
        print(f"  No projected calls at {parquet_path}")
        return pd.DataFrame()

    proj = pd.read_parquet(parquet_path)

    rows = []
    for marker in FUNCTIONAL_MARKERS:
        for ct in ACTIVE_TYPES.get(marker, []):
            score_col = f"func_{marker}_{ct}_score"
            call_col = f"func_{marker}_{ct}_call"
            if score_col not in proj.columns:
                continue

            scores = proj[score_col].dropna()
            calls = proj[call_col].dropna() if call_col in proj.columns else pd.Series(dtype=float)

            rows.append(
                {
                    "pair_id": f"{ct}__{marker}",
                    "cell_type": ct,
                    "marker": marker,
                    "n_projected_cells": len(scores),
                    "mean_score": float(scores.mean()),
                    "median_score": float(scores.median()),
                    "frac_positive_predicted": float(calls.mean()) if len(calls) > 0 else np.nan,
                    "score_p25": float(scores.quantile(0.25)),
                    "score_p75": float(scores.quantile(0.75)),
                    "score_p90": float(scores.quantile(0.90)),
                }
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Section 3: Module 2 marker interest signal
# ---------------------------------------------------------------------------


def section3_module2_signals() -> pd.DataFrame:
    """
    Load the pseudo-Visium spot-level protein data and compute Module 1/2-style
    marker interest signals for the functional markers.
    Shows whether Module 2 would flag these markers as interesting.
    """
    cite_path = PSEUDOVISIUM_DIR / "h5ad_objects" / f"Xenium_region_{REGION_ID}_CITE.h5ad"
    if not cite_path.exists():
        print(f"  No CITE h5ad at {cite_path}")
        return pd.DataFrame()

    adata = sc.read_h5ad(cite_path)
    ab_df = adata.to_df()

    rows = []
    for marker in FUNCTIONAL_MARKERS:
        if marker not in ab_df.columns:
            # try suffix version
            candidates = [c for c in ab_df.columns if c.startswith(marker)]
            if not candidates:
                print(f"  marker {marker} not in CITE data, columns: {list(ab_df.columns)[:10]}...")
                continue
            marker_col = candidates[0]
        else:
            marker_col = marker

        values = ab_df[marker_col].values.astype(float)

        # Kurtosis
        from scipy.stats import kurtosis

        kurt = float(kurtosis(values, fisher=True))

        # Fraction zero at spot level
        frac_zero_spots = float((values == 0).mean())

        # 2-comp GMM at spot level (same as Module 1)
        log_vals = np.log1p(values)
        gmm_calls, gmm_info = gmm_binary_calls(values)

        # Coefficient of variation
        cv = float(np.std(values) / (np.mean(values) + 1e-8))

        rows.append(
            {
                "marker": marker,
                "n_spots": len(values),
                "mean_spot_expr": float(values.mean()),
                "frac_zero_spots": frac_zero_spots,
                "kurtosis": kurt,
                "kurtosis_gate": kurt > 2.0,
                "cv": cv,
                "gmm_bimodal": gmm_info.get("gmm_bimodal", False),
                "gmm_separation_ratio": gmm_info.get("gmm_separation_ratio", np.nan),
                "gmm_frac_positive": gmm_info.get("gmm_frac_positive", np.nan),
            }
        )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    print("=" * 70)
    print("Diagnostic: GT calling validity and Module 2 / 3.5 interaction")
    print(f"Region: {REGION_ID}")
    print("=" * 70)

    # Load single-cell protein data
    print("\nLoading Xenium single-cell data for region 0...")
    _, adata_protein = load_xenium_singlecell(region_id=REGION_ID)

    # Load cell type assignments
    mapping = pd.read_csv(PSEUDOVISIUM_DIR / "cell_to_spot_mapping.csv", index_col=0)
    labels = pd.read_csv(PSEUDOVISIUM_DIR / "cell_type_assignments.csv", index_col=0)
    mapping = mapping[(mapping["region_id"] == REGION_ID) & (mapping["spot_idx"] != -1)].copy()
    mapping = mapping.join(labels, how="left")
    mapping = mapping[mapping["cell_type"].isin(GT_TO_6TYPE)].copy()
    mapping["cell_type_6"] = mapping["cell_type"].map(GT_TO_6TYPE)

    protein_df = adata_protein.to_df()
    protein_df["cell_type"] = mapping["cell_type_6"].reindex(protein_df.index)
    protein_df = protein_df.dropna(subset=["cell_type"]).copy()
    cell_types = protein_df["cell_type"]

    print(f"  {len(protein_df):,} cells with GT labels")
    print(f"  Cell type distribution: {cell_types.value_counts().to_dict()}")
    print(f"  Protein columns: {[c for c in protein_df.columns if c != 'cell_type']}")

    # Section 1
    print("\n" + "=" * 70)
    print("SECTION 1: Cell-level protein distribution per (cell_type, marker)")
    print("=" * 70)
    dist_df = section1_protein_distributions(protein_df, cell_types)
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    pd.set_option("display.float_format", "{:.3f}".format)
    print(
        dist_df[
            [
                "pair_id",
                "n_cells",
                "frac_zero",
                "log1p_q80",
                "q80_is_zero",
                "quantile_frac_positive",
                "gmm_bimodal",
                "gmm_separation_ratio",
                "gmm_frac_positive",
            ]
        ].to_string(index=False)
    )

    print("\n--- KEY FINDING for single_class_gt pairs ---")
    single_class = dist_df[dist_df["q80_is_zero"]]
    if len(single_class) > 0:
        print("Pairs where quantile(0.8) threshold collapses to 0 → 100% positive GT:")
        for _, row in single_class.iterrows():
            print(
                f"  {row['pair_id']}: frac_zero={row['frac_zero']:.3f}, "
                f"quantile_frac_pos={row['quantile_frac_positive']:.3f}, "
                f"gmm_bimodal={row['gmm_bimodal']}, "
                f"gmm_frac_pos={row['gmm_frac_positive']:.3f}"
            )
    else:
        print("No pairs with collapsed threshold found for this region.")

    # Section 2
    print("\n" + "=" * 70)
    print("SECTION 2: Module 3.5 spot-level ratio signal (predicted scores)")
    print("=" * 70)
    spot_df = section2_spot_ratio_signal(REGION_ID)
    if not spot_df.empty:
        print(spot_df.to_string(index=False))
    else:
        print("  No projected functional calls found. Run the benchmark first.")

    # Section 3
    print("\n" + "=" * 70)
    print("SECTION 3: Module 2 marker interest signals (spot-level)")
    print("=" * 70)
    mod2_df = section3_module2_signals()
    if not mod2_df.empty:
        print(
            mod2_df[
                [
                    "marker",
                    "n_spots",
                    "mean_spot_expr",
                    "frac_zero_spots",
                    "kurtosis",
                    "kurtosis_gate",
                    "gmm_bimodal",
                    "gmm_separation_ratio",
                ]
            ].to_string(index=False)
        )
    else:
        print("  Could not load CITE spot data.")

    # Section 4: Cross-level consistency summary
    print("\n" + "=" * 70)
    print("SECTION 4: Cross-level consistency summary")
    print("=" * 70)
    if not dist_df.empty and not mod2_df.empty:
        merged = dist_df.merge(
            mod2_df[["marker", "kurtosis_gate", "gmm_bimodal"]], on="marker", how="left", suffixes=("_cell", "_spot")
        )
        merged = merged.rename(
            columns={
                "kurtosis_gate": "spot_kurtosis_gate",
                "gmm_bimodal_spot": "spot_gmm_bimodal",
                "gmm_bimodal_cell": "cell_gmm_bimodal",
            }
        )
        # Highlight mismatches
        print("\nPer-pair: cell-level GMM bimodality vs spot-level GMM bimodality")
        for _, row in merged.iterrows():
            cell_bio = row.get("cell_gmm_bimodal", "?")
            spot_bio = row.get("spot_gmm_bimodal", "?")
            mismatch = " *** MISMATCH" if cell_bio != spot_bio else ""
            print(
                f"  {row['pair_id']:<30}  cell_bimodal={cell_bio}  "
                f"spot_bimodal={spot_bio}  q80_zero={row['q80_is_zero']}{mismatch}"
            )

    print("\n" + "=" * 70)
    print("DIAGNOSTIC COMPLETE")
    print("=" * 70)
    print("""
SUMMARY OF ISSUES TO INVESTIGATE:
  1. Quantile collapse: if frac_zero > 0.8 for a pair, quantile(0.8)=0 → all cells positive
  2. GMM validity at cell level: does the single-cell protein data show bimodal structure?
  3. Level mismatch: Module 3.5 gates at SPOT level, but its Module 2 candidates
     were generated from single cells; GT calls at CELL level
  4. If spot-level signals are bimodal but cell-level is unimodal → spot aggregation is key
""")


if __name__ == "__main__":
    main()
