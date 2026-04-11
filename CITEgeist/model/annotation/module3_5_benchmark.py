from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    precision_score,
    recall_score,
    roc_auc_score,
)


@dataclass(frozen=True)
class PairBenchmarkResult:
    pair_id: str
    cell_type: str
    marker: str
    n_positive: int
    n_negative: int
    n_supporting_spots: int
    auroc: float
    auprc: float
    balanced_accuracy: float
    precision: float
    recall: float
    gt_positive_fraction: float
    pred_positive_fraction: float
    evaluable: bool
    exclusion_reason: str | None
    headline: bool = False


def build_gt_binary_calls(
    cell_table: pd.DataFrame,
    cell_type: str,
    marker: str,
) -> pd.Series:
    """Build binary ground-truth calls using a 2-component GMM on log1p protein values.

    Fits a GMM to identify the high-expression component.  Falls back to
    nonzero gating when the distribution is unimodal (separation < 0.5 * pooled std),
    which handles ubiquitous markers correctly.

    The previous quantile(0.8) approach collapsed to threshold=0 for sparse markers
    (>80% zeros), labelling every cell positive (single_class_gt), and forced an
    artificial 20/80 split for ubiquitous markers like Vimentin where the biology is
    >95% positive.
    """
    from sklearn.mixture import GaussianMixture  # pylint: disable=import-outside-toplevel

    required_columns = {"cell_type", marker}
    missing_columns = required_columns - set(cell_table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise KeyError(f"Missing required columns: {missing}")

    subset = cell_table.loc[cell_table["cell_type"] == cell_type, [marker]].copy()
    if subset.empty:
        return pd.Series(dtype=int, name=marker)

    values = subset[marker].astype(float).values
    log_vals = np.log1p(values).reshape(-1, 1)

    try:
        gmm = GaussianMixture(n_components=2, covariance_type="full", random_state=42, n_init=5)
        gmm.fit(log_vals)

        means = gmm.means_.flatten()
        stds = np.sqrt(gmm.covariances_.flatten())
        pooled_std = np.sqrt(gmm.weights_[0] * stds[0] ** 2 + gmm.weights_[1] * stds[1] ** 2)
        separation = abs(means[0] - means[1])
        high_comp = int(np.argmax(means))

        if separation > 0.5 * pooled_std:
            posteriors = gmm.predict_proba(log_vals)
            calls_arr = (posteriors[:, high_comp] > 0.5).astype(int)
        else:
            # Unimodal / flat — fall back to nonzero gating
            calls_arr = (values > 0).astype(int)

    except Exception:  # pylint: disable=broad-except
        # Safe fallback: nonzero
        calls_arr = (values > 0).astype(int)

    calls: pd.Series = pd.Series(calls_arr, index=subset.index, name=f"gt_{marker}")
    return calls


def score_pair_predictions(
    pair_id: str,
    cell_type: str,
    marker: str,
    gt_calls: pd.Series,
    pred_scores: pd.Series,
    *,
    pred_calls: pd.Series,
    n_supporting_spots: int,
    headline: bool,
) -> PairBenchmarkResult:
    """Score a functional pair against ground-truth and predicted calls."""
    common_index = gt_calls.index.intersection(pred_scores.index).intersection(pred_calls.index)
    if len(common_index) == 0:
        return PairBenchmarkResult(
            pair_id=pair_id,
            cell_type=cell_type,
            marker=marker,
            n_positive=0,
            n_negative=0,
            n_supporting_spots=n_supporting_spots,
            auroc=0.0,
            auprc=0.0,
            balanced_accuracy=0.0,
            precision=0.0,
            recall=0.0,
            gt_positive_fraction=0.0,
            pred_positive_fraction=0.0,
            evaluable=False,
            exclusion_reason="no_shared_cells",
            headline=headline,
        )

    gt = gt_calls.loc[common_index].astype(int)
    scores = pred_scores.loc[common_index].astype(float)
    calls = pred_calls.loc[common_index].astype(int)

    n_positive = int(gt.sum())
    n_negative = int(len(gt) - n_positive)
    gt_positive_fraction = float(gt.mean()) if len(gt) else 0.0
    pred_positive_fraction = float(calls.mean()) if len(calls) else 0.0

    evaluable = bool(len(common_index) >= 4 and n_positive > 0 and n_negative > 0 and n_supporting_spots >= 20)

    if not evaluable:
        exclusion_reason = "insufficient_support"
        if n_positive == 0 or n_negative == 0:
            exclusion_reason = "single_class_gt"

        return PairBenchmarkResult(
            pair_id=pair_id,
            cell_type=cell_type,
            marker=marker,
            n_positive=n_positive,
            n_negative=n_negative,
            n_supporting_spots=n_supporting_spots,
            auroc=0.0,
            auprc=0.0,
            balanced_accuracy=0.0,
            precision=0.0,
            recall=0.0,
            gt_positive_fraction=gt_positive_fraction,
            pred_positive_fraction=pred_positive_fraction,
            evaluable=False,
            exclusion_reason=exclusion_reason,
            headline=headline,
        )

    auroc = float(roc_auc_score(gt, scores))
    auprc = float(average_precision_score(gt, scores))
    balanced_acc = float(balanced_accuracy_score(gt, calls))
    precision = float(precision_score(gt, calls, zero_division=0))
    recall = float(recall_score(gt, calls, zero_division=0))

    return PairBenchmarkResult(
        pair_id=pair_id,
        cell_type=cell_type,
        marker=marker,
        n_positive=n_positive,
        n_negative=n_negative,
        n_supporting_spots=n_supporting_spots,
        auroc=auroc,
        auprc=auprc,
        balanced_accuracy=balanced_acc,
        precision=precision,
        recall=recall,
        gt_positive_fraction=gt_positive_fraction,
        pred_positive_fraction=pred_positive_fraction,
        evaluable=True,
        exclusion_reason=None,
        headline=headline,
    )


def _pair_validated(pair: PairBenchmarkResult, thresholds: dict[str, float]) -> bool:
    if not pair.evaluable:
        return False
    return (
        pair.auroc >= thresholds["auroc"]
        and pair.auprc >= thresholds["auprc"]
        and pair.balanced_accuracy >= thresholds["balanced_accuracy"]
    )


def aggregate_module3_5_results(
    pairs: list[PairBenchmarkResult],
    thresholds: dict[str, float],
) -> dict[str, Any]:
    pair_rows: list[dict[str, Any]] = []
    validated_pairs: list[str] = []
    evaluable_pairs: list[str] = []
    headline_failures: list[str] = []

    for pair in pairs:
        validated = _pair_validated(pair, thresholds)
        row = asdict(pair)
        row["validated"] = validated
        pair_rows.append(row)

        if pair.evaluable:
            evaluable_pairs.append(pair.pair_id)
        if validated:
            validated_pairs.append(pair.pair_id)
        if pair.headline and pair.evaluable and not validated:
            headline_failures.append(pair.pair_id)

    majority_pass = len(validated_pairs) >= ((len(evaluable_pairs) // 2) + 1) if evaluable_pairs else False
    benchmark_passed = majority_pass and not headline_failures

    return {
        "thresholds": thresholds,
        "pairs": pair_rows,
        "n_pairs": len(pairs),
        "n_evaluable_pairs": len(evaluable_pairs),
        "validated_pairs": validated_pairs,
        "headline_failures": headline_failures,
        "benchmark_passed": benchmark_passed,
    }


def score_spot_attribution(  # pylint: disable=too-many-positional-arguments
    cell_type: str,
    marker: str,
    lambda_df: pd.DataFrame,
    proportions_df: pd.DataFrame,
    protein_df: pd.DataFrame,
    mapping: pd.DataFrame,
) -> dict | None:
    """Compute SACE-protein spot attribution accuracy in GT-positive spots.

    For spots containing ≥1 GT-positive cell of `cell_type`, measures the
    fraction of the marker signal attributed to `cell_type` by SACE-protein λ:

        share[spot, t] = p[spot, t] * λ[t, marker] / Σ_t' p[spot, t'] * λ[t', marker]

    This is distinct from cell-level AUROC — it asks "in spots where the signal
    exists, does λ correctly identify which type drives it?" which is robust to
    marker sparsity.

    Returns None when λ or GT data is unavailable for this pair, or a dict with:
        n_gt_positive_cells  – GT-positive cells of this type
        n_hot_spots          – spots containing ≥1 GT-positive cell
        mean_attribution     – mean attribution share for target type in hot spots
        baseline_attribution – mean cell-type proportion in hot spots (λ-blind baseline)
        fraction_dominant    – fraction of hot spots where target type has highest share
    """
    # Require non-zero λ for this (type, marker) pair
    if cell_type not in lambda_df.index or marker not in lambda_df.columns:
        return None
    if float(lambda_df.at[cell_type, marker]) == 0.0:
        return None

    # GT-positive cells of this type using GMM gating
    try:
        gt_calls = build_gt_binary_calls(protein_df, cell_type=cell_type, marker=marker)
    except (KeyError, ValueError, RuntimeError):
        return None

    n_positive = int(gt_calls.sum())
    if n_positive == 0:
        return {
            "n_gt_positive_cells": 0,
            "n_hot_spots": 0,
            "mean_attribution": None,
            "baseline_attribution": None,
            "fraction_dominant": None,
        }

    # Map positive cells to spots
    positive_cell_ids = gt_calls.index[gt_calls == 1]
    valid_cells = mapping.index.intersection(positive_cell_ids)
    if valid_cells.empty:
        return {
            "n_gt_positive_cells": n_positive,
            "n_hot_spots": 0,
            "mean_attribution": None,
            "baseline_attribution": None,
            "fraction_dominant": None,
        }

    hot_spot_ids = mapping.loc[valid_cells, "spot_id"].unique()
    hot_spot_ids_in_prop = [s for s in hot_spot_ids if s in proportions_df.index]
    if not hot_spot_ids_in_prop:
        return {
            "n_gt_positive_cells": n_positive,
            "n_hot_spots": int(len(hot_spot_ids)),
            "mean_attribution": None,
            "baseline_attribution": None,
            "fraction_dominant": None,
            "note": "hot spots not in proportions index",
        }

    prop_hot = proportions_df.loc[hot_spot_ids_in_prop]  # (H, T)

    # λ vector across all types for this marker
    lambda_vec = lambda_df[marker].reindex(proportions_df.columns).fillna(0.0).values

    # Expected signal contribution per type per spot: p[spot,t] * λ[t,m]
    weighted = prop_hot.values * lambda_vec[np.newaxis, :]  # (H, T)
    total_signal = weighted.sum(axis=1, keepdims=True)
    total_signal = np.maximum(total_signal, 1e-12)
    attribution_share = weighted / total_signal  # (H, T) normalized

    if cell_type not in proportions_df.columns:
        return None
    type_idx = list(proportions_df.columns).index(cell_type)

    target_attribution = attribution_share[:, type_idx]
    baseline_prop = prop_hot.iloc[:, type_idx].values
    fraction_dominant = float((attribution_share.argmax(axis=1) == type_idx).mean())

    return {
        "n_gt_positive_cells": n_positive,
        "n_hot_spots": len(hot_spot_ids_in_prop),
        "mean_attribution": float(target_attribution.mean()),
        "baseline_attribution": float(baseline_prop.mean()),
        "fraction_dominant": fraction_dominant,
    }
