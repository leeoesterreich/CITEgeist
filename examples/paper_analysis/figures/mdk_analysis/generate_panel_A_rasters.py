#!/usr/bin/env python -u
"""Generate Fig 4 Panel A: spatial maps for two progressing patients (P1, P4).

Layout: 2 rows x 4 columns (cancer proportion / late E2 directional response).
Uses SACE cancer-layer spot-level GEX (not per-cell assignment) so that
low-cellularity biopsies like P4-S1 are evaluable.

E2 scoring: directional UP − DN contrast of MSigDB Li 2023
LI_ESTROGENE_LATE_E2_RESPONSE_{UP,DN} on library-size-normalized,
log1p-transformed SACE cancer-compartment counts per spot.
Spots with cancer proportion ≤ 0.1 or total cancer-layer counts < 50
are NaN-filled (rendered as translucent gray).

H&E background via squidpy spatial_scatter with crop_coord to spot extent.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from matplotlib.colors import Normalize

SCRIPT_DIR = Path(__file__).resolve().parent

# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
UNIFIED_DIR = PROJECT / "CITEgeist/output/unified_pipeline"
MODULE3_DIR = PROJECT / "output/module3_unified"
OUT_DIR = PROJECT / "mdk_analysis/outputs/if_analysis"
QC_DIR = PROJECT / "figures/mdk_analysis/outputs"
QC_PANEL_DIR = QC_DIR / "panel_A_qc"

sys.path.insert(0, str(SCRIPT_DIR.parent))
sys.path.insert(0, str(SCRIPT_DIR))
import _e2_scoring as e2  # noqa: E402
from _shared.spatial_utils import inject_spatial_metadata  # noqa: E402

SAMPLES: list[tuple[str, str, str, str, str]] = [
    # (m3_basename, specimen_id, display_label, he_specimen, unified_specimen)
    (
        "HCC22-088-P1-S1_module3_results.h5ad",
        "HCC22-088-P1-S1",
        "P1-S1\n(biopsy)",
        "HCC22-088-P1-S1",
        "HCC22-088-P1-S1",
    ),
    (
        "HCC22-088-P1-S2_module3_results.h5ad",
        "HCC22-088-P1-S2",
        "P1-S2\n(surgical)",
        "HCC22-088-P1-S2",
        "HCC22-088-P1-S2",
    ),
    (
        "HCC22-088-P4-S1_module3_results.h5ad",
        "HCC22-088-P4-S1",
        "P4-S1\n(biopsy)",
        "HCC22-088-P4-S1",
        "HCC22-088-P4-S1",
    ),
    (
        "HCC22-088-P4-S2_1i_rep_module3_results.h5ad",
        "HCC22-088-P4-S2_1i_rep",
        "P4-S2\n(surgical, D538G)",
        "HCC22-088-P4-S2_1i_rep",
        "HCC22-088-P4-S2_1i_rep",
    ),
]

CANCER_PROP_THRESHOLD = 0.1
MIN_CANCER_LAYER_COUNTS = 50

CBAR_MARGIN_IN = 0.6
FONT_TITLE = 18
FONT_ROWLABEL = 18
FONT_CBAR_LABEL = 16
FONT_CBAR_TICK = 14
SQ_SPOT_SIZE = 1.3
NA_RGBA = (0.83, 0.83, 0.83, 0.4)


def _spot_crop(adata, pad_frac=0.03):
    c = np.asarray(adata.obsm["spatial"], dtype=float)
    xn, xx = float(c[:, 0].min()), float(c[:, 0].max())
    yn, yx = float(c[:, 1].min()), float(c[:, 1].max())
    px, py = (xx - xn) * pad_frac, (yx - yn) * pad_frac
    return (max(0, int(xn - px)), max(0, int(yn - py)), int(xx + px), int(yx + py))


def _load_cancer_layer(specimen: str) -> tuple[ad.AnnData, pd.Series]:
    """Load SACE cancer-layer CSV into AnnData, normalize, log1p.

    Returns (normalized_adata, raw_total_counts_per_spot).
    """
    layer_path = (
        UNIFIED_DIR / specimen / "module3" / f"{specimen}_pass1" / "layers" / "pass1" / "Cancer_layer_pass1.csv"
    )
    if not layer_path.exists():
        raise FileNotFoundError(f"Cancer layer not found: {layer_path}")

    df = pd.read_csv(layer_path, index_col=0)
    raw_total_counts = df.sum(axis=1)
    adata = ad.AnnData(
        X=df.values.astype(np.float32), obs=pd.DataFrame(index=df.index), var=pd.DataFrame(index=df.columns)
    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata, raw_total_counts


def _compute_layer_score(
    cancer_layer: ad.AnnData, cancer_prop: pd.Series, raw_total_counts: pd.Series, specimen: str
) -> tuple[pd.Series | None, dict]:
    """Score late E2 on SACE cancer-layer spots. Returns (scores, qc_dict)."""
    total_counts = raw_total_counts

    common = cancer_layer.obs_names.intersection(cancer_prop.index)
    prop_aligned = cancer_prop.reindex(common)
    counts_aligned = total_counts.reindex(common)

    mask = (prop_aligned > CANCER_PROP_THRESHOLD) & (counts_aligned >= MIN_CANCER_LAYER_COUNTS)
    n_evaluable = int(mask.sum())

    qc = {
        "n_spots_total": len(common),
        "n_above_prop_threshold": int((prop_aligned > CANCER_PROP_THRESHOLD).sum()),
        "n_above_count_threshold": int((counts_aligned >= MIN_CANCER_LAYER_COUNTS).sum()),
        "n_evaluable": n_evaluable,
        "mean_cancer_prop_evaluable": float(prop_aligned[mask].mean()) if n_evaluable > 0 else 0.0,
        "median_cancer_layer_counts": float(counts_aligned[mask].median()) if n_evaluable > 0 else 0.0,
    }

    late_set = e2.SIGNATURE_GENES["late_up"]
    late_dn_set = e2.SIGNATURE_GENES["late_dn"]
    up_overlap = sum(1 for g in late_set if g in cancer_layer.var_names)
    dn_overlap = sum(1 for g in late_dn_set if g in cancer_layer.var_names)
    qc["overlap"] = {"late_up": up_overlap, "late_dn": dn_overlap}

    print(
        f"  {specimen}: {n_evaluable} evaluable spots "
        f"(prop>{CANCER_PROP_THRESHOLD}: {qc['n_above_prop_threshold']}, "
        f"counts>={MIN_CANCER_LAYER_COUNTS}: {qc['n_above_count_threshold']}) "
        f"late_up={up_overlap} late_dn={dn_overlap}"
    )

    if n_evaluable == 0:
        return None, qc

    evaluable_barcodes = common[mask]
    layer_sub = cancer_layer[evaluable_barcodes].copy()
    scores = e2.compute_directional_score_layer(layer_sub, "late")

    if scores is not None:
        # Reindex to full spot set — non-evaluable spots get NaN
        full_scores = pd.Series(np.nan, index=common, name="late_e2")
        full_scores.loc[scores.index] = scores.values
        return full_scores, qc
    return None, qc


def _load_panel_dims(panel_sizes_path: Path | None) -> tuple[float, float] | None:
    if panel_sizes_path is None or not panel_sizes_path.exists():
        return None
    data = json.loads(panel_sizes_path.read_text())
    a = data.get("A", {})
    w = float(a.get("w_mm", 0))
    h = float(a.get("h_mm", 0))
    if w > 0 and h > 0:
        return w, h
    return None


def _render(loaded, e2_vmax, e2_tag, cl_vmax, out_path: Path, *, w_in: float, h_in: float) -> None:
    """Render the 2-row 4-col panel (cancer prop + late E2)."""
    n_cols = len(loaded)
    n_rows = 2
    fig = plt.figure(figsize=(w_in, h_in))
    gs = fig.add_gridspec(
        n_rows,
        n_cols + 1,
        width_ratios=[1.0] * n_cols + [0.07],
        wspace=0.06,
        hspace=0.22,
        left=0.045,
        right=0.955,
        top=0.92,
        bottom=0.02,
    )

    cl_norm = Normalize(vmin=0, vmax=cl_vmax)
    e2_norm = Normalize(vmin=-e2_vmax, vmax=+e2_vmax)
    cl_mappable = None
    e2_mappable = None

    for col, (label, adata, lib_id) in enumerate(loaded):
        crop = _spot_crop(adata)

        ax0 = fig.add_subplot(gs[0, col])
        sq.pl.spatial_scatter(
            adata,
            color="cancer_prop",
            library_id=lib_id,
            img=True,
            img_res_key="hires",
            crop_coord=crop,
            cmap="Blues",
            norm=cl_norm,
            size=SQ_SPOT_SIZE,
            shape="hex",
            alpha=0.7,
            na_color=NA_RGBA,
            colorbar=False,
            frameon=False,
            ax=ax0,
            return_ax=True,
            title="",
        )
        ax0.set_title(label, fontsize=FONT_TITLE, fontfamily="Arial", pad=4)
        if col == 0:
            ax0.set_ylabel("Cancer\nproportion", fontsize=FONT_ROWLABEL, fontfamily="Arial", labelpad=10)
            cl_mappable = mcm.ScalarMappable(norm=cl_norm, cmap="Blues")
            cl_mappable.set_array([])

        ax1 = fig.add_subplot(gs[1, col])
        sq.pl.spatial_scatter(
            adata,
            color="late_e2",
            library_id=lib_id,
            img=True,
            img_res_key="hires",
            crop_coord=crop,
            cmap="RdBu_r",
            norm=e2_norm,
            size=SQ_SPOT_SIZE,
            shape="hex",
            alpha=0.7,
            na_color=NA_RGBA,
            colorbar=False,
            frameon=False,
            ax=ax1,
            return_ax=True,
            title="",
        )
        if col == 0:
            ax1.set_ylabel("Late E2\n(UP − DN)", fontsize=FONT_ROWLABEL, fontfamily="Arial", labelpad=10)
            e2_mappable = mcm.ScalarMappable(norm=e2_norm, cmap="RdBu_r")
            e2_mappable.set_array([])

    cax_top = fig.add_subplot(gs[0, n_cols])
    cbar_top = fig.colorbar(cl_mappable, cax=cax_top)
    cbar_top.set_label("Proportion", fontsize=FONT_CBAR_LABEL, fontfamily="Arial")
    cbar_top.ax.tick_params(labelsize=FONT_CBAR_TICK)

    cax_e2 = fig.add_subplot(gs[1, n_cols])
    cbar_e2 = fig.colorbar(e2_mappable, cax=cax_e2)
    cbar_e2.set_label("Late E2 (UP − DN)", fontsize=FONT_CBAR_LABEL, fontfamily="Arial")
    cbar_e2.ax.tick_params(labelsize=FONT_CBAR_TICK)

    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.40)
    plt.close(fig)
    print(f"Saved: {out_path}  ({out_path.stat().st_size / 1024:.1f} KB)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--panel-sizes", type=Path, default=None, help="Path to panel_sizes.json (drives raster figsize for panel A)"
    )
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_PANEL_DIR.mkdir(parents=True, exist_ok=True)

    qc: dict = {"per_sample": {}, "global": {}}
    loaded = []
    for m3_basename, specimen, label, he_specimen, unified_specimen in SAMPLES:
        print(f"Loading {m3_basename} ...")
        adata = ad.read_h5ad(MODULE3_DIR / m3_basename)

        cancer_cols = [c for c in adata.obs.columns if "Cancer" in c]
        if cancer_cols:
            adata.obs["cancer_prop"] = adata.obs[cancer_cols].sum(axis=1).astype(float)
        else:
            adata.obs["cancer_prop"] = 0.0

        print(f"  Loading SACE cancer layer for {unified_specimen} ...")
        cancer_layer, raw_total_counts = _load_cancer_layer(unified_specimen)
        cancer_prop = adata.obs["cancer_prop"].reindex(cancer_layer.obs_names)

        late_scores, sample_qc = _compute_layer_score(
            cancer_layer, adata.obs["cancer_prop"], raw_total_counts, unified_specimen
        )
        qc["per_sample"][specimen] = sample_qc

        if late_scores is not None:
            adata.obs["late_e2"] = adata.obs.index.map(lambda bc: late_scores.get(bc, np.nan)).astype(float)
        else:
            adata.obs["late_e2"] = np.nan

        print(f"  Injecting spatial metadata for {he_specimen} ...")
        lib_id = inject_spatial_metadata(adata, he_specimen)
        loaded.append((label, adata, lib_id))

    cl_vmax = float(max(a.obs["cancer_prop"].max() for _, a, _ in loaded))

    e2_vals = np.concatenate([a.obs["late_e2"].values for _, a, _ in loaded])
    e2_finite = np.abs(e2_vals[np.isfinite(e2_vals)])
    e2_vmax, e2_tag = e2.decide_percentile(e2_finite)
    p99 = float(np.nanpercentile(e2_finite, 99)) if e2_finite.size else 1.0
    p95 = float(np.nanpercentile(e2_finite, 95)) if e2_finite.size else 1.0
    qc["global"] = {
        "e2_vmax": e2_vmax,
        "e2_tag": e2_tag,
        "p99": p99,
        "p95": p95,
        "ratio": p99 / max(p95, 1e-9),
        "cl_vmax": cl_vmax,
        "scoring_method": "SACE_cancer_layer_spot_level",
        "cancer_prop_threshold": CANCER_PROP_THRESHOLD,
        "min_cancer_layer_counts": MIN_CANCER_LAYER_COUNTS,
    }
    print(f"[QC] e2_vmax={e2_vmax:.4f} tag={e2_tag} P99={p99:.4f} P95={p95:.4f} ratio={p99/max(p95,1e-9):.3f}")

    # Per-sample r between cancer_prop and late_e2 (for caption)
    for (label, a, _), (_, spec, _, _, _) in zip(loaded, SAMPLES):
        m = np.isfinite(a.obs["late_e2"]) & (a.obs["cancer_prop"] > CANCER_PROP_THRESHOLD)
        n_eval = int(m.sum())
        qc["per_sample"][spec]["n_evaluable_final"] = n_eval
        if n_eval >= 3:
            vals = a.obs.loc[m, "late_e2"].values
            if np.std(vals) > 0:
                qc["per_sample"][spec]["late_e2_mean"] = float(np.mean(vals))
                qc["per_sample"][spec]["late_e2_std"] = float(np.std(vals))
        print(f"[QC] {spec}: n_evaluable_final={n_eval}")

    # Persist QC JSON
    qc_path = QC_DIR / "panel_A_e2_qc.json"
    qc_path.write_text(json.dumps(qc, indent=2))
    print(f"[QC] wrote {qc_path}")

    # Figsize from panel_sizes.json if provided
    dims = _load_panel_dims(args.panel_sizes)
    if dims:
        w_mm, h_mm = dims
        cbar_margin_in = CBAR_MARGIN_IN
        w_in = (w_mm / 25.4) + cbar_margin_in
        h_in = h_mm / 25.4
    else:
        n_cols = len(loaded)
        w_in = 4.0 * n_cols + 1.4
        h_in = 7.2
    print(f"[QC] figsize_in=({w_in:.2f}, {h_in:.2f}) from_panel_sizes={dims is not None}")

    out_path = OUT_DIR / "panel_A_grid.png"
    _render(loaded, e2_vmax, e2_tag, cl_vmax, out_path, w_in=w_in, h_in=h_in)

    # Alternate percentile QC
    alt_vmax = p95 if e2_tag == "p99" else p99
    alt_tag = "p95" if e2_tag == "p99" else "p99"
    _render(loaded, alt_vmax, alt_tag, cl_vmax, QC_PANEL_DIR / "p99_vs_p95.png", w_in=w_in, h_in=h_in)
    print(f"[QC] wrote {QC_PANEL_DIR}/p99_vs_p95.png (alt percentile={alt_tag})")


if __name__ == "__main__":
    main()
