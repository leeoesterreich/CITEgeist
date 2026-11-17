#!/usr/bin/env python -u
"""Generate Fig 5 Panel A (v42): MDK and secretory-score spatial maps of the v2 eligible cancer compartments.

Rows: MDK expression (N_P) and the secretory score (mean of per-gene z-scores of the 11 secretory genes), the two variables
of the bivariate Moran's I. Columns: P1-S1 (biopsy), P1-S2 (surgical), P4-S2 (surgical, D538G); P4-S1 is not estimable
(too few eligible compartments) and is not drawn. Only eligible Cancer_Cells compartments are coloured; other spots are not
drawn. Source: stamped B5 compartments (`_shared/v2_midkine.py`), the same objects the B9 Moran's I used. No keratin and no
estrogen-response score appear here. The v41 panel (cancer proportion + late E2 from output/module3_unified, a pre-v2
pipeline) is not read.

H&E background via squidpy spatial_scatter with crop_coord to the spot extent.
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
import squidpy as sq
from matplotlib.colors import Normalize

# Root of your local copy of the CITEgeist analysis outputs (see README).
PROJECT = Path("/path/to/CITEgeist_analysis")
OUT_DIR = Path(__file__).resolve().parent / "outputs"
OUT_PNG = OUT_DIR / "panel_A_grid_v42.png"

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _shared import sidecar  # noqa: E402
from _shared.spatial_utils import inject_spatial_metadata  # noqa: E402
from _shared.v2_midkine import load_mdk_secretory  # noqa: E402

SAMPLES: list[tuple[str, str]] = [
    # (specimen_id, display_label)
    ("HCC22-088-P1-S1", "P1-S1\n(biopsy)"),
    ("HCC22-088-P1-S2", "P1-S2\n(surgical)"),
    ("HCC22-088-P4-S2_1i_rep", "P4-S2\n(surgical, D538G)"),
]

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


def _load_panel_dims(panel_sizes_path: Path | None) -> tuple[float, float] | None:
    if panel_sizes_path is None or not panel_sizes_path.exists():
        return None
    a = json.loads(panel_sizes_path.read_text()).get("A", {})
    w, h = float(a.get("w_mm", 0)), float(a.get("h_mm", 0))
    return (w, h) if w > 0 and h > 0 else None


def _load_specimen(specimen: str):
    """AnnData of the eligible cancer spots with obs mdk / secretory and the H&E metadata injected; plus (frame, h5ad path)."""
    df, path = load_mdk_secretory(specimen)
    adata = ad.AnnData(X=np.zeros((len(df), 1), dtype=np.float32))
    adata.obs_names = df.index
    adata.obs["mdk"] = df["mdk"].to_numpy(float)
    adata.obs["secretory"] = df["secretory"].to_numpy(float)
    adata.obsm["spatial"] = df[["x", "y"]].to_numpy(float)
    lib_id = inject_spatial_metadata(adata, specimen)
    return adata, lib_id, path


def _scatter(adata, lib_id, color, cmap, norm, ax):
    sq.pl.spatial_scatter(
        adata,
        color=color,
        library_id=lib_id,
        img=True,
        img_res_key="hires",
        crop_coord=_spot_crop(adata),
        cmap=cmap,
        norm=norm,
        size=SQ_SPOT_SIZE,
        shape="hex",
        alpha=0.7,
        na_color=NA_RGBA,
        colorbar=False,
        frameon=False,
        ax=ax,
        return_ax=True,
        title="",
    )


def _render(loaded, mdk_vmax, sec_lim, out_path: Path, *, w_in: float, h_in: float) -> None:
    n_cols = len(loaded)
    fig = plt.figure(figsize=(w_in, h_in))
    gs = fig.add_gridspec(
        2,
        n_cols + 1,
        width_ratios=[1.0] * n_cols + [0.07],
        wspace=0.04,
        hspace=0.22,
        left=0.10,
        right=0.955,
        top=0.92,
        bottom=0.02,
    )
    mdk_norm = Normalize(vmin=0, vmax=mdk_vmax)
    sec_norm = Normalize(vmin=-sec_lim, vmax=sec_lim)
    mdk_cmap, sec_cmap = "YlOrRd", "PuBu"
    for col, (label, adata, lib_id) in enumerate(loaded):
        ax0 = fig.add_subplot(gs[0, col])
        _scatter(adata, lib_id, "mdk", mdk_cmap, mdk_norm, ax0)
        ax0.set_title(label, fontsize=FONT_TITLE, fontfamily="Arial", pad=4)
        if col == 0:
            ax0.text(
                -0.03,
                0.5,
                "MDK\nexpression",
                transform=ax0.transAxes,
                rotation=90,
                ha="right",
                va="center",
                fontsize=FONT_ROWLABEL,
                fontfamily="Arial",
            )
        ax1 = fig.add_subplot(gs[1, col])
        _scatter(adata, lib_id, "secretory", sec_cmap, sec_norm, ax1)
        if col == 0:
            ax1.text(
                -0.03,
                0.5,
                "Secretory\nscore",
                transform=ax1.transAxes,
                rotation=90,
                ha="right",
                va="center",
                fontsize=FONT_ROWLABEL,
                fontfamily="Arial",
            )
    for row, (norm, cmap, text) in enumerate(((mdk_norm, mdk_cmap, "MDK"), (sec_norm, sec_cmap, "Secretory z-score"))):
        m = mcm.ScalarMappable(norm=norm, cmap=cmap)
        m.set_array([])
        cax = fig.add_subplot(gs[row, n_cols])
        cb = fig.colorbar(m, cax=cax)
        cb.set_label(text, fontsize=FONT_CBAR_LABEL, fontfamily="Arial")
        cb.ax.tick_params(labelsize=FONT_CBAR_TICK)
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

    with sidecar.ReadTracer() as tracer:
        loaded = []
        for specimen, label in SAMPLES:
            print(f"Loading {specimen} (B5 compartments) ...")
            adata, lib_id, _ = _load_specimen(specimen)
            print(f"  {adata.n_obs} eligible cancer spots; MDK>0 in {int((adata.obs['mdk'] > 0).sum())}")
            loaded.append((label, adata, lib_id))

        # shared colour ranges across the three specimens
        mdk_all = np.concatenate([a.obs["mdk"].to_numpy() for _, a, _ in loaded])
        sec_all = np.concatenate([a.obs["secretory"].to_numpy() for _, a, _ in loaded])
        mdk_vmax = float(np.percentile(mdk_all[mdk_all > 0], 99))
        sec_lim = float(np.percentile(np.abs(sec_all), 99))
        print(f"[QC] mdk_vmax(p99 of >0)={mdk_vmax:.3f} secretory_lim(p99 |z|)={sec_lim:.3f}")

        dims = _load_panel_dims(args.panel_sizes)
        if dims:
            w_in, h_in = dims[0] / 25.4 + CBAR_MARGIN_IN, dims[1] / 25.4
        else:
            w_in, h_in = 4.0 * len(loaded) + 1.4, 7.2
        print(f"[QC] figsize_in=({w_in:.2f}, {h_in:.2f}) from_panel_sizes={dims is not None}")
        _render(loaded, mdk_vmax, sec_lim, OUT_PNG, w_in=w_in, h_in=h_in)

    sidecar.stamp_outputs([OUT_PNG], stages=("B5", "B9"), inputs=tracer.inputs(), generator=Path(__file__))
    print(f"[sidecar] {sidecar.sidecar_path(OUT_PNG).name}")


if __name__ == "__main__":
    main()
