#!/usr/bin/env python
"""
analyze_hw120_if.py — Quantitative immunofluorescence analysis of HW120 MDK IF images.

Analyzes MDK protein in MCF7/T47D cell lines (WT vs D538G ESR1 mutant).
Biological question: D538G causes opposite effects on secretory machinery
in MCF7 vs T47D. Quantifies MDK intensity and subcellular localization
(perinuclear enrichment).

Pipeline:
    1. Image discovery and metadata parsing
    2. Image loading (with z-stack max-projection)
    3. Rolling-ball background subtraction
    4. Cellpose nuclear segmentation
    5. Compartment definition (nuclear / perinuclear / cytoplasmic)
    6. Per-cell intensity measurement
    7. Statistical analysis (two-way ANOVA + Tukey HSD)
    8. Figure generation

Usage:
    python analyze_hw120_if.py --data-root /path/to/data --output-dir if_analysis
    python analyze_hw120_if.py --dry-run  # just discover images, no analysis
"""

import matplotlib
matplotlib.use('Agg')

import argparse
import logging
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import tifffile
from scipy import ndimage
from skimage.restoration import rolling_ball
from skimage.segmentation import expand_labels
from skimage.measure import regionprops

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────

DATA_ROOT_DEFAULT = None  # required: the folder of raw z-stack TIFFs (one sub-folder per condition), see wetlab/README.md
OUTPUT_DIR_DEFAULT = str(Path(__file__).resolve().parent / "if_analysis")

# Cell lines we want to analyze (skip controls like 293T, HCC1187)
VALID_CELL_LINES = {"MCF7", "T47D"}

# Folder name patterns to explicitly skip
SKIP_PATTERNS = [
    r"293T",
    r"HCC1187",
    r"No Primary Control",
    r"antibody test",
    r"Antibody Testing",
]

CHANNEL_MAP = {
    "_Hoechst342.tif": "hoechst",
    "_AF555.tif": "ecadherin",   # E-cadherin
    "_AF647.tif": "mdk",         # MDK
}

ROLLING_BALL_RADIUS = 50
CELLPOSE_DIAMETER_10X = 30
CELLPOSE_DIAMETER_60X = 120
PERINUCLEAR_DILATION_PX = 15
MIN_NUCLEUS_AREA = 50  # minimum nucleus area in pixels

FIGURE_DPI = 300
FONT_SIZE_MIN = 10


# ──────────────────────────────────────────────────────────────────────
# 1. Image Discovery
# ──────────────────────────────────────────────────────────────────────

def parse_condition_folder(folder_name: str) -> Optional[Dict[str, str]]:
    """
    Parse experimental condition from folder name.

    Returns dict with keys: cell_line, genotype, condition
    or None if the folder should be skipped.

    Examples:
        'MCF7 WT 5 Day' -> {cell_line: MCF7, genotype: WT, condition: 5_Day}
        'T47D D538G 5 day +E2' -> {cell_line: T47D, genotype: D538G, condition: 5_Day_E2}
        'MCF7 WT Full MDK rb 1-100 ...' -> {cell_line: MCF7, genotype: WT, condition: Full}
    """
    # Skip known control / non-experimental folders
    for pat in SKIP_PATTERNS:
        if re.search(pat, folder_name, re.IGNORECASE):
            logger.debug(f"Skipping control folder: {folder_name}")
            return None

    # Try to match cell_line and genotype at the start
    m = re.match(
        r'^(MCF7|T47D)\s+(WT|D538G)\s+(.+)$',
        folder_name,
        re.IGNORECASE,
    )
    if not m:
        logger.debug(f"Could not parse folder: {folder_name}")
        return None

    cell_line = m.group(1).upper()
    genotype = m.group(2).upper()
    rest = m.group(3).strip()

    if cell_line not in VALID_CELL_LINES:
        return None

    # Parse the condition from the rest of the folder name
    # "Full MDK rb 1-100 647 ecad ms 1-100 555" -> "Full"
    # "5 Day" -> "5_Day"
    # "5 day +E2" -> "5_Day_E2"
    # "1 Day" -> "1_Day"
    # "3 Day" -> "3_Day"

    # Check for "Full" first
    if rest.lower().startswith("full"):
        condition = "Full"
    else:
        # Try to match N Day [+E2]
        day_match = re.match(r'^(\d+)\s+[Dd]ay\s*(\+E2)?', rest)
        if day_match:
            n_days = day_match.group(1)
            e2 = day_match.group(2)
            condition = f"{n_days}_Day"
            if e2:
                condition += "_E2"
        else:
            # Fallback: use rest as-is, cleaned
            condition = re.sub(r'\s+', '_', rest)

    return {
        "cell_line": cell_line,
        "genotype": genotype,
        "condition": condition,
    }


def discover_fovs(folder_path: Path) -> List[Dict]:
    """
    Discover FOVs within a condition folder.

    Returns list of dicts with keys:
        fov_name, magnification, is_zstack, channels: {channel_name: filepath}
    """
    fovs = defaultdict(lambda: {"channels": {}})

    for fname in sorted(os.listdir(folder_path)):
        fpath = folder_path / fname

        # Skip directories, .nd2, _RGB files
        if fpath.is_dir():
            continue
        if not fname.endswith(".tif"):
            continue
        if "_RGB" in fname:
            continue

        # Identify channel
        channel_key = None
        base_name = fname
        for suffix, ch_name in CHANNEL_MAP.items():
            if fname.endswith(suffix):
                channel_key = ch_name
                base_name = fname[: -len(suffix)]
                break

        if channel_key is None:
            # This is the composite TIF (e.g., "10x 1.tif") — skip
            continue

        # Parse magnification and FOV id
        # "10x 1" or "60x z stack 2"
        is_zstack = "z stack" in base_name.lower() or "zstack" in base_name.lower()
        if "60x" in base_name.lower():
            mag = "60x"
        elif "10x" in base_name.lower():
            mag = "10x"
        else:
            mag = "unknown"

        fov_key = base_name.strip()
        fovs[fov_key]["fov_name"] = fov_key
        fovs[fov_key]["magnification"] = mag
        fovs[fov_key]["is_zstack"] = is_zstack
        fovs[fov_key]["channels"][channel_key] = str(fpath)

    # Filter to FOVs that at least have hoechst + mdk
    results = []
    for fov_key, fov_data in fovs.items():
        chs = fov_data["channels"]
        if "hoechst" in chs and "mdk" in chs:
            results.append(fov_data)
        else:
            logger.debug(
                f"Skipping FOV {fov_key} — missing channels "
                f"(have: {list(chs.keys())})"
            )
    return results


def discover_images(data_root: str) -> pd.DataFrame:
    """
    Walk the data directory, parse metadata, and return a DataFrame of all FOVs.

    Columns: cell_line, genotype, condition, fov_name, magnification, is_zstack,
             hoechst_path, mdk_path, ecadherin_path, condition_folder
    """
    root = Path(data_root)
    if not root.is_dir():
        logger.error(f"Data root does not exist: {data_root}")
        sys.exit(1)

    records = []
    for entry in sorted(os.listdir(root)):
        folder_path = root / entry
        if not folder_path.is_dir():
            continue

        meta = parse_condition_folder(entry)
        if meta is None:
            continue

        fovs = discover_fovs(folder_path)
        for fov in fovs:
            records.append({
                "cell_line": meta["cell_line"],
                "genotype": meta["genotype"],
                "condition": meta["condition"],
                "condition_folder": entry,
                "fov_name": fov["fov_name"],
                "magnification": fov["magnification"],
                "is_zstack": fov["is_zstack"],
                "hoechst_path": fov["channels"].get("hoechst"),
                "mdk_path": fov["channels"].get("mdk"),
                "ecadherin_path": fov["channels"].get("ecadherin"),
            })

    df = pd.DataFrame(records)
    logger.info(f"Discovered {len(df)} FOVs across {df['condition_folder'].nunique()} conditions")
    if len(df) > 0:
        for mag in df["magnification"].unique():
            n = (df["magnification"] == mag).sum()
            logger.info(f"  {mag}: {n} FOVs")
    return df


# ──────────────────────────────────────────────────────────────────────
# 2. Image Loading
# ──────────────────────────────────────────────────────────────────────

def load_image(path: str) -> np.ndarray:
    """
    Load a TIFF image, max-project z-stacks, convert to float32.

    Divides by dtype max (255 for uint8, 65535 for uint16) to normalize
    to [0, 1] range without per-FOV min-max normalization.
    """
    img = tifffile.imread(path)

    # Max-project z-stacks (3D arrays with shape Z, H, W)
    if img.ndim == 3:
        logger.debug(f"  Z-stack detected ({img.shape[0]} slices), max-projecting")
        img = img.max(axis=0)

    # Handle multi-channel composites that might slip through
    if img.ndim == 3:
        # Could be H, W, C — take first channel
        logger.warning(f"  Unexpected 3D shape {img.shape} after projection, taking channel 0")
        img = img[:, :, 0]

    # Convert to float32
    if img.dtype == np.uint8:
        img = img.astype(np.float32) / 255.0
    elif img.dtype == np.uint16:
        img = img.astype(np.float32) / 65535.0
    else:
        img = img.astype(np.float32)

    return img


# ──────────────────────────────────────────────────────────────────────
# 3. Background Subtraction
# ──────────────────────────────────────────────────────────────────────

def subtract_background(img: np.ndarray, radius: int = ROLLING_BALL_RADIUS) -> np.ndarray:
    """
    Rolling-ball background subtraction.

    skimage.restoration.rolling_ball returns the estimated background.
    """
    bg = rolling_ball(img, radius=radius)
    corrected = img - bg
    corrected = np.clip(corrected, 0, None)
    return corrected


# ──────────────────────────────────────────────────────────────────────
# 4. Nuclear Segmentation
# ──────────────────────────────────────────────────────────────────────

def segment_nuclei(
    hoechst: np.ndarray,
    magnification: str,
    gpu: bool = True,
) -> np.ndarray:
    """
    Segment nuclei using Cellpose nuclei model.

    Args:
        hoechst: 2D float32 Hoechst image (0-1 range).
        magnification: '10x' or '60x' to select diameter.
        gpu: Whether to use GPU.

    Returns:
        Integer label mask (0 = background).
    """
    from cellpose import models

    diameter = CELLPOSE_DIAMETER_10X if magnification == "10x" else CELLPOSE_DIAMETER_60X

    model = models.CellposeModel(model_type='nuclei', gpu=gpu)

    # Cellpose expects (H, W) for single-channel
    # Scale to uint16-like range for better segmentation
    img_for_cp = (hoechst * 65535).astype(np.float32)

    masks, flows, styles = model.eval(
        img_for_cp,
        diameter=diameter,
        channels=[0, 0],  # grayscale
        flow_threshold=0.4,
        cellprob_threshold=0.0,
    )

    # Filter small objects
    props = regionprops(masks)
    for p in props:
        if p.area < MIN_NUCLEUS_AREA:
            masks[masks == p.label] = 0

    # Relabel
    unique_labels = np.unique(masks)
    unique_labels = unique_labels[unique_labels > 0]
    new_masks = np.zeros_like(masks)
    for i, lbl in enumerate(unique_labels, start=1):
        new_masks[masks == lbl] = i

    n_nuclei = len(unique_labels)
    logger.info(f"    Segmented {n_nuclei} nuclei (diameter={diameter})")
    return new_masks


def save_segmentation_qc(
    hoechst: np.ndarray,
    masks: np.ndarray,
    output_path: str,
) -> None:
    """Save a QC overlay: Hoechst in gray with nuclear contours in red."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax.imshow(hoechst, cmap='gray', vmin=0, vmax=np.percentile(hoechst, 99.5))

    # Draw contours
    from skimage.segmentation import find_boundaries
    boundaries = find_boundaries(masks, mode='outer')
    overlay = np.zeros((*hoechst.shape, 4), dtype=np.float32)
    overlay[boundaries, 0] = 1.0  # red
    overlay[boundaries, 3] = 0.8  # alpha
    ax.imshow(overlay)

    ax.set_title(f"Nuclear Segmentation QC (n={masks.max()})", fontsize=12)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


# ──────────────────────────────────────────────────────────────────────
# 5. Compartment Definition
# ──────────────────────────────────────────────────────────────────────

def build_compartments(
    nuclear_mask: np.ndarray,
    magnification: str,
) -> Dict[str, np.ndarray]:
    """
    Build compartment masks from nuclear segmentation.

    10x: nuclear + whole-cell (Voronoi)
    60x: nuclear + perinuclear ring + cytoplasmic

    Returns dict of boolean masks keyed by compartment name.
    Each mask has shape (n_cells, H, W) — but we return label masks instead
    for efficiency: {compartment_name: label_array} where label IDs match
    nuclear_mask IDs.
    """
    h, w = nuclear_mask.shape
    max_dim = max(h, w)

    # Voronoi territories via expand_labels
    voronoi = expand_labels(nuclear_mask, distance=max_dim)

    compartments = {
        "nuclear": nuclear_mask,
        "whole_cell": voronoi,
    }

    if magnification == "60x":
        # Perinuclear ring: dilate each nucleus, subtract nucleus, clip to Voronoi
        struct = np.ones((3, 3), dtype=bool)

        # Dilate entire nuclear mask, then assign cell IDs
        # We need per-cell dilation, so iterate labels
        perinuclear = np.zeros_like(nuclear_mask)
        cytoplasmic = np.zeros_like(nuclear_mask)

        unique_ids = np.unique(nuclear_mask)
        unique_ids = unique_ids[unique_ids > 0]

        for cell_id in unique_ids:
            nuc_binary = nuclear_mask == cell_id
            vor_binary = voronoi == cell_id

            # Dilate nucleus
            dilated = ndimage.binary_dilation(
                nuc_binary, structure=struct, iterations=PERINUCLEAR_DILATION_PX
            )

            # Perinuclear = dilated minus nucleus, clipped to Voronoi territory
            peri = dilated & ~nuc_binary & vor_binary
            perinuclear[peri] = cell_id

            # Cytoplasmic = Voronoi minus nucleus minus perinuclear
            cyto = vor_binary & ~nuc_binary & ~peri
            cytoplasmic[cyto] = cell_id

        compartments["perinuclear"] = perinuclear
        compartments["cytoplasmic"] = cytoplasmic

    return compartments


# ──────────────────────────────────────────────────────────────────────
# 6. Per-Cell Measurement
# ──────────────────────────────────────────────────────────────────────

def measure_cells(
    mdk: np.ndarray,
    ecadherin: Optional[np.ndarray],
    compartments: Dict[str, np.ndarray],
    magnification: str,
) -> pd.DataFrame:
    """
    Measure per-cell intensities in each compartment.

    Returns DataFrame with one row per cell.
    """
    nuclear_mask = compartments["nuclear"]
    cell_ids = np.unique(nuclear_mask)
    cell_ids = cell_ids[cell_ids > 0]

    records = []
    for cell_id in cell_ids:
        row = {"cell_id": int(cell_id)}

        # Nuclear MDK
        nuc_mask = nuclear_mask == cell_id
        nuc_pixels = mdk[nuc_mask]
        if len(nuc_pixels) == 0:
            continue
        row["mdk_nuclear_mean"] = float(nuc_pixels.mean())
        row["mdk_nuclear_median"] = float(np.median(nuc_pixels))
        row["nuclear_area"] = int(nuc_mask.sum())

        # Whole-cell MDK
        wc_mask = compartments["whole_cell"] == cell_id
        wc_pixels = mdk[wc_mask]
        row["mdk_wholecell_mean"] = float(wc_pixels.mean()) if len(wc_pixels) > 0 else np.nan
        row["wholecell_area"] = int(wc_mask.sum())

        if magnification == "60x":
            # Perinuclear MDK
            peri_mask = compartments.get("perinuclear")
            if peri_mask is not None:
                pm = peri_mask == cell_id
                peri_pixels = mdk[pm]
                row["mdk_perinuclear_mean"] = float(peri_pixels.mean()) if len(peri_pixels) > 0 else np.nan
                row["perinuclear_area"] = int(pm.sum())
            else:
                row["mdk_perinuclear_mean"] = np.nan
                row["perinuclear_area"] = 0

            # Cytoplasmic MDK
            cyto_mask = compartments.get("cytoplasmic")
            if cyto_mask is not None:
                cm = cyto_mask == cell_id
                cyto_pixels = mdk[cm]
                row["mdk_cytoplasmic_mean"] = float(cyto_pixels.mean()) if len(cyto_pixels) > 0 else np.nan
                row["cytoplasmic_area"] = int(cm.sum())
            else:
                row["mdk_cytoplasmic_mean"] = np.nan
                row["cytoplasmic_area"] = 0

            # Perinuclear enrichment ratios
            if row.get("mdk_perinuclear_mean") and row.get("mdk_nuclear_mean"):
                nuc_val = row["mdk_nuclear_mean"]
                peri_val = row["mdk_perinuclear_mean"]
                cyto_val = row.get("mdk_cytoplasmic_mean", np.nan)

                row["peri_nuc_ratio"] = peri_val / nuc_val if nuc_val > 0 else np.nan
                row["peri_cyto_ratio"] = peri_val / cyto_val if cyto_val and cyto_val > 0 else np.nan

        # E-cadherin
        if ecadherin is not None:
            row["ecad_nuclear_mean"] = float(ecadherin[nuc_mask].mean())
            row["ecad_wholecell_mean"] = float(ecadherin[wc_mask].mean()) if len(wc_pixels) > 0 else np.nan

        records.append(row)

    return pd.DataFrame(records)


# ──────────────────────────────────────────────────────────────────────
# 7. Statistics
# ──────────────────────────────────────────────────────────────────────

def run_statistics(
    df_10x: pd.DataFrame,
    df_60x: pd.DataFrame,
    output_dir: str,
) -> None:
    """
    Aggregate to FOV-level means, run two-way ANOVA and Tukey HSD.

    Saves statistical_summary.txt.
    """
    import statsmodels.api as sm
    from statsmodels.formula.api import ols
    from statsmodels.stats.multicomp import pairwise_tukeyhsd

    out_path = os.path.join(output_dir, "statistical_summary.txt")

    with open(out_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("HW120 MDK IF — Statistical Summary\n")
        f.write("=" * 70 + "\n\n")

        # ── 10x Analysis ──
        if len(df_10x) > 0:
            f.write("─" * 50 + "\n")
            f.write("10x ANALYSIS: MDK Nuclear Intensity\n")
            f.write("─" * 50 + "\n\n")

            # Aggregate to FOV-level means
            fov_means = (
                df_10x.groupby(["cell_line", "genotype", "condition", "fov_name"])
                ["mdk_nuclear_mean"]
                .mean()
                .reset_index()
            )
            fov_means.rename(columns={"mdk_nuclear_mean": "fov_mean_mdk"}, inplace=True)

            f.write(f"N FOVs: {len(fov_means)}\n")
            f.write(f"Groups: {fov_means.groupby(['cell_line', 'genotype', 'condition']).size().to_string()}\n\n")

            # Two-way ANOVA per condition
            for cond in sorted(fov_means["condition"].unique()):
                subset = fov_means[fov_means["condition"] == cond]
                if len(subset) < 4:
                    f.write(f"\nCondition {cond}: too few FOVs ({len(subset)}) for ANOVA\n")
                    continue

                f.write(f"\n--- Condition: {cond} ---\n")
                try:
                    model = ols(
                        'fov_mean_mdk ~ C(cell_line) * C(genotype)', data=subset
                    ).fit()
                    anova_table = sm.stats.anova_lm(model, typ=2)
                    f.write(f"\nTwo-way ANOVA:\n{anova_table.to_string()}\n")

                    # Tukey HSD on combined group label
                    subset = subset.copy()
                    subset["group"] = subset["cell_line"] + "_" + subset["genotype"]
                    if subset["group"].nunique() > 1:
                        tukey = pairwise_tukeyhsd(
                            subset["fov_mean_mdk"], subset["group"], alpha=0.05
                        )
                        f.write(f"\nTukey HSD:\n{tukey.summary()}\n")
                except Exception as e:
                    f.write(f"\nANOVA failed: {e}\n")

        # ── 60x Analysis ──
        if len(df_60x) > 0:
            f.write("\n\n")
            f.write("─" * 50 + "\n")
            f.write("60x ANALYSIS: Perinuclear Enrichment\n")
            f.write("─" * 50 + "\n\n")

            # Aggregate to FOV-level means
            ratio_cols = [c for c in ["peri_nuc_ratio", "peri_cyto_ratio"] if c in df_60x.columns]
            if ratio_cols:
                fov_60x = (
                    df_60x.groupby(["cell_line", "genotype", "condition", "fov_name"])
                    [ratio_cols + ["mdk_nuclear_mean", "mdk_perinuclear_mean"]]
                    .mean()
                    .reset_index()
                )

                f.write(f"N FOVs (60x): {len(fov_60x)}\n\n")

                for metric in ratio_cols:
                    f.write(f"\n=== {metric} ===\n")
                    for cond in sorted(fov_60x["condition"].unique()):
                        subset = fov_60x[fov_60x["condition"] == cond].dropna(subset=[metric])
                        if len(subset) < 4:
                            f.write(f"\nCondition {cond}: too few FOVs ({len(subset)})\n")
                            continue

                        f.write(f"\n--- Condition: {cond} ---\n")
                        try:
                            model = ols(
                                f'{metric} ~ C(cell_line) * C(genotype)', data=subset
                            ).fit()
                            anova_table = sm.stats.anova_lm(model, typ=2)
                            f.write(f"\nTwo-way ANOVA:\n{anova_table.to_string()}\n")

                            subset = subset.copy()
                            subset["group"] = subset["cell_line"] + "_" + subset["genotype"]
                            if subset["group"].nunique() > 1:
                                tukey = pairwise_tukeyhsd(
                                    subset[metric], subset["group"], alpha=0.05
                                )
                                f.write(f"\nTukey HSD:\n{tukey.summary()}\n")
                        except Exception as e:
                            f.write(f"\nANOVA failed: {e}\n")

    logger.info(f"Statistical summary saved to {out_path}")


# ──────────────────────────────────────────────────────────────────────
# 8. Figures
# ──────────────────────────────────────────────────────────────────────

def _setup_figure_style():
    """Set consistent matplotlib/seaborn style."""
    plt.rcParams.update({
        'font.size': FONT_SIZE_MIN,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.dpi': FIGURE_DPI,
    })
    sns.set_style("whitegrid")


def _aggregate_fov(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    """Aggregate per-cell measurements to FOV-level means."""
    return (
        df.groupby(["cell_line", "genotype", "condition", "fov_name"])
        [value_col]
        .mean()
        .reset_index()
    )


def plot_mdk_intensity_barplot(df: pd.DataFrame, output_dir: str) -> None:
    """Grouped bar by condition, hue=cell_line+genotype."""
    _setup_figure_style()

    fov_df = _aggregate_fov(df, "mdk_nuclear_mean")
    fov_df["group"] = fov_df["cell_line"] + " " + fov_df["genotype"]

    # Sort conditions logically
    cond_order = []
    for c in ["1_Day", "3_Day", "5_Day", "Full", "1_Day_E2", "3_Day_E2", "5_Day_E2"]:
        if c in fov_df["condition"].values:
            cond_order.append(c)
    # Add any remaining
    for c in sorted(fov_df["condition"].unique()):
        if c not in cond_order:
            cond_order.append(c)

    group_order = ["MCF7 WT", "MCF7 D538G", "T47D WT", "T47D D538G"]
    group_order = [g for g in group_order if g in fov_df["group"].values]

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.barplot(
        data=fov_df, x="condition", y="mdk_nuclear_mean", hue="group",
        order=cond_order, hue_order=group_order,
        errorbar="se", capsize=0.05, ax=ax,
        palette="Set2",
    )
    ax.set_xlabel("Condition")
    ax.set_ylabel("MDK Nuclear Intensity (mean)")
    ax.set_title("MDK Intensity by Condition and Cell Line/Genotype")
    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()

    path = os.path.join(output_dir, "mdk_intensity_barplot.png")
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved {path}")


def plot_mdk_timecourse(df: pd.DataFrame, output_dir: str) -> None:
    """Line plot of MDK intensity at 1/3/5 day timepoints, 4 groups."""
    _setup_figure_style()

    # Filter to Day conditions (no E2, no Full)
    day_df = df[df["condition"].isin(["1_Day", "3_Day", "5_Day"])].copy()
    if len(day_df) == 0:
        logger.warning("No Day-only conditions found for timecourse plot")
        return

    fov_df = _aggregate_fov(day_df, "mdk_nuclear_mean")
    fov_df["group"] = fov_df["cell_line"] + " " + fov_df["genotype"]
    fov_df["day"] = fov_df["condition"].str.extract(r'(\d+)').astype(int)

    fig, ax = plt.subplots(figsize=(7, 5))

    group_order = ["MCF7 WT", "MCF7 D538G", "T47D WT", "T47D D538G"]
    palette = dict(zip(group_order, sns.color_palette("Set2", 4)))
    markers = {"MCF7 WT": "o", "MCF7 D538G": "s", "T47D WT": "^", "T47D D538G": "D"}

    for group in group_order:
        gdata = fov_df[fov_df["group"] == group]
        if len(gdata) == 0:
            continue
        # Compute mean and SEM per day
        summary = gdata.groupby("day")["mdk_nuclear_mean"].agg(["mean", "sem"]).reset_index()
        ax.errorbar(
            summary["day"], summary["mean"], yerr=summary["sem"],
            marker=markers.get(group, "o"), label=group,
            color=palette.get(group), linewidth=2, capsize=4,
        )

    ax.set_xlabel("Days (E2 deprivation)")
    ax.set_ylabel("MDK Nuclear Intensity (mean)")
    ax.set_title("MDK Timecourse: E2 Deprivation")
    ax.set_xticks([1, 3, 5])
    ax.legend(title="Group")
    fig.tight_layout()

    path = os.path.join(output_dir, "mdk_timecourse.png")
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved {path}")


def plot_perinuclear_enrichment(df_60x: pd.DataFrame, output_dir: str) -> None:
    """Bar plot of perinuclear/cytoplasmic ratio (60x data)."""
    _setup_figure_style()

    if "peri_cyto_ratio" not in df_60x.columns:
        logger.warning("No peri_cyto_ratio in 60x data")
        return

    fov_df = _aggregate_fov(df_60x.dropna(subset=["peri_cyto_ratio"]), "peri_cyto_ratio")
    if len(fov_df) == 0:
        logger.warning("No valid 60x FOV means for perinuclear enrichment plot")
        return

    fov_df["group"] = fov_df["cell_line"] + " " + fov_df["genotype"]
    group_order = ["MCF7 WT", "MCF7 D538G", "T47D WT", "T47D D538G"]
    group_order = [g for g in group_order if g in fov_df["group"].values]

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(
        data=fov_df, x="condition", y="peri_cyto_ratio", hue="group",
        hue_order=group_order, errorbar="se", capsize=0.05,
        palette="Set2", ax=ax,
    )
    ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.8, label='No enrichment')
    ax.set_xlabel("Condition")
    ax.set_ylabel("Perinuclear / Cytoplasmic Ratio")
    ax.set_title("MDK Perinuclear Enrichment (60x)")
    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()

    path = os.path.join(output_dir, "perinuclear_enrichment.png")
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved {path}")


def plot_interaction(df: pd.DataFrame, output_dir: str) -> None:
    """Interaction plot: cell_line x genotype at 5-day condition."""
    _setup_figure_style()

    day5 = df[df["condition"] == "5_Day"]
    if len(day5) == 0:
        logger.warning("No 5_Day condition data for interaction plot")
        return

    fov_df = _aggregate_fov(day5, "mdk_nuclear_mean")

    fig, ax = plt.subplots(figsize=(6, 5))

    # Plot means with error bars
    summary = fov_df.groupby(["cell_line", "genotype"])["mdk_nuclear_mean"].agg(
        ["mean", "sem"]
    ).reset_index()

    for cl in ["MCF7", "T47D"]:
        sub = summary[summary["cell_line"] == cl].sort_values("genotype")
        if len(sub) == 0:
            continue
        ax.errorbar(
            sub["genotype"], sub["mean"], yerr=sub["sem"],
            marker='o', linewidth=2, capsize=5, markersize=8,
            label=cl,
        )

    # Overlay individual FOV data points
    for cl, marker in [("MCF7", "o"), ("T47D", "^")]:
        sub = fov_df[fov_df["cell_line"] == cl]
        if len(sub) == 0:
            continue
        jitter = np.random.normal(0, 0.03, len(sub))
        x_pos = sub["genotype"].map({"D538G": 0, "WT": 1}).values + jitter
        ax.scatter(x_pos, sub["mdk_nuclear_mean"], alpha=0.4, s=30, marker=marker)

    ax.set_xlabel("Genotype")
    ax.set_ylabel("MDK Nuclear Intensity (mean)")
    ax.set_title("Cell Line x Genotype Interaction (5-Day E2 Deprivation)")
    ax.legend(title="Cell Line")
    fig.tight_layout()

    path = os.path.join(output_dir, "interaction_plot.png")
    fig.savefig(path, dpi=FIGURE_DPI, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved {path}")


# ──────────────────────────────────────────────────────────────────────
# 9. Main Pipeline
# ──────────────────────────────────────────────────────────────────────

def process_fov(
    row: pd.Series,
    output_dir: str,
    gpu: bool = True,
) -> Optional[pd.DataFrame]:
    """
    Process a single FOV: load, background subtract, segment, measure.

    Returns per-cell measurement DataFrame with metadata columns added.
    """
    fov_name = row["fov_name"]
    mag = row["magnification"]
    logger.info(f"  Processing {fov_name} ({mag})...")

    # Load channels
    hoechst = load_image(row["hoechst_path"])
    mdk_raw = load_image(row["mdk_path"])
    ecadherin = None
    if pd.notna(row.get("ecadherin_path")):
        ecadherin = load_image(row["ecadherin_path"])

    # Background subtraction on signal channels
    logger.debug("    Background subtraction...")
    mdk = subtract_background(mdk_raw)
    if ecadherin is not None:
        ecadherin = subtract_background(ecadherin)

    # Nuclear segmentation
    logger.debug("    Nuclear segmentation...")
    masks = segment_nuclei(hoechst, mag, gpu=gpu)

    if masks.max() == 0:
        logger.warning(f"    No nuclei detected in {fov_name}")
        return None

    # Save QC overlay
    qc_dir = os.path.join(output_dir, "qc_overlays")
    os.makedirs(qc_dir, exist_ok=True)
    safe_name = re.sub(r'[^\w\-.]', '_', f"{row['condition_folder']}_{fov_name}")
    save_segmentation_qc(hoechst, masks, os.path.join(qc_dir, f"{safe_name}_qc.png"))

    # Build compartments
    compartments = build_compartments(masks, mag)

    # Measure per-cell intensities
    cell_df = measure_cells(mdk, ecadherin, compartments, mag)

    if len(cell_df) == 0:
        logger.warning(f"    No valid measurements in {fov_name}")
        return None

    # Add metadata
    cell_df["cell_line"] = row["cell_line"]
    cell_df["genotype"] = row["genotype"]
    cell_df["condition"] = row["condition"]
    cell_df["condition_folder"] = row["condition_folder"]
    cell_df["fov_name"] = fov_name
    cell_df["magnification"] = mag

    return cell_df


def main():
    parser = argparse.ArgumentParser(
        description="HW120 MDK IF quantitative analysis pipeline"
    )
    parser.add_argument(
        "--data-root", type=str, default=DATA_ROOT_DEFAULT, required=True,
        help="Root directory containing experimental condition folders",
    )
    parser.add_argument(
        "--output-dir", type=str, default=OUTPUT_DIR_DEFAULT,
        help="Output directory for results and figures",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Discover images and print summary without running analysis",
    )
    parser.add_argument(
        "--no-gpu", action="store_true",
        help="Disable GPU for Cellpose segmentation",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ── Step 1: Image Discovery ──
    logger.info("Step 1: Discovering images...")
    fov_df = discover_images(args.data_root)

    if len(fov_df) == 0:
        logger.error("No FOVs discovered. Check --data-root path.")
        sys.exit(1)

    # Save discovery summary
    discovery_path = os.path.join(args.output_dir, "discovered_fovs.csv")
    fov_df.to_csv(discovery_path, index=False)
    logger.info(f"Discovery saved to {discovery_path}")

    # Print summary
    logger.info("\nDiscovery Summary:")
    summary = fov_df.groupby(["cell_line", "genotype", "condition", "magnification"]).size()
    logger.info(f"\n{summary.to_string()}")

    if args.dry_run:
        logger.info("\n--dry-run: stopping after discovery.")
        return

    # ── Steps 2-6: Process each FOV ──
    logger.info("\nSteps 2-6: Processing FOVs...")
    use_gpu = not args.no_gpu

    all_measurements_10x = []
    all_measurements_60x = []

    for idx, row in fov_df.iterrows():
        logger.info(
            f"\n[{idx + 1}/{len(fov_df)}] {row['condition_folder']} / {row['fov_name']}"
        )
        try:
            cell_df = process_fov(row, args.output_dir, gpu=use_gpu)
            if cell_df is not None:
                if row["magnification"] == "10x":
                    all_measurements_10x.append(cell_df)
                else:
                    all_measurements_60x.append(cell_df)
        except Exception as e:
            logger.error(f"  FAILED: {e}", exc_info=True)
            continue

    # Combine measurements
    df_10x = pd.concat(all_measurements_10x, ignore_index=True) if all_measurements_10x else pd.DataFrame()
    df_60x = pd.concat(all_measurements_60x, ignore_index=True) if all_measurements_60x else pd.DataFrame()

    # Save measurement CSVs
    if len(df_10x) > 0:
        path_10x = os.path.join(args.output_dir, "measurements_10x.csv")
        df_10x.to_csv(path_10x, index=False)
        logger.info(f"Saved {len(df_10x)} cells (10x) to {path_10x}")

    if len(df_60x) > 0:
        path_60x = os.path.join(args.output_dir, "measurements_60x.csv")
        df_60x.to_csv(path_60x, index=False)
        logger.info(f"Saved {len(df_60x)} cells (60x) to {path_60x}")

    # ── Step 7: Statistics ──
    logger.info("\nStep 7: Running statistics...")
    run_statistics(df_10x, df_60x, args.output_dir)

    # ── Step 8: Figures ──
    logger.info("\nStep 8: Generating figures...")

    # Use 10x data for overview plots (more FOVs), fall back to 60x
    df_all = pd.concat([df_10x, df_60x], ignore_index=True) if len(df_10x) > 0 or len(df_60x) > 0 else pd.DataFrame()

    if len(df_all) > 0:
        plot_mdk_intensity_barplot(df_all, args.output_dir)
        plot_mdk_timecourse(df_all, args.output_dir)
        plot_interaction(df_all, args.output_dir)

    if len(df_60x) > 0:
        plot_perinuclear_enrichment(df_60x, args.output_dir)

    logger.info("\nPipeline complete!")
    logger.info(f"Results in: {args.output_dir}")


if __name__ == "__main__":
    main()
