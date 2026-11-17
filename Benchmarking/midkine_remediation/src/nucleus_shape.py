"""Per-nucleus shape table for one breast specimen (diagnostic E1; not part of any production path).

Replicates the production call `compute_spot_nuclei_counts_patchwise` (model/morphology/segmentation.py)
per spot patch, using the same StarDistSegmenter (2D_versatile_he, prob_thresh 0.005, scale = px/0.25, 20-500 um^2 filter),
same patch size (55 um spot + 25 um margin each side), same blank-patch skip, same in-spot test on StarDist points.
The only addition: skimage regionprops on the patch label image, kept per nucleus.
Package code is NOT modified. Output: <out_root>/<specimen>/nuclei_shape.csv (+ spot_counts.csv, meta.json).

Usage: nucleus_shape.py SPECIMEN OUT_ROOT [DATA_ROOT]
"""
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
import _paths as _cg  # noqa: E402

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

REPO = Path(str(_cg.REPO_ROOT))
DATA_ROOT = Path(str(_cg.root("DATA") / "GSE289326" / "processed_files"))
SPOT_UM = 55.0
MARGIN_UM = 25.0


def load_seg():
    # load segmentation.py by path so the heavy package __init__ is not imported
    p = REPO / "CITEgeist/model/morphology/segmentation.py"
    spec = importlib.util.spec_from_file_location("seg_mod", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules["seg_mod"] = m
    spec.loader.exec_module(m)
    return m


def main(specimen, out_root, data_root=DATA_ROOT, max_spots=None):
    from skimage.color import rgb2hed
    from skimage.measure import regionprops

    Image.MAX_IMAGE_PIXELS = None
    seg = load_seg()
    outs = Path(data_root) / specimen / "outs"
    sf = json.loads((outs / "spatial/scalefactors_json.json").read_text())
    px = SPOT_UM / float(sf["spot_diameter_fullres"])  # same as estimate_pixel_size_um
    pos = pd.read_csv(outs / "spatial/tissue_positions.csv")
    pos = pos[pos["in_tissue"] == 1].reset_index(drop=True)
    if max_spots:
        pos = pos.iloc[:max_spots]
    xy = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(float)  # (x, y) as obsm["spatial"]
    names = pos["barcode"].to_numpy()
    img = np.array(Image.open(outs / "spatial/tissue_fullres_image.tif").convert("RGB"))
    H, W = img.shape[:2]
    segr = seg.StarDistSegmenter(modality="he")

    r_um = SPOT_UM / 2
    half = int(np.ceil((r_um + MARGIN_UM) / px))
    r_px = r_um / px
    rows, counts = [], np.zeros(len(names), dtype=int)
    t0 = time.time()
    for i, (cx, cy) in enumerate(xy):
        if not (np.isfinite(cx) and np.isfinite(cy)):
            continue
        cxi, cyi = int(round(cx)), int(round(cy))
        y0, y1 = max(0, cyi - half), min(H, cyi + half)
        x0, x1 = max(0, cxi - half), min(W, cxi + half)
        if y1 - y0 < 10 or x1 - x0 < 10:
            continue
        patch = img[y0:y1, x0:x1]
        if patch.std() < 3.0:
            continue
        masks, cdf = segr.segment(patch, pixel_size_um=px)
        if len(cdf) == 0:
            continue
        props = regionprops(masks, intensity_image=rgb2hed(patch)[..., 0])
        if len(props) != len(cdf):
            print(f"WARN spot {names[i]}: {len(props)} regions vs {len(cdf)} points; skipped", flush=True)
            continue
        # kept labels (ascending) correspond one-to-one, in order, to the filtered StarDist points
        for pr, (py, pxx) in zip(props, cdf[["y_pixel", "x_pixel"]].to_numpy()):
            fx, fy = pxx + x0, py + y0
            if np.hypot(fx - cx, fy - cy) > r_px:
                continue
            counts[i] += 1
            rows.append(
                dict(
                    barcode=names[i], x_fullres=fx, y_fullres=fy, area_px=pr.area, area_um2=pr.area * px * px,
                    perimeter_px=pr.perimeter, eccentricity=pr.eccentricity, solidity=pr.solidity,
                    major_axis_px=pr.major_axis_length, minor_axis_px=pr.minor_axis_length,
                    mean_hema=pr.mean_intensity,
                )
            )
        if i % 500 == 0:
            print(f"{i}/{len(names)} spots, {len(rows)} nuclei, {time.time()-t0:.0f}s", flush=True)

    out = Path(out_root) / specimen
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "nuclei_shape.csv", index=False)
    pd.Series(counts, index=names, name="n_nuclei").rename_axis("barcode").to_csv(out / "spot_counts.csv")
    (out / "meta.json").write_text(json.dumps(dict(specimen=specimen, pixel_size_um=px, n_spots=len(names),
                                                   n_nuclei=len(rows), seconds=time.time() - t0, image_hw=[H, W])))
    print(f"done {specimen}: {len(rows)} nuclei in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    ms = int(sys.argv[4]) if len(sys.argv) > 4 else None
    main(sys.argv[1], sys.argv[2], Path(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != "-" else DATA_ROOT, ms)
