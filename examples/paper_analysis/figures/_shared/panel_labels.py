"""Shared panel-letter-label stamping for composites and supplementals.

Stamps <text id="panel-label-{L}"> onto any SVG root, converting physical-mm
panel rects into that SVG's user units. Pure: lxml + constants only.
"""
from __future__ import annotations

import re
from typing import Dict, List, Optional

from lxml import etree

SVG_NS = "http://www.w3.org/2000/svg"
MIN_LABEL_X_MM = 1.0  # leftmost label must sit >= 1mm inside the canvas

_DIM_RE = re.compile(r"^([\d.]+)\s*(mm|pt|px|in|cm)?$")
_UNIT_MM = {"mm": 1.0, "pt": 25.4 / 72.0, "px": 25.4 / 96.0, "in": 25.4, "cm": 10.0, "": 1.0}
_PT_PER_MM = 72.0 / 25.4


def mm_per_user_unit(svg_root: etree._Element) -> float:
    """Physical mm per SVG user unit, from the root width attr vs viewBox width.

    Composite canvases use width="Wmm" + viewBox "0 0 W H" -> 1.0.
    Matplotlib supp SVGs use width="Wpt" + viewBox "0 0 W H" -> 25.4/72.
    Returns 1.0 when metadata is missing (assume user units already mm).
    """
    vb = (svg_root.get("viewBox") or "").split()
    m = _DIM_RE.match((svg_root.get("width") or "").strip())
    if not m or len(vb) != 4:
        return 1.0
    width_mm = float(m.group(1)) * _UNIT_MM.get(m.group(2) or "", 1.0)
    vb_w = float(vb[2])
    return width_mm / vb_w if vb_w else 1.0


def _detect_grid(panels, tolerance_mm=8.0):
    """Detect row/column structure from panel positions.

    Returns (rows, columns) where each is a list of lists of panel dicts.
    Panels within tolerance_mm of each other's y (rows) or x (columns) are grouped.
    """
    # Detect rows
    sorted_by_y = sorted(panels, key=lambda p: p["y_mm"])
    rows = []
    current_row = [sorted_by_y[0]]
    for p in sorted_by_y[1:]:
        if abs(p["y_mm"] - current_row[0]["y_mm"]) < tolerance_mm:
            current_row.append(p)
        else:
            rows.append(sorted(current_row, key=lambda p: p["x_mm"]))
            current_row = [p]
    rows.append(sorted(current_row, key=lambda p: p["x_mm"]))

    # Detect columns
    sorted_by_x = sorted(panels, key=lambda p: p["x_mm"])
    cols = []
    current_col = [sorted_by_x[0]]
    for p in sorted_by_x[1:]:
        if abs(p["x_mm"] - current_col[0]["x_mm"]) < tolerance_mm:
            current_col.append(p)
        else:
            cols.append(current_col)
            current_col = [p]
    cols.append(current_col)

    return rows, cols


def stamp_panel_labels(
    svg_root: etree._Element,
    panels: List[dict],
    label_config: dict,
    label_positions: Optional[List[dict]] = None,
    mm_per_unit: Optional[float] = None,
) -> None:
    """Stamp <text id="panel-label-{L}"> for each panel onto svg_root.

    ``panels`` rects are physical mm. Coordinates and font-size are converted to
    the SVG's user units via ``mm_per_unit`` (auto-derived from svg_root if None).
    """
    style = label_config.get("style", "lowercase_bold")
    if style == "none":
        return
    font_pt = label_config.get("font_size_pt", 10)
    off_x = label_config.get("offset_x_mm", 0.0)
    off_y = label_config.get("offset_y_mm", 0.0)

    mmpu = mm_per_user_unit(svg_root) if mm_per_unit is None else mm_per_unit
    to_u = lambda mm: mm / mmpu  # noqa: E731
    font_u = (font_pt / _PT_PER_MM) / mmpu
    min_x_u = to_u(MIN_LABEL_X_MM)
    min_y_u = font_u + to_u(0.5)

    precomputed: Dict[str, tuple] = {p["label"]: (p["x_mm"], p["y_mm"]) for p in (label_positions or [])}
    col_x: Dict[str, float] = {}
    row_y: Dict[str, float] = {}
    if not precomputed:
        rows, cols = _detect_grid(panels)
        for col in cols:
            min_x = min(p["x_mm"] for p in col)
            for p in col:
                col_x[p["label"]] = min_x + off_x
        for row in rows:
            min_y = min(p["y_mm"] for p in row)
            for p in row:
                row_y[p["label"]] = min_y + off_y

    for panel_def in panels:
        label = panel_def["label"]
        if "_cont" in label.lower() or panel_def.get("skip_label"):
            continue
        if precomputed:
            if label not in precomputed:
                continue
            x_mm, y_mm = precomputed[label]
        else:
            x_mm = col_x.get(label, panel_def["x_mm"] + off_x)
            y_mm = row_y.get(label, panel_def["y_mm"] + off_y)

        x_u = max(min_x_u, to_u(x_mm))
        y_u = max(min_y_u, to_u(y_mm))
        display = label.upper() if style == "uppercase_bold" else label.lower()

        text = etree.SubElement(svg_root, f"{{{SVG_NS}}}text")
        text.set("id", f"panel-label-{label}")
        text.set("x", f"{x_u:.2f}")
        text.set("y", f"{y_u:.2f}")
        text.set(
            "style",
            f"font-family:Arial;font-size:{font_u:.2f};" f"font-weight:bold;dominant-baseline:alphabetic",
        )
        text.text = display
