"""Post-process matplotlib SVGs to add CSS role classes to text elements.

Classification uses matplotlib's SVG group IDs as the primary signal:
  xtick_*, ytick_* -> fig-tick-label
  legend_*         -> fig-annotation
  matplotlib.axis_* with rotated text -> fig-axis-label
  bold text in text_* groups -> fig-title
  everything else  -> fig-data-label

Text inside <defs>, <clipPath>, <mask>, <pattern>, <symbol> is skipped.
Existing class attributes are preserved (role classes are appended).
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

from lxml import etree

SVG_NS = "http://www.w3.org/2000/svg"

log = logging.getLogger(__name__)

_SKIP_CONTAINERS = frozenset(
    [
        f"{{{SVG_NS}}}defs",
        f"{{{SVG_NS}}}clipPath",
        f"{{{SVG_NS}}}mask",
        f"{{{SVG_NS}}}pattern",
        f"{{{SVG_NS}}}symbol",
    ]
)

_ROTATE_RE = re.compile(r"rotate\(\s*(-?\d+\.?\d*)")


def _is_inside_skip_container(el) -> bool:
    """Check if element is inside a defs/clipPath/mask/pattern/symbol."""
    parent = el.getparent()
    while parent is not None:
        if parent.tag in _SKIP_CONTAINERS:
            return True
        parent = parent.getparent()
    return False


def _ancestor_group_ids(el) -> list[str]:
    """Collect all ancestor group IDs from element up to root."""
    ids = []
    parent = el.getparent()
    while parent is not None:
        gid = parent.get("id", "")
        if gid:
            ids.append(gid)
        parent = parent.getparent()
    return ids


def _has_rotation(el) -> bool:
    """Check if element or any ancestor has a rotate transform."""
    current = el
    while current is not None:
        transform = current.get("transform", "")
        m = _ROTATE_RE.search(transform)
        if m and abs(float(m.group(1))) > 45:
            return True
        current = current.getparent()
    return False


def _is_bold(el) -> bool:
    """Check if text element has bold font-weight."""
    style = el.get("style", "")
    m = re.search(r"font-weight\s*:\s*(\w+)", style)
    if m:
        val = m.group(1)
        return val == "bold" or (val.isdigit() and int(val) >= 700)
    weight = el.get("font-weight", "")
    return weight == "bold" or (weight.isdigit() and int(weight) >= 700)


def _classify_text(text_el) -> str:
    """Classify a single <text> element by its ancestor group IDs."""
    ancestor_ids = _ancestor_group_ids(text_el)

    # Strip panel prefixes (e.g., "panel-A-xtick_1" -> "xtick_1")
    stripped = []
    for gid in ancestor_ids:
        parts = gid.split("-")
        if len(parts) >= 3 and parts[0] == "panel":
            stripped.append("-".join(parts[2:]))
        else:
            stripped.append(gid)

    # Rule 1: inside xtick_* or ytick_* -> tick label
    for gid in stripped:
        if re.match(r"[xy]tick_\d+", gid):
            return "fig-tick-label"

    # Rule 2: inside legend_* -> legend
    for gid in stripped:
        if gid.startswith("legend"):
            return "fig-legend"

    # Rule 3: inside matplotlib.axis_* with rotation -> axis label
    in_axis = any(gid.startswith("matplotlib.axis") for gid in stripped)
    if in_axis:
        if _has_rotation(text_el):
            return "fig-axis-label"
        in_tick = any(re.match(r"[xy]tick_\d+", gid) for gid in stripped)
        if not in_tick:
            return "fig-axis-label"

    # Rule 4: bold text in text_* group -> title
    in_text_group = any(re.match(r"text_\d+", gid) for gid in stripped)
    if in_text_group and _is_bold(text_el):
        return "fig-title"

    # Rule 5: bold text anywhere -> annotation
    if _is_bold(text_el):
        return "fig-annotation"

    # Default: data label
    return "fig-data-label"


def _append_class(el, cls: str):
    """Append a CSS class to an element, preserving existing classes."""
    existing = el.get("class", "").strip()
    if existing:
        el.set("class", f"{existing} {cls}")
    else:
        el.set("class", cls)


def classify_svg_text(svg_path: Path) -> Path:
    """Post-process an SVG to add CSS role classes to <text> elements.

    Modifies the SVG file in-place. Returns the path.
    """
    svg_path = Path(svg_path)
    svg = etree.parse(str(svg_path)).getroot()

    texts = list(svg.iter(f"{{{SVG_NS}}}text"))

    renderable = [t for t in texts if not _is_inside_skip_container(t)]
    if not renderable:
        paths = list(svg.iter(f"{{{SVG_NS}}}path"))
        if paths:
            log.warning(
                "%s: no <text> elements found but %d <path> elements exist — "
                "text may be rendered as paths. Set matplotlib rcParams['svg.fonttype'] = 'none'.",
                svg_path.name,
                len(paths),
            )
        return svg_path

    classified = 0
    for text_el in renderable:
        role = _classify_text(text_el)
        _append_class(text_el, role)
        classified += 1

    tree = etree.ElementTree(svg)
    tree.write(str(svg_path), xml_declaration=True, encoding="utf-8", pretty_print=True)

    log.info("%s: classified %d text elements", svg_path.name, classified)
    return svg_path
