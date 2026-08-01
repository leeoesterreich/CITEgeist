"""Shared utilities for figure panel generation.

Provides:
- panel_argparser()    — argparse setup with --panel-sizes and --panels-only
- load_panel_sizes()   — load JSON panel-size manifest
- PanelContext         — context manager: style, figure creation, save, restore
- save_all_panels()    — convenience loop over (label, draw_fn) specs
"""

from __future__ import annotations

import argparse
import json
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Callable, Dict, Generator, List, Optional, Sequence, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# ---------------------------------------------------------------------------
# Sibling imports (same _shared directory is on sys.path at call time, but
# guard against being imported from elsewhere by temporarily adjusting path).
# ---------------------------------------------------------------------------
_SHARED_DIR = Path(__file__).resolve().parent


def _import_shared(name: str):
    """Import a module from the _shared directory."""
    if str(_SHARED_DIR) not in sys.path:
        sys.path.insert(0, str(_SHARED_DIR))
    import importlib

    return importlib.import_module(name)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def panel_argparser(description: str = "") -> argparse.ArgumentParser:
    """Return an ArgumentParser pre-loaded with common panel generation flags.

    Flags added:
    - ``--panel-sizes PATH``  Path to a JSON file with per-panel size dicts
                              (keys ``w_mm``, ``h_mm``).  Default: None.
    - ``--panels-only``       When set, skip any post-panel assembly steps.
                              Default: False.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--panel-sizes",
        type=Path,
        default=None,
        metavar="PATH",
        help="Path to JSON file mapping panel labels to {w_mm, h_mm} dicts.",
    )
    parser.add_argument(
        "--panels-only",
        action="store_true",
        default=False,
        help="Generate individual panels only; skip composite assembly.",
    )
    return parser


def load_panel_sizes(args) -> Optional[Dict]:
    """Read the panel-sizes JSON manifest, defaulting to a generator-local file.

    Resolution order:
      1. ``--panel-sizes PATH`` if explicitly given.
      2. ``panel_sizes.json`` sitting next to the running generator script
         (``sys.argv[0]``).  It ships next to each generator, so a bare
         ``python generate.py`` run picks up the correct target sizes
         automatically.
      3. None — but only after a LOUD warning, because falling through to
         PanelContext's ``default_figsize`` produces wrong-aspect panels that the
         assembly silently shrinks (the failure this guard exists to prevent).
    """
    if args.panel_sizes is not None:
        path = Path(args.panel_sizes)
        with path.open() as fh:
            return json.load(fh)

    # Idiot-proof default: auto-discover the generator's sibling panel_sizes.json
    # so omitting --panel-sizes can never silently fall back to default figsizes.
    local = Path(sys.argv[0]).resolve().parent / "panel_sizes.json"
    if local.exists():
        print(
            f"[panel-sizes] --panel-sizes not given; auto-loading {local}",
            flush=True,
        )
        with local.open() as fh:
            return json.load(fh)

    print(
        "WARNING: no --panel-sizes given and no panel_sizes.json found next to "
        f"{Path(sys.argv[0]).name}. Panels will render at PanelContext's "
        "default_figsize and will NOT match the assembly layout (wrong aspect "
        "ratio -> shrunk panels + whitespace). Pass --panel-sizes <path>.",
        file=sys.stderr,
        flush=True,
    )
    return None


def place_annotation(ax, x, y, text, margin: float = 0.02, **kwargs):
    """Place a text annotation in axes coords, clamped to a safe interior zone.

    ``x``, ``y`` are axes-fraction coordinates (0..1). They are clamped to
    ``[margin, 1 - margin]`` so the annotation never sits within ``margin``
    of an axes edge — eliminating bbox clipping past the panel rectangle.

    Defaults: anchor is chosen so the annotation grows toward the nearest
    inside corner (e.g. top-left ``(0.02, 0.97)`` → ha="left", va="top"),
    fontsize 8, transform ``ax.transAxes``.

    Returns the matplotlib Text object.
    """
    cx = max(margin, min(1.0 - margin, x))
    cy = max(margin, min(1.0 - margin, y))
    ha = kwargs.pop("ha", "left" if cx < 0.5 else "right")
    va = kwargs.pop("va", "bottom" if cy < 0.5 else "top")
    fontsize = kwargs.pop("fontsize", 8)
    transform = kwargs.pop("transform", ax.transAxes)
    return ax.text(
        cx,
        cy,
        text,
        ha=ha,
        va=va,
        fontsize=fontsize,
        transform=transform,
        **kwargs,
    )


def _auto_position(handles, w_mm, h_mm) -> str:
    """Pick legend position from {above, below, right} based on panel size + item count.

    Heuristic:
      - aspect >= 1.3 and n > 6 → above (wide, many items)
      - aspect < 0.8            → right (tall panel)
      - n <= 6                  → right (few items)
      - else                    → below (wide/square with many items)

    'left' is never selected automatically. Returns 'below' if dims or handles missing.
    """
    n = len(handles) if handles else 0
    if w_mm is None or h_mm is None or n == 0:
        return "below"
    aspect = w_mm / h_mm
    if aspect >= 1.3 and n > 6:
        return "above"
    if aspect < 0.8:
        return "right"
    if n <= 6:
        return "right"
    return "below"


def _fit_outside_legend(
    fig,
    handles,
    labels,
    loc,
    kw,
    *,
    auto_ncol: bool,
    n: int,
    floor_pt: float,
    axes_reserve_frac: float = 0.55,
    width_margin_px: float = 2.0,
):
    """Build a fig-level outside legend that fits BOTH the figure width and a
    height budget. Searches ncol (only when auto_ncol) then steps fontsize down
    to floor_pt. Returns the final legend; warns (does not crash) if it cannot
    fit at the floor. Its own rebuilds leave exactly one legend in fig.legends.
    """
    import warnings

    # Resolve the starting fontsize to POINTS. kw["fontsize"] may be a named
    # size ("small"/"large"/...) — resolve via FontProperties, never float() it.
    raw_fs = kw.get("fontsize", None)
    if raw_fs is None:
        raw_fs = plt.rcParams.get("legend.fontsize", 7)
    try:
        base_fs = float(raw_fs)
    except (TypeError, ValueError):
        base_fs = float(FontProperties(size=raw_fs).get_size_in_points())

    # Fontsize ladder: base down to floor in 1pt steps (floor always included).
    fontsizes = []
    fs = base_fs
    while fs > floor_pt:
        fontsizes.append(round(fs, 1))
        fs -= 1.0
    fontsizes.append(round(float(floor_pt), 1))

    # ncol search range: the FULL [1, n], because the height budget may require
    # MORE than min(n,4) columns (a short, wide panel with many labels needs
    # extra columns to stay short). Order candidates by closeness to the
    # conventional min(n,4), tie-breaking toward more columns (shorter legend).
    if auto_ncol:
        start = min(n, 4)
        ncol_candidates = sorted(range(1, n + 1), key=lambda c: (abs(c - start), -c))
    else:
        ncol_candidates = [kw.get("ncol", 1)]

    current = None

    def _rebuild(ncol, fontsize):
        nonlocal current
        if current is not None:
            current.remove()
            current = None
        k = dict(kw)
        k["ncol"] = ncol
        k["fontsize"] = fontsize
        current = fig.legend(handles=handles, labels=labels, loc=loc, **k)
        fig.canvas.draw()  # constrained layout repositions; refresh extents
        return current

    for fontsize in fontsizes:
        for ncol in ncol_candidates:
            leg = _rebuild(ncol, fontsize)
            bb = leg.get_window_extent()
            width_budget = fig.bbox.width - width_margin_px
            height_budget = fig.bbox.height * (1.0 - axes_reserve_frac)
            if bb.width <= width_budget and bb.height <= height_budget:
                return leg  # fits both dimensions

    # No (fontsize, ncol) combination satisfied both budgets.  Build the
    # least-bad fallback: floor fontsize with the MAXIMUM number of columns
    # (= min height), which is the last ncol in ncol_candidates order that
    # we haven't already settled on.  Rebuild to ensure fig.legends contains
    # exactly the one legend we return.
    max_ncol = n if auto_ncol else ncol_candidates[0]
    fallback = _rebuild(max_ncol, round(float(floor_pt), 1))
    warnings.warn(
        "place_legend: legend does not fit panel width/height even at the "
        f"font floor ({floor_pt}pt); it may clip. Use fewer/shorter labels.",
        UserWarning,
        stacklevel=3,
    )
    return fallback


# ---------------------------------------------------------------------------
# Text overflow auto-fit guard
# ---------------------------------------------------------------------------


def fit_text_to_figure(
    fig,
    floor_pt: float = 6.0,
    margin_px: float = 2.0,
) -> int:
    """Detect text artists that overflow the figure bbox and wrap/shrink them.

    For each axes title, axis label, ``ax.texts`` entry, and ``fig.texts``
    entry whose ``get_window_extent()`` lies outside the figure pixel box
    (allowing ``margin_px`` slack), the function tries:

    1. **Wrap** (titles/labels with >1 word): iteratively insert line-breaks
       across 2, 3, ... lines using :func:`textwrap.fill`.  If wrapping makes
       it fit, move on.
    2. **Shrink** in 0.5 pt steps down to ``floor_pt``, redrawing after each
       step.  If a wrapped form was tried but did not fully fix the overflow,
       the text is reset to the original single-line string before shrinking
       so the two strategies don't compound.
    3. If still overflowing at the floor: emit a :class:`UserWarning` and
       leave the artist as-is.

    Always calls ``fig.canvas.draw()`` before measuring, and again after every
    mutation so constrained layout repositions properly.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to inspect.
    floor_pt : float
        Minimum font size (points) the function will shrink to.  Default 6.0.
    margin_px : float
        Slack (pixels) allowed outside the figure box on each edge.  Default 2.0.

    Returns
    -------
    int
        Number of text artists that were adjusted (wrapped or shrunk).
    """
    import textwrap
    import warnings

    fig.canvas.draw()
    fig_w = fig.bbox.width
    fig_h = fig.bbox.height

    def _overflows(artist) -> bool:
        bb = artist.get_window_extent()
        return bb.x0 < -margin_px or bb.x1 > fig_w + margin_px or bb.y0 < -margin_px or bb.y1 > fig_h + margin_px

    # Collect all text artists: axes titles, axis labels, ax.texts, fig.texts.
    artists = []
    for ax in fig.get_axes():
        artists.append(ax.title)
        artists.append(ax.xaxis.label)
        artists.append(ax.yaxis.label)
        artists.extend(ax.texts)
    artists.extend(fig.texts)

    # Filter to only non-empty artists that actually overflow.
    overflowing = [a for a in artists if a.get_text().strip() and _overflows(a)]

    n_adjusted = 0

    for artist in overflowing:
        original_text = artist.get_text()
        original_fs = artist.get_fontsize()
        words = original_text.split()

        fixed = False

        # --- Strategy 1: wrap (only when there are at least 2 words) -------
        if len(words) > 1:
            for nlines in range(2, len(words) + 1):
                char_width = max(1, -(-len(original_text) // nlines))  # ceil division
                wrapped = textwrap.fill(original_text, width=char_width)
                artist.set_text(wrapped)
                fig.canvas.draw()
                if not _overflows(artist):
                    fixed = True
                    break

        # --- Strategy 2: shrink fontsize ------------------------------------
        if not fixed:
            # Reset to single-line before shrinking (deterministic baseline).
            artist.set_text(original_text)
            fs = original_fs
            while fs - 0.5 >= floor_pt:
                fs = round(fs - 0.5, 1)
                artist.set_fontsize(fs)
                fig.canvas.draw()
                if not _overflows(artist):
                    fixed = True
                    break
            # Ensure we hit the floor exactly if we haven't crossed it yet.
            if not fixed and fs > floor_pt:
                artist.set_fontsize(floor_pt)
                fig.canvas.draw()
                fixed = not _overflows(artist)

        if fixed:
            n_adjusted += 1
        else:
            # Floor reached and still overflowing — warn, leave as-is.
            warnings.warn(
                f"PanelContext: text {original_text!r} overflows the panel even at the "
                f"{floor_pt}pt floor; it may clip.",
                UserWarning,
                stacklevel=2,
            )

    return n_adjusted


def place_legend(
    ax,
    handles=None,
    labels=None,
    position: Optional[str] = None,
    title=None,
    fontsize=None,
    ncol: Optional[int] = None,
    w_mm: Optional[float] = None,
    h_mm: Optional[float] = None,
    **kwargs,
):
    """Place a legend using ``fig.legend(loc='outside ...')`` for clip-free placement.

    Constrained layout automatically shrinks axes to make room for the legend
    within the fixed panel canvas.  For panels without constrained layout,
    falls back to inside placement with a warning.

    Parameters
    ----------
    ax : matplotlib Axes
        Target axes whose handles/labels are collected.  ``fig`` is obtained
        via ``ax.figure``.
    handles, labels : optional
        Explicit legend entries.  When *None*, collected from *ax*.
    position : str
        One of ``"right"`` (default), ``"above"``, ``"below"``, ``"bottom"``
        (alias for ``"below"``), ``"inside"``, ``"none"``.
    title, fontsize, ncol : forwarded to the legend call.
    **kwargs : dict
        Extra keyword arguments forwarded to the legend call (e.g.
        ``handlelength``, ``handletextpad``, ``borderpad``).

    Returns
    -------
    matplotlib.legend.Legend or None
    """
    import warnings

    fig = ax.figure

    # -- resolve position ------------------------------------------------
    if position is None:
        position = "auto"
    if position == "auto":
        # Resolve panel dims from active context if not provided
        if w_mm is None or h_mm is None:
            ctx = PanelContext._active
            if ctx is not None:
                if w_mm is None:
                    w_mm = getattr(ctx, "_w_mm", None)
                if h_mm is None:
                    h_mm = getattr(ctx, "_h_mm", None)
        # Need handles to count items; collect early if needed
        if handles is None:
            _h_tmp, _l_tmp = ax.get_legend_handles_labels()
        else:
            _h_tmp = handles
        position = _auto_position(_h_tmp, w_mm, h_mm)
    if position == "bottom":
        position = "below"
    if position == "none":
        return None

    # -- collect handles/labels from the *target* axes only ---------------
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
    elif labels is None:
        labels = [h.get_label() for h in handles]

    if not handles:
        return None

    # -- constrained-layout guard ----------------------------------------
    if position in ("right", "above", "below"):
        try:
            has_cl = fig.get_constrained_layout()
        except AttributeError:
            has_cl = False
        if not has_cl:
            warnings.warn(
                f"place_legend: position='{position}' requires constrained_layout "
                f"but this figure doesn't have it. Falling back to 'inside'.",
                UserWarning,
                stacklevel=2,
            )
            position = "inside"

    # -- build common kwargs ---------------------------------------------
    kw: Dict = {"frameon": False}
    if title is not None:
        kw["title"] = title
    if fontsize is not None:
        kw["fontsize"] = fontsize
    kw.update(kwargs)

    n = len(handles)
    loc_map = {
        "right": "outside right",
        "above": "outside upper center",
        "below": "outside lower center",
    }

    if position in ("above", "below"):
        # Resolve the per-panel legend font floor from the active context.
        floor_pt = 7.0
        ctx = PanelContext._active
        if ctx is not None and getattr(ctx, "_fonts", None):
            floor_pt = float(ctx._fonts.get("legend_pt", 7.0))
        auto_ncol = ncol is None
        if not auto_ncol:
            kw["ncol"] = ncol
        return _fit_outside_legend(
            fig,
            handles,
            labels,
            loc_map[position],
            kw,
            auto_ncol=auto_ncol,
            n=n,
            floor_pt=floor_pt,
        )

    # right: single-column outside legend (vertical overflow handled by CL)
    if position == "right":
        kw["ncol"] = ncol if ncol is not None else 1
        return fig.legend(handles=handles, labels=labels, loc=loc_map["right"], **kw)

    # "inside" — use ax.legend so it attaches to the target axes
    kw["ncol"] = ncol if ncol is not None else 1
    return ax.legend(handles=handles, labels=labels, loc="best", **kw)


# ---------------------------------------------------------------------------
# compute_scaled_fonts
# ---------------------------------------------------------------------------


def compute_scaled_fonts(
    w_mm: float,
    h_mm: float,
    floor: float = 7.0,
    cap: float = 10.0,
    coefficient: float = 0.17,
) -> Dict[str, float]:
    """Compute font sizes scaled to panel physical dimensions.

    Uses the smaller dimension as the governing size — narrow panels get
    smaller fonts regardless of height, and vice versa.
    """
    base = min(w_mm, h_mm)
    tick = max(floor, min(base * coefficient, cap))
    return {
        "tick_label_pt": round(tick, 1),
        "axis_title_pt": round(tick + 1, 1),
        "legend_pt": round(tick, 1),
        "title_pt": round(tick + 1, 1),
        "annotation_pt": round(tick - 0.5, 1),
    }


# ---------------------------------------------------------------------------
# PanelContext
# ---------------------------------------------------------------------------


class PanelContext:
    """Context manager for generating a single figure panel.

    Parameters
    ----------
    label : str
        Panel label (e.g. ``"A"``).  Used to name output files
        ``panel_{label}.{fmt}``.
    panel_sizes : dict or None
        Mapping of label → ``{"w_mm": float, "h_mm": float}``.
        If the label is found the figsize is derived from it; otherwise
        *default_figsize* is used.
    output_dir : Path-like
        Directory where panel files are written.  Created if absent.
    default_figsize : (float, float)
        Fallback figure size in inches when *panel_sizes* is None or the
        label is absent.
    journal : str
        Journal profile passed to ``apply_style()`` (``"nature"`` / …).
    figonly : bool
        When True the context yields the ``Figure`` alone; when False it
        yields ``(fig, ax)``.
    constrained_layout : bool
        If False overrides the ``figure.constrained_layout.use`` rcParam
        set by ``apply_style()``.
    formats : sequence of str
        File formats to save (e.g. ``("jpg", "svg")``).

    Usage
    -----
    ::

        with PanelContext("A", panel_sizes, output_dir) as (fig, ax):
            ax.plot(...)
    """

    # Active context for legend placement to read panel dimensions.
    # Set in __enter__, cleared in __exit__ finally.
    _active: "PanelContext | None" = None

    def __init__(
        self,
        label: str,
        panel_sizes: Optional[Dict],
        output_dir: Union[str, Path],
        default_figsize: Tuple[float, float] = (5, 4),
        journal: str = "nature",
        figonly: bool = False,
        constrained_layout: bool = True,
        formats: Sequence[str] = ("jpg", "svg"),
        basename: Optional[str] = None,
        fit_overflow_text: bool = False,
        panel_label_config: Optional[Dict] = None,
    ) -> None:
        self.label = label
        self.panel_sizes = panel_sizes or {}
        self.output_dir = Path(output_dir)
        self.default_figsize = default_figsize
        self.journal = journal
        self.figonly = figonly
        self.constrained_layout = constrained_layout
        self.formats = list(formats)
        self.basename = basename
        self.fit_overflow_text = fit_overflow_text
        self._saved_rc: Optional[Dict] = None
        self._fig: Optional[plt.Figure] = None
        self._panel_title: Optional[str] = None
        self._fonts: Optional[Dict] = None
        self._panel_label_axes: Dict = {}
        self._panel_label_config = panel_label_config or {
            "style": "uppercase_bold",
            "font_size_pt": 14,
            "offset_x_mm": -4.0,
            "offset_y_mm": -1.0,
        }

        if self.panel_sizes:
            panel_spec = self.panel_sizes.get(self.label, {})
            if isinstance(panel_spec, dict):
                self._panel_title = panel_spec.get("title")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _compute_figsize(self) -> Tuple[float, float]:
        spec = self.panel_sizes.get(self.label)
        if spec is not None:
            return spec["w_mm"] / 25.4, spec["h_mm"] / 25.4
        return self.default_figsize

    def _apply_style(self) -> None:
        style = _import_shared("style")
        style.apply_style(journal=self.journal)

    def _classify_svg(self, path: Path) -> None:
        svg_post = _import_shared("svg_postprocess")
        svg_post.classify_svg_text(path)

    # ------------------------------------------------------------------
    # Context manager protocol
    # ------------------------------------------------------------------

    def __enter__(self):
        # 1. Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # 2. Save current rcParams
        self._saved_rc = dict(plt.rcParams)

        # 3. Apply journal style
        self._apply_style()

        # 4. Compute and apply scaled fonts
        spec = self.panel_sizes.get(self.label) if self.panel_sizes else None
        if spec is not None and isinstance(spec, dict):
            w_mm = spec.get("w_mm", self.default_figsize[0] * 25.4)
            h_mm = spec.get("h_mm", self.default_figsize[1] * 25.4)
        else:
            w_mm = self.default_figsize[0] * 25.4
            h_mm = self.default_figsize[1] * 25.4

        if spec is not None and isinstance(spec, dict) and "fonts" in spec and spec["fonts"]:
            fonts = spec["fonts"]
        else:
            fonts = compute_scaled_fonts(w_mm, h_mm)
        self._fonts = fonts

        plt.rcParams["axes.labelsize"] = fonts.get("axis_title_pt", 8)
        plt.rcParams["xtick.labelsize"] = fonts.get("tick_label_pt", 7)
        plt.rcParams["ytick.labelsize"] = fonts.get("tick_label_pt", 7)
        plt.rcParams["legend.fontsize"] = fonts.get("legend_pt", 7)
        plt.rcParams["axes.titlesize"] = fonts.get("title_pt", 8)
        plt.rcParams["font.size"] = fonts.get("annotation_pt", 7)
        # Per-panel constrained_layout pad scaled by title font (0.5×/0.25× coeffs).
        # Replaces static 0.04/0.02 in style.py — sized to actual title height needing
        # clearance from outside legends.
        title_pt = fonts.get("title_pt", 8)
        plt.rcParams["figure.constrained_layout.h_pad"] = (title_pt * 0.5) / 72.0
        plt.rcParams["figure.constrained_layout.w_pad"] = (title_pt * 0.25) / 72.0

        # 5. Override constrained_layout if requested
        if not self.constrained_layout:
            plt.rcParams["figure.constrained_layout.use"] = False

        # 6. Compute figsize and create figure
        figsize = self._compute_figsize()
        self._figsize = figsize  # pinned slot size — restored in _save_files

        # Expose panel dimensions and register as active context (after all setup,
        # so a partial __enter__ never leaves a half-initialized active).
        self._w_mm = w_mm
        self._h_mm = h_mm
        PanelContext._active = self

        if self.figonly:
            self._fig = plt.figure(figsize=figsize)
            self._primary_ax = None
            return self._fig
        else:
            self._fig, ax = plt.subplots(figsize=figsize)
            self._primary_ax = ax
            return self._fig, ax

    def __exit__(self, exc_type, exc_val, exc_tb):
        try:
            if exc_type is None:
                # Apply manifest panel title on the primary axes.
                # Uses ax.set_title (not fig.suptitle) so constrained_layout
                # manages it below any outside legend — no collision.
                if self._panel_title and self._fig is not None:
                    ax = self._primary_ax
                    if ax is not None and not ax.get_title():
                        title_size = self._fonts.get("title_pt", 8) if self._fonts else 8
                        ax.set_title(self._panel_title, fontsize=title_size, pad=4)
                elif self._fig is not None:
                    has_any_title = any(ax.get_title() for ax in self._fig.get_axes())
                    if not has_any_title:
                        import warnings

                        warnings.warn(
                            f"Panel {self.label}: no title set (add to manifest panel_titles or call ax.set_title())",
                            UserWarning,
                            stacklevel=2,
                        )

                # Opt-in text overflow guard: wrap/shrink before save.
                if self.fit_overflow_text and self._fig is not None:
                    fit_text_to_figure(self._fig)

                # Success path — save files
                self._save_files()
        finally:
            # Always close the figure and restore rcParams
            if self._fig is not None:
                plt.close(self._fig)
                self._fig = None
            if self._saved_rc is not None:
                plt.rcParams.update(self._saved_rc)
                self._saved_rc = None
            PanelContext._active = None

        # Re-raise any exception
        return False

    def set_panel_label_axes(self, mapping: Dict) -> None:
        """Register {label: Axes} for post-save panel-label stamping.

        Uses the same PanelContext._active pattern as place_legend; call inside
        the with-block after the sub-axes exist.
        """
        self._panel_label_axes = dict(mapping)

    def _save_files(self) -> None:
        # Pin canvas to the slot size before saving.  Constrained layout can
        # silently expand the figure to accommodate out-of-axes artists (e.g.
        # legends placed below the axes with large bbox_to_anchor y-offsets).
        # Restoring the slot size here ensures SVGs are exactly the right
        # physical size without re-adding bbox_inches='tight'.
        if hasattr(self, "_figsize") and self._figsize is not None:
            self._fig.set_size_inches(self._figsize)

        stem = self.basename if self.basename else f"panel_{self.label}"

        # 1. SVG first (canonical; gets stamped).  Write it if requested or if
        #    panel labels are present and need to be burned into raster formats.
        want_svg = "svg" in self.formats
        svg_out = self.output_dir / f"{stem}.svg"
        need_svg_for_stamp = bool(self._panel_label_axes)
        if want_svg or need_svg_for_stamp:
            self._fig.savefig(svg_out, dpi=300, facecolor="white")
            if want_svg:
                try:
                    self._classify_svg(svg_out)
                except Exception as exc:  # noqa: BLE001
                    print(f"[generate_utils] SVG classification warning: {exc}")

        # 2. Stamp panel letter labels onto the SVG (post-draw axes positions).
        if self._panel_label_axes:
            from lxml import etree as _etree

            _stamp = _import_shared("panel_labels").stamp_panel_labels
            w_in, h_in = self._fig.get_size_inches()
            W_mm, H_mm = w_in * 25.4, h_in * 25.4
            # Use each axes' TIGHT bbox (includes its title + tick labels), not
            # get_position() (content box only). The panel letter must clear the
            # centered subplot title that sits ABOVE the content box, so we
            # anchor the label to the top-left of the full panel footprint.
            renderer = self._fig.canvas.get_renderer()
            inv = self._fig.transFigure.inverted()
            panels = []
            for L, ax in self._panel_label_axes.items():
                bb = ax.get_tightbbox(renderer).transformed(inv)  # figure fraction
                panels.append(
                    {
                        "label": L,
                        "x_mm": bb.x0 * W_mm,
                        "y_mm": (1.0 - bb.y1) * H_mm,  # SVG top-down: top edge
                        "width_mm": bb.width * W_mm,
                        "height_mm": bb.height * H_mm,
                    }
                )
            root = _etree.parse(str(svg_out)).getroot()
            _stamp(root, panels, self._panel_label_config)  # auto mm_per_unit
            _etree.ElementTree(root).write(str(svg_out), xml_declaration=True, encoding="utf-8")

        # 3. Other formats.  When labels were stamped, rasterize the STAMPED
        #    SVG so PNG/PDF carry labels.  Otherwise save from matplotlib
        #    directly (unchanged behavior).
        raster = _import_shared("svg_raster")
        for fmt in self.formats:
            if fmt == "svg":
                continue
            out_path = self.output_dir / f"{stem}.{fmt}"
            if self._panel_label_axes and fmt in ("png", "pdf"):
                if fmt == "png":
                    raster.svg_to_png(svg_out, out_path, dpi=300)
                else:
                    raster.svg_to_pdf(svg_out, out_path)
            else:
                self._fig.savefig(out_path, dpi=300, facecolor="white")

        # 4. If svg was only needed for stamping but not requested, remove it.
        if need_svg_for_stamp and not want_svg:
            svg_out.unlink(missing_ok=True)

        fmts_str = "+".join(self.formats)
        print(f"  {stem}: saved ({fmts_str}) -> {self.output_dir}")


# ---------------------------------------------------------------------------
# save_all_panels
# ---------------------------------------------------------------------------


def save_all_panels(
    panel_specs: List[Tuple[str, Callable]],
    panel_sizes: Optional[Dict],
    output_dir: Union[str, Path],
    default_figsize: Tuple[float, float] = (5, 4),
    journal: str = "nature",
    formats: Sequence[str] = ("jpg", "svg"),
) -> None:
    """Generate and save all panels from a list of (label, draw_fn) specs.

    Parameters
    ----------
    panel_specs : list of (str, callable)
        Each tuple is ``(label, draw_fn)`` where ``draw_fn(fig, ax)`` draws
        into the provided axes.
    panel_sizes : dict or None
        Passed to :class:`PanelContext`.
    output_dir : Path-like
        Destination directory for all panels.
    default_figsize : (float, float)
        Fallback size in inches.
    journal : str
        Journal style profile.
    formats : sequence of str
        Output formats for every panel.
    """
    for label, draw_fn in panel_specs:
        with PanelContext(
            label=label,
            panel_sizes=panel_sizes,
            output_dir=output_dir,
            default_figsize=default_figsize,
            journal=journal,
            figonly=False,
            formats=formats,
        ) as (fig, ax):
            draw_fn(fig, ax)
