"""Rasterize an SVG file to PNG/PDF (used to re-derive raster formats from a
label-stamped SVG so PNG/PDF carry the same labels as the SVG).

Primary backend: rsvg-convert (librsvg).
Fallback: cairosvg, if installed.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def _rsvg_available() -> bool:
    return shutil.which("rsvg-convert") is not None


def svg_to_png(svg_path, out_path, dpi: int = 300) -> None:
    """Convert SVG to PNG at the given DPI."""
    svg_path = Path(svg_path)
    out_path = Path(out_path)
    if _rsvg_available():
        subprocess.run(
            [
                "rsvg-convert",
                "--format=png",
                f"--dpi-x={dpi}",
                f"--dpi-y={dpi}",
                "--background-color=white",
                "-o",
                str(out_path),
                str(svg_path),
            ],
            check=True,
        )
    else:
        import cairosvg

        cairosvg.svg2png(url=str(svg_path), write_to=str(out_path), dpi=dpi)


def svg_to_pdf(svg_path, out_path) -> None:
    """Convert SVG to PDF."""
    svg_path = Path(svg_path)
    out_path = Path(out_path)
    if _rsvg_available():
        subprocess.run(
            [
                "rsvg-convert",
                "--format=pdf",
                "-o",
                str(out_path),
                str(svg_path),
            ],
            check=True,
        )
    else:
        import cairosvg

        cairosvg.svg2pdf(url=str(svg_path), write_to=str(out_path))
