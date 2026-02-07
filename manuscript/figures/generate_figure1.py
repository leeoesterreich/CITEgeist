#!/usr/bin/env python3
"""
Figure 1: CITEgeist Pipeline Overview

ALL PANELS ARE SCHEMATICS - use SVG files from output/schematics/
- Panel A: figure1_panel_a_pipeline.svg
- Panel B: figure1_panel_b_spatial_stats.svg
- Panel C: figure1_panel_c_resolution.svg

Generate SVGs with: python svg_schematics.py
Combine in Illustrator for final figure.
"""

from pathlib import Path

OUTPUT_DIR = Path(__file__).parent / "output"
SCHEMATIC_DIR = OUTPUT_DIR / "schematics"

def main():
    print("=" * 60)
    print("Figure 1: CITEgeist Pipeline Overview")
    print("=" * 60)
    print("\nThis figure contains ONLY schematic panels.")
    print("Use the SVG files for professional quality:")
    print()

    svgs = [
        ("Panel A", "figure1_panel_a_pipeline.svg", "CITEgeist Modular Pipeline"),
        ("Panel B", "figure1_panel_b_spatial_stats.svg", "Spatial Statistics Foundation"),
        ("Panel C", "figure1_panel_c_resolution.svg", "Resolution-Agnostic Design"),
    ]

    for panel, filename, desc in svgs:
        filepath = SCHEMATIC_DIR / filename
        status = "✓ EXISTS" if filepath.exists() else "✗ MISSING"
        print(f"  {panel}: {filename}")
        print(f"          {desc}")
        print(f"          {status}")
        print()

    print("To regenerate SVGs: python svg_schematics.py")
    print("\nTo create final figure:")
    print("  1. Open SVGs in Illustrator/Inkscape")
    print("  2. Arrange panels A, B, C in desired layout")
    print("  3. Export as figure1_pipeline_overview.pdf")


if __name__ == "__main__":
    main()
