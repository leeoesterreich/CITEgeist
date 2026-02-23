# Pseudo-Visium Simulation Method Comparison

Date: 2026-02-11

## Scope
This note compares:
- **Current CITEgeist Xenium -> pseudo-Visium generation** in this repository.
- **External method** (from provided paper excerpt): 65 um x 65 um bins with 35 um gaps, then count aggregation.

Goal of comparison: determine which approach is more correct for **simulating 10x Visium-like spots**.

## Method A (Current CITEgeist implementation)
Implementation files:
- `Benchmarking/xenium_pseudovisium/src/create_pseudo_spots.py`
- `Benchmarking/xenium_pseudovisium/src/split_regions.py`

Key geometry/assignment details:
- Spot diameter: **55 um** (`VISIUM_SPOT_DIAMETER = 55.0`).
- Center spacing (pitch): **100 um** (`VISIUM_CENTER_SPACING = 100.0`).
- Grid type: **hexagonal/honeycomb**, with odd-row x-offset.
- Cell assignment: nearest spot center, retained if within radius **27.5 um**.
- Spot counts: sum cell-level counts per spot.
- Additional practical filter: spots with fewer than 3 cells are removed (`min_cells=3` default in region generation pipeline).

Interpretation:
- This is a Visium-like circular capture model on a hex lattice with realistic pitch.

## Method B (External paper excerpt)
Provided description:
- Divide tissue into **65 um x 65 um square bins**.
- Use **35 um gaps** between bins (effective 100 um pitch).
- Aggregate detected counts in each bin to simulate low-resolution spots.

Interpretation:
- This is a square-bin downsampling strategy with matched pitch, but not Visium circular/hex geometry.

## Judgment
For the specific goal of simulating **Visium-like spots**, **Method A (current CITEgeist)** is more correct because it matches:
- Circular spot footprint (55 um diameter).
- Hexagonal spatial arrangement.
- Visium-like center spacing.

Method B is reasonable for generic low-resolution spatial aggregation, but is less geometrically faithful to Visium capture design.

## Caveat
Current CITEgeist generation includes a quality-control style filter (`min_cells=3`), which is useful for benchmark stability but is not a literal property of Visium chemistry. This should be disclosed when describing benchmark construction.

## Suggested manuscript wording
"Pseudo-Visium spots were generated from Xenium single-cell coordinates using a Visium-like hexagonal grid (55 um diameter, 100 um center spacing). Cells were assigned to the nearest spot center within a 27.5 um radius and counts were aggregated per spot. Spots with fewer than 3 cells were excluded."
