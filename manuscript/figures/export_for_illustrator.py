#!/usr/bin/env python3
"""
Export all figures as PNG and SVG for Illustrator editing.
"""
import subprocess
import shutil
from pathlib import Path

OUTPUT_DIR = Path(__file__).parent / "for_illustrator"
SOURCE_DIR = Path(__file__).parent / "output"

# Figure mapping (v5 manuscript - 5 main figures + supplementary)
FIGURES = {
    1: "figure1_pipeline_overview",
    2: "figure2_profile_discovery",
    3: "figure3_benchmarking",
    4: "figure4_midkine_esr1",
    5: "figure5_full_pipeline",
}

SUPP_FIGURES = {
    "S2": "supp_figure2_de_pathway",
}

def main():
    # Ensure output directories exist for main figures
    for i in range(1, 6):
        (OUTPUT_DIR / f"Figure{i}").mkdir(parents=True, exist_ok=True)

    # Ensure output directories exist for supplementary figures
    for supp_id in SUPP_FIGURES.keys():
        (OUTPUT_DIR / f"Supp{supp_id}").mkdir(parents=True, exist_ok=True)

    # First, regenerate all figures to get SVGs
    print("Regenerating figures with SVG output...\n")

    scripts_dir = Path(__file__).parent

    # Generate main figures
    for fig_num in range(1, 6):
        script = scripts_dir / f"generate_figure{fig_num}.py"
        if script.exists():
            print(f"Running {script.name}...")
            subprocess.run(["python", str(script)], cwd=str(scripts_dir))

    # Generate supplementary figures
    for supp_id in SUPP_FIGURES.keys():
        script = scripts_dir / f"generate_supp_figure{supp_id.lower()}.py"
        if script.exists():
            print(f"Running {script.name}...")
            subprocess.run(["python", str(script)], cwd=str(scripts_dir))

    print("\n" + "="*50)
    print("Organizing files for Illustrator...")
    print("="*50 + "\n")

    # Copy main figure files to organized structure
    for fig_num, base_name in FIGURES.items():
        fig_dir = OUTPUT_DIR / f"Figure{fig_num}"

        # Copy PNG
        png_src = SOURCE_DIR / f"{base_name}.png"
        if png_src.exists():
            png_dst = fig_dir / f"{base_name}.png"
            shutil.copy2(png_src, png_dst)
            print(f"  Figure {fig_num}: {png_dst.name}")

        # Copy PDF (vector, can be opened in Illustrator)
        pdf_src = SOURCE_DIR / f"{base_name}.pdf"
        if pdf_src.exists():
            pdf_dst = fig_dir / f"{base_name}.pdf"
            shutil.copy2(pdf_src, pdf_dst)
            print(f"  Figure {fig_num}: {pdf_dst.name}")

        # Copy SVG if exists
        svg_src = SOURCE_DIR / f"{base_name}.svg"
        if svg_src.exists():
            svg_dst = fig_dir / f"{base_name}.svg"
            shutil.copy2(svg_src, svg_dst)
            print(f"  Figure {fig_num}: {svg_dst.name}")

    # Copy supplementary figure files
    for supp_id, base_name in SUPP_FIGURES.items():
        fig_dir = OUTPUT_DIR / f"Supp{supp_id}"

        # Copy PNG
        png_src = SOURCE_DIR / f"{base_name}.png"
        if png_src.exists():
            png_dst = fig_dir / f"{base_name}.png"
            shutil.copy2(png_src, png_dst)
            print(f"  Supp {supp_id}: {png_dst.name}")

        # Copy PDF
        pdf_src = SOURCE_DIR / f"{base_name}.pdf"
        if pdf_src.exists():
            pdf_dst = fig_dir / f"{base_name}.pdf"
            shutil.copy2(pdf_src, pdf_dst)
            print(f"  Supp {supp_id}: {pdf_dst.name}")

        # Copy SVG
        svg_src = SOURCE_DIR / f"{base_name}.svg"
        if svg_src.exists():
            svg_dst = fig_dir / f"{base_name}.svg"
            shutil.copy2(svg_src, svg_dst)
            print(f"  Supp {supp_id}: {svg_dst.name}")

    print(f"\nFiles organized in: {OUTPUT_DIR}")
    print("\nFolder structure:")

    # List main figures
    for fig_num in range(1, 6):
        fig_dir = OUTPUT_DIR / f"Figure{fig_num}"
        files = list(fig_dir.glob("*"))
        print(f"  Figure{fig_num}/")
        for f in files:
            print(f"    - {f.name}")

    # List supplementary figures
    for supp_id in SUPP_FIGURES.keys():
        fig_dir = OUTPUT_DIR / f"Supp{supp_id}"
        if fig_dir.exists():
            files = list(fig_dir.glob("*"))
            print(f"  Supp{supp_id}/")
            for f in files:
                print(f"    - {f.name}")

if __name__ == "__main__":
    main()
