#!/bin/bash
#SBATCH --job-name=gen_figures
#SBATCH --output=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/generate_all_%j.out
#SBATCH --error=/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/slurm_log/generate_all_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --cluster=htc
#SBATCH --partition=htc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alc376@pitt.edu

# Generate all manuscript figures

set -e

REPO="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist"
FIG_DIR="${REPO}/manuscript/figures"

export PATH="/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env/bin:$PATH"
source activate /ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/envs/CITEgeist_env

cd "${REPO}"

export CITEGEIST_DATA_ROOT="${REPO}"
export CITEGEIST_OUTPUT_ROOT="${REPO}/output"
export CITEGEIST_LICENSE_FILE="${REPO}/LICENSE"

echo "=============================================="
echo "Generating Manuscript Figures"
echo "Start time: $(date)"
echo "=============================================="
echo "DEPRECATION NOTICE: canonical entrypoint is now"
echo "  python -m repro.cli.run_figures --set v5_all --config repro/config/example.paths.yaml"

python -m repro.cli.run_figures --set v5_all

echo ""
echo "--- Exporting individual panels ---"
python "${FIG_DIR}/export_panels.py"
echo "Panel export complete"

echo ""
echo "=============================================="
echo "All figures generated (v5: Figures 1-5)"
echo "End time: $(date)"
echo "=============================================="

echo ""
echo "--- Re-rendering SVG schematics to PNG ---"
python3 << 'PYEOF'
try:
    from svglib.svglib import svg2rlg
    from reportlab.graphics import renderPM
    from pathlib import Path
    
    schematic_dir = Path("/ix1/alee/LO_LAB/Personal/Alexander_Chang/alc376/CITEgeist/manuscript/figures/output/schematics")
    rendered_dir = schematic_dir / "rendered"
    rendered_dir.mkdir(exist_ok=True)
    
    for svg_file in schematic_dir.glob("*.svg"):
        png_file = rendered_dir / f"{svg_file.stem}.png"
        drawing = svg2rlg(str(svg_file))
        if drawing:
            renderPM.drawToFile(drawing, str(png_file), fmt="PNG", dpi=300)
            print(f"  Rendered: {svg_file.name} -> {png_file.name}")
except ImportError:
    print("  svglib not available, skipping SVG rendering")
except Exception as e:
    print(f"  SVG rendering error: {e}")
PYEOF
echo "SVG rendering complete"
