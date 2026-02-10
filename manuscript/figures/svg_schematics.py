#!/usr/bin/env python3
"""
SVG Schematic Generator for CITEgeist Manuscript Figures

Generates clean SVG schematics for diagram panels. These can be:
1. Combined with matplotlib data panels in Illustrator
2. Rendered to PNG for embedding in matplotlib figures

SVG is written as raw XML for zero dependencies.
"""

from pathlib import Path

# Output directory
OUTPUT_DIR = Path(__file__).parent / "output" / "schematics"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Color palette (matching figure_style.py)
COLORS = {
    'primary': '#3498db',
    'secondary': '#9b59b6',
    'accent1': '#e67e22',
    'accent2': '#27ae60',
    'highlight': '#e91e63',
    'neutral': '#7f8c8d',
    'background': '#f8f9fa',
    'border': '#dee2e6',
    'text': '#2c3e50',
    # Module colors
    'module1': '#3498db',
    'module2': '#9b59b6',
    'module3': '#27ae60',
    'module4': '#e67e22',
    'module5': '#e91e63',
    # Cell type colors
    'tcells': '#e74c3c',
    'macrophages': '#9b59b6',
    'fibroblasts': '#3498db',
    'epithelial': '#27ae60',
    'endothelial': '#f39c12',
    'bcells': '#1abc9c',
}


def svg_header(width, height, viewbox=None):
    """Generate SVG header."""
    vb = viewbox or f"0 0 {width} {height}"
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="{vb}">
  <defs>
    <style>
      .title {{ font-family: Arial, sans-serif; font-size: 14px; font-weight: bold; }}
      .subtitle {{ font-family: Arial, sans-serif; font-size: 11px; }}
      .label {{ font-family: Arial, sans-serif; font-size: 10px; }}
      .small {{ font-family: Arial, sans-serif; font-size: 9px; }}
      .italic {{ font-style: italic; }}
      .bold {{ font-weight: bold; }}
      .white {{ fill: white; }}
      .neutral {{ fill: {COLORS['neutral']}; }}
    </style>
  </defs>
'''


def svg_footer():
    return '</svg>\n'


def rounded_rect(x, y, w, h, rx, fill, stroke='#2c3e50', stroke_width=1.5, opacity=1):
    """Generate rounded rectangle."""
    return f'  <rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" opacity="{opacity}"/>\n'


def text(x, y, content, css_class='label', anchor='middle', fill=None):
    """Generate text element."""
    fill_attr = f' fill="{fill}"' if fill else ''
    return f'  <text x="{x}" y="{y}" class="{css_class}" text-anchor="{anchor}"{fill_attr}>{content}</text>\n'


def arrow(x1, y1, x2, y2, color='#2c3e50', width=2):
    """Generate arrow line with arrowhead."""
    # Calculate arrowhead
    import math
    angle = math.atan2(y2 - y1, x2 - x1)
    arrow_len = 8
    arrow_angle = math.pi / 6

    ax1 = x2 - arrow_len * math.cos(angle - arrow_angle)
    ay1 = y2 - arrow_len * math.sin(angle - arrow_angle)
    ax2 = x2 - arrow_len * math.cos(angle + arrow_angle)
    ay2 = y2 - arrow_len * math.sin(angle + arrow_angle)

    return f'''  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>
  <polygon points="{x2},{y2} {ax1},{ay1} {ax2},{ay2}" fill="{color}"/>
'''


def circle(cx, cy, r, fill, stroke='none', stroke_width=1):
    """Generate circle."""
    return f'  <circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}"/>\n'


# =============================================================================
# Figure 1 Schematics
# =============================================================================

def figure1_panel_a():
    """Figure 1A: CITEgeist Modular Pipeline."""
    width, height = 800, 280
    svg = svg_header(width, height)

    # Title
    svg += text(400, 25, "CITEgeist Modular Pipeline", "title")

    # Module boxes
    modules = [
        (1, "Marker", "Interest", "Kurtosis/GMM/Moran's I", COLORS['module1']),
        (2, "Profile", "Discovery", "Colocalization + Clustering", COLORS['module2']),
        (3, "Deconvolution", "", "Proportions + GEX", COLORS['module3']),
        (4, "Program", "Discovery", "NMF + Spatial Coherence", COLORS['module4']),
        (5, "Integration", "", "Cross-sample Alignment", COLORS['module5']),
    ]

    box_w, box_h = 110, 120
    start_x, spacing = 95, 135
    y = 70

    for i, (num, line1, line2, subtitle, color) in enumerate(modules):
        x = start_x + i * spacing

        # Box
        svg += rounded_rect(x, y, box_w, box_h, 8, color, stroke_width=2)

        # Number circle
        svg += circle(x + 20, y + 20, 15, 'white', color, 2)
        svg += text(x + 20, y + 25, str(num), 'label bold', fill=color)

        # Title
        if line2:
            svg += text(x + box_w/2, y + 55, line1, 'label bold white')
            svg += text(x + box_w/2, y + 70, line2, 'label bold white')
        else:
            svg += text(x + box_w/2, y + 62, line1, 'label bold white')

        # Subtitle
        svg += text(x + box_w/2, y + box_h - 15, subtitle, 'small italic white')

        # Arrow to next
        if i < 4:
            svg += arrow(x + box_w + 5, y + box_h/2, x + spacing - 5, y + box_h/2)

    # Input label
    svg += rounded_rect(10, 100, 70, 60, 5, COLORS['background'], COLORS['border'])
    svg += text(45, 122, "Spatial", "small")
    svg += text(45, 135, "Transcriptomics", "small")
    svg += text(45, 148, "+ CITE-seq", "small")

    # Output label
    svg += rounded_rect(720, 100, 70, 60, 5, COLORS['background'], COLORS['border'])
    svg += text(755, 122, "Spatial", "small")
    svg += text(755, 135, "Programs +", "small")
    svg += text(755, 148, "Relationships", "small")

    # Bottom annotation
    svg += text(400, 220, "Protein-anchored deconvolution with automatic profile discovery", "subtitle italic neutral")

    svg += svg_footer()
    return svg


def figure1_panel_b():
    """Figure 1B: Spatial Statistics Foundation."""
    width, height = 400, 320
    svg = svg_header(width, height)

    # Title
    svg += text(200, 25, "Spatial Statistics Foundation", "title")

    # Section 1: Moran's I
    svg += text(70, 55, "Moran's I", "label bold")
    # Formula (simplified)
    svg += text(70, 90, "I = Σᵢⱼ wᵢⱼ(xᵢ-x̄)(xⱼ-x̄)", "small")
    svg += text(70, 110, "      / Σᵢ(xᵢ-x̄)²", "small")
    svg += text(70, 145, "Spatial autocorrelation", "small italic neutral")
    svg += text(70, 160, "(marker clustering)", "small italic neutral")

    # Section 2: Colocalization
    svg += text(200, 55, "Colocalization", "label bold")
    # Venn diagram
    svg += circle(180, 100, 30, COLORS['primary'], stroke_width=0)
    svg += f'  <circle cx="180" cy="100" r="30" fill="{COLORS["primary"]}" opacity="0.5"/>\n'
    svg += f'  <circle cx="220" cy="100" r="30" fill="{COLORS["highlight"]}" opacity="0.5"/>\n'
    svg += text(170, 105, "M1", "small bold")
    svg += text(230, 105, "M2", "small bold")
    svg += text(200, 145, "Bivariate Moran's I", "small italic neutral")
    svg += text(200, 160, "Co-occurrence scoring", "small italic neutral")

    # Section 3: Spatial Graph
    svg += text(330, 55, "Spatial Graph", "label bold")
    # Hexagonal graph nodes
    nodes = [(310, 90), (350, 90), (330, 120), (310, 120), (350, 120), (330, 90)]
    for nx, ny in nodes:
        svg += circle(nx, ny, 10, COLORS['secondary'])
    # Edges
    edges = [((310,90),(350,90)), ((310,90),(330,120)), ((350,90),(330,120)),
             ((310,120),(330,120)), ((350,120),(330,120)), ((310,90),(310,120)), ((350,90),(350,120))]
    for (x1,y1), (x2,y2) in edges:
        svg += f'  <line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{COLORS["text"]}" stroke-width="1" opacity="0.5"/>\n'
    svg += text(330, 155, "k-NN / Radius", "small italic neutral")
    svg += text(330, 170, "Neighbor weighting", "small italic neutral")

    # Bottom annotation
    svg += rounded_rect(50, 220, 300, 50, 5, COLORS['background'], COLORS['border'])
    svg += text(200, 245, "All modules leverage spatial context through", "small italic neutral")
    svg += text(200, 262, "neighbor-weighted statistics", "small italic neutral")

    svg += svg_footer()
    return svg


def figure1_panel_c():
    """Figure 1C: Resolution-Agnostic Design."""
    width, height = 400, 280
    svg = svg_header(width, height)

    # Title
    svg += text(200, 25, "Resolution-Agnostic Design", "title")

    # Left: Spot Resolution
    svg += text(100, 55, "Spot Resolution", "label bold", fill=COLORS['primary'])
    svg += text(100, 72, "(Visium, ~10-30 cells/spot)", "small italic neutral")

    # Large spots with mixed cells
    import random
    random.seed(42)
    spot_positions = [(60, 130), (100, 145), (140, 125), (80, 170)]
    cell_colors = [COLORS['tcells'], COLORS['macrophages'], COLORS['fibroblasts'], COLORS['epithelial']]

    for sx, sy in spot_positions:
        svg += f'  <circle cx="{sx}" cy="{sy}" r="25" fill="{random.choice(cell_colors)}" opacity="0.3" stroke="{random.choice(cell_colors)}" stroke-width="1.5"/>\n'
        # Inner cells
        for _ in range(4):
            dx, dy = random.randint(-12, 12), random.randint(-12, 12)
            svg += circle(sx + dx, sy + dy, 5, random.choice(cell_colors))

    svg += text(100, 215, "Mixed cell types per spot", "small italic neutral")

    # Center arrow
    svg += f'  <line x1="170" y1="140" x2="230" y2="140" stroke="{COLORS["text"]}" stroke-width="2.5" marker-end="url(#arrowhead)"/>\n'
    svg += f'  <line x1="230" y1="140" x2="170" y2="140" stroke="{COLORS["text"]}" stroke-width="2.5" marker-end="url(#arrowhead)"/>\n'
    svg += text(200, 120, "Same", "label bold")
    svg += text(200, 135, "Algorithm", "label bold")

    # Right: Cell Resolution
    svg += text(300, 55, "Cell Resolution", "label bold", fill=COLORS['accent2'])
    svg += text(300, 72, "(Xenium, CosMx, MERFISH)", "small italic neutral")

    # Individual cells
    for _ in range(25):
        cx = random.randint(250, 350)
        cy = random.randint(100, 180)
        svg += circle(cx, cy, 6, random.choice(cell_colors), 'white', 0.5)

    svg += text(300, 215, "Direct cell type observation", "small italic neutral")

    # Legend
    legend_items = [("T cells", COLORS['tcells']), ("Macrophage", COLORS['macrophages']),
                    ("Fibroblast", COLORS['fibroblasts']), ("Epithelial", COLORS['epithelial']),
                    ("Endothelial", COLORS['endothelial'])]
    for i, (label, color) in enumerate(legend_items):
        x = 30 + i * 78
        svg += circle(x, 255, 6, color)
        svg += text(x + 12, 259, label, "small", "start")

    svg += svg_footer()
    return svg


# =============================================================================
# Figure 2 Schematics
# =============================================================================

def figure2_panel_a():
    """Figure 2A: Marker Interest Detection (Module 1)."""
    width, height = 480, 280
    svg = svg_header(width, height)

    # Title
    svg += text(width / 2, 25, "Module 1: Marker Interest Detection", "title")

    # Input: Raw antibody data
    svg += rounded_rect(20, 55, 80, 50, 5, COLORS['background'], COLORS['border'])
    svg += text(60, 75, "Raw", "small bold")
    svg += text(60, 90, "Antibodies", "small")

    # Arrow
    svg += arrow(105, 80, 135, 80)

    # Three parallel test boxes
    tests = [
        ("Kurtosis", "≥ 2.0 = peaked", COLORS['primary'], 60),
        ("GMM", "2-comp SNR", COLORS['secondary'], 115),
        ("Moran's I", "Spatial clustering", COLORS['accent2'], 170),
    ]
    for label, desc, color, y in tests:
        svg += rounded_rect(140, y, 90, 40, 5, color, 'black', 1.5)
        svg += text(185, y + 18, label, "label bold white")
        svg += text(185, y + 32, desc, "small white")

    # Arrows from tests
    for y in [80, 135, 190]:
        svg += arrow(235, y, 265, y)

    # Decision logic box
    svg += rounded_rect(270, 90, 85, 65, 5, '#fff2cc', '#d6b656', 2)
    svg += text(312, 110, "Decision", "label bold")
    svg += text(312, 125, "Kurtosis OR", "small")
    svg += text(312, 140, "Moran's I", "small")

    # Arrows to outputs
    svg += arrow(360, 110, 390, 90)
    svg += arrow(360, 135, 390, 155)

    # Output: Interesting/Boring
    svg += rounded_rect(395, 70, 70, 35, 5, COLORS['accent2'], 'black', 1.5)
    svg += text(430, 92, "Interesting", "small bold white")

    svg += rounded_rect(395, 135, 70, 35, 5, COLORS['neutral'], 'black', 1.5)
    svg += text(430, 157, "Boring", "small bold white")

    # Bottom annotation
    svg += text(width / 2, 250, "Automatic identification of spatially-variable protein markers", "small italic neutral")

    svg += svg_footer()
    return svg


def figure2_panel_b():
    """Figure 2B: Profile Discovery Workflow (Module 2)."""
    width, height = 400, 280
    svg = svg_header(width, height)

    # Title
    svg += text(200, 25, "Module 2: Profile Discovery", "title")

    # Stage 1: Colocalization
    svg += rounded_rect(20, 55, 100, 55, 5, COLORS['primary'], '#2980b9', 1.5)
    svg += text(70, 75, "2a: Colocalization", "small bold white")
    svg += text(70, 92, "Same-spot, Adjacent", "small white")

    # Arrow
    svg += arrow(125, 82, 155, 82)

    # Stage 2: Graph
    svg += rounded_rect(160, 55, 100, 55, 5, COLORS['secondary'], '#7d3c98', 1.5)
    svg += text(210, 75, "2b: Profile Graph", "small bold white")
    svg += text(210, 92, "Clustering", "small white")

    # Arrow
    svg += arrow(265, 82, 295, 82)

    # Stage 3: Selection
    svg += rounded_rect(300, 55, 95, 55, 5, COLORS['accent2'], '#1e8449', 1.5)
    svg += text(347, 75, "2c: Selection", "small bold white")
    svg += text(347, 92, "Reconstruction", "small white")

    # Show example profiles below
    svg += text(200, 140, "Discovered Profiles:", "label bold")

    # Profile boxes with markers
    profiles = [
        ("Epithelial", ["PanCK", "EPCAM"], COLORS['epithelial']),
        ("Macrophage", ["CD68", "CD163"], COLORS['macrophages']),
        ("T cells", ["CD3E", "CD4/8"], COLORS['tcells']),
        ("Fibroblast", ["aSMA", "VIM"], COLORS['fibroblasts']),
    ]
    for i, (name, markers, color) in enumerate(profiles):
        x = 30 + i * 95
        svg += rounded_rect(x, 160, 85, 60, 5, color, 'black', 1.5)
        svg += text(x + 42, 178, name, "small bold white")
        svg += text(x + 42, 195, markers[0], "small white")
        svg += text(x + 42, 210, markers[1], "small white")

    # Bottom annotation
    svg += text(200, 250, "Automatic marker-to-cell-type profile assembly", "small italic neutral")

    svg += svg_footer()
    return svg


def figure2_panel_c():
    """Figure 2C: Xenium RCC Single-Cell Demonstration."""
    width, height = 400, 280
    svg = svg_header(width, height)

    # Title
    svg += text(200, 25, "Xenium RCC Demonstration", "title")
    svg += text(200, 45, "Single-Cell Resolution Validation", "subtitle neutral")

    # Left: Ground truth (protein)
    svg += text(100, 75, "Protein Ground Truth", "label bold", fill=COLORS['primary'])

    # Simulated tissue section with colored cells
    import random
    random.seed(123)
    cell_colors = [COLORS['epithelial'], COLORS['macrophages'], COLORS['tcells'],
                   COLORS['fibroblasts'], COLORS['endothelial']]
    for _ in range(35):
        cx = random.randint(45, 155)
        cy = random.randint(90, 180)
        svg += circle(cx, cy, 6, random.choice(cell_colors), 'white', 0.5)

    # Right: CITEgeist output
    svg += text(300, 75, "CITEgeist Output", "label bold", fill=COLORS['accent2'])

    # Similar pattern showing matching
    random.seed(123)  # Same seed for matching pattern
    for _ in range(35):
        cx = random.randint(245, 355)
        cy = random.randint(90, 180)
        svg += circle(cx, cy, 6, random.choice(cell_colors), 'white', 0.5)

    # Connecting arrows
    svg += f'  <line x1="165" y1="135" x2="235" y2="135" stroke="{COLORS["text"]}" stroke-width="2" stroke-dasharray="4,2"/>\n'
    svg += text(200, 125, "vs", "label bold")

    # Accuracy metric
    svg += rounded_rect(150, 200, 100, 40, 5, '#d5e8d4', '#82b366', 2)
    svg += text(200, 218, "Acc: 94.2%", "label bold")
    svg += text(200, 233, "5 regions", "small neutral")

    # Legend
    legend = [("Epi", COLORS['epithelial']), ("Mac", COLORS['macrophages']),
              ("T", COLORS['tcells']), ("Fib", COLORS['fibroblasts']), ("Endo", COLORS['endothelial'])]
    for i, (label, color) in enumerate(legend):
        x = 40 + i * 70
        svg += circle(x, 265, 6, color)
        svg += text(x + 12, 269, label, "small", "start")

    svg += svg_footer()
    return svg


# =============================================================================
# Figure 3 Schematic
# =============================================================================

def figure3_panel_a():
    """Figure 3A: Two-Pass Deconvolution."""
    width, height = 500, 200
    svg = svg_header(width, height)

    # Title
    svg += text(250, 25, "Two-Pass Deconvolution", "title")

    # Input: Mixed spot (pie chart style)
    svg += text(50, 55, "Mixed", "small")
    svg += text(50, 68, "Spot", "small")
    colors = [COLORS['tcells'], COLORS['macrophages'], COLORS['fibroblasts'], COLORS['epithelial']]
    # Simplified pie as colored segments
    svg += f'''  <path d="M50,110 L50,85 A25,25 0 0,1 75,110 Z" fill="{colors[0]}"/>
  <path d="M50,110 L75,110 A25,25 0 0,1 50,135 Z" fill="{colors[1]}"/>
  <path d="M50,110 L50,135 A25,25 0 0,1 25,110 Z" fill="{colors[2]}"/>
  <path d="M50,110 L25,110 A25,25 0 0,1 50,85 Z" fill="{colors[3]}"/>
'''

    # Arrow
    svg += arrow(85, 110, 120, 110)

    # Pass 1 box
    svg += rounded_rect(125, 75, 130, 70, 8, '#d5e8d4', '#82b366', 2)
    svg += text(190, 100, "Pass 1", "label bold")
    svg += text(190, 118, "Protein → Proportions", "small")
    svg += text(190, 135, "QP + Spatial smoothing", "small italic neutral")

    # Arrow
    svg += arrow(260, 110, 295, 110)

    # Pass 2 box
    svg += rounded_rect(300, 75, 130, 70, 8, '#dae8fc', '#6c8ebf', 2)
    svg += text(365, 100, "Pass 2", "label bold")
    svg += text(365, 118, "RNA → GEX Layers", "small")
    svg += text(365, 135, "Per-cell-type expression", "small italic neutral")

    # Arrow
    svg += arrow(435, 110, 470, 110)

    # Outputs
    svg += text(480, 60, "Outputs", "label bold", "start")
    # Proportion bar
    svg += rounded_rect(470, 75, 25, 15, 2, colors[0])
    svg += rounded_rect(495, 75, 30, 15, 2, colors[1])
    svg += rounded_rect(525, 75, 25, 15, 2, colors[2])
    svg += rounded_rect(550, 75, 20, 15, 2, colors[3])
    svg += text(520, 105, "Proportions", "small", "start")

    # GEX layers
    for i, c in enumerate(colors[:3]):
        svg += rounded_rect(475 + i*25, 115, 20, 35, 2, c, 'white', 0.5)
    svg += text(520, 165, "GEX Layers", "small", "start")

    # Bottom note
    svg += rounded_rect(100, 175, 300, 25, 5, COLORS['background'], COLORS['border'])
    svg += text(250, 192, "Same-slide protein anchors deconvolution", "small italic neutral")

    svg += svg_footer()
    return svg


# =============================================================================
# Figure 4 Schematic
# =============================================================================

def figure4_panel_a():
    """Figure 4A: NMF Program Discovery."""
    width, height = 400, 200
    svg = svg_header(width, height)

    # Title
    svg += text(200, 25, "Program Discovery", "title")

    # Expression matrix
    svg += rounded_rect(20, 50, 80, 100, 5, '#dae8fc', '#6c8ebf', 2)
    svg += text(60, 95, "Expression", "label bold")
    svg += text(60, 112, "Matrix", "label bold")
    svg += text(60, 160, "Spots × Genes", "small neutral")

    # Arrow
    svg += arrow(105, 100, 135, 100)

    # NMF box
    svg += rounded_rect(140, 65, 70, 70, 5, '#fff2cc', '#d6b656', 2)
    svg += text(175, 105, "NMF", "label bold")

    # Arrow
    svg += arrow(215, 100, 245, 100)

    # W matrix
    svg += rounded_rect(250, 45, 60, 70, 5, '#d5e8d4', '#82b366', 2)
    svg += text(280, 85, "W", "title bold")
    svg += text(280, 130, "Gene", "small neutral")
    svg += text(280, 145, "Loadings", "small neutral")

    # H matrix
    svg += rounded_rect(250, 125, 90, 35, 5, '#d5e8d4', '#82b366', 2)
    svg += text(295, 147, "H", "label bold")
    svg += text(350, 140, "Program", "small neutral", "start")
    svg += text(350, 155, "Activities", "small neutral", "start")

    # Arrow down to validation
    svg += arrow(295, 165, 295, 190)

    # Moran's I validation
    svg += rounded_rect(250, 195, 90, 35, 5, '#f8cecc', '#b85450', 2)
    svg += text(295, 210, "Moran's I", "label bold")
    svg += text(295, 225, "Validation", "small")

    svg += svg_footer()
    return svg


# =============================================================================
# Figure 5 Schematic
# =============================================================================

def figure5_panel_a():
    """Figure 5A: Cross-Sample Integration."""
    width, height = 400, 200
    svg = svg_header(width, height)

    # Title
    svg += text(200, 25, "Cross-Sample Integration", "title")

    # Sample boxes (left)
    sample_colors = ['#2ecc71', '#2ecc71', '#e74c3c', '#e74c3c']  # responder/progressor
    for i, c in enumerate(sample_colors):
        svg += rounded_rect(20, 50 + i*35, 50, 30, 3, c, 'black', 1)
        svg += text(45, 70 + i*35, f"S{i+1}", "label bold")
    svg += text(45, 45, "12 Samples", "small bold")

    # Arrow
    svg += arrow(75, 110, 105, 110)

    # Module 4 box
    svg += rounded_rect(110, 75, 70, 70, 5, '#dae8fc', '#6c8ebf', 2)
    svg += text(145, 105, "Module 4", "small bold")
    svg += text(145, 120, "Programs", "small")

    # Arrow
    svg += arrow(185, 110, 215, 110)

    # Harmony box
    svg += rounded_rect(220, 75, 70, 70, 5, '#fff2cc', '#d6b656', 2)
    svg += text(255, 105, "Harmony", "small bold")
    svg += text(255, 120, "Alignment", "small")

    # Arrow
    svg += arrow(295, 110, 325, 110)

    # Outputs
    svg += rounded_rect(330, 55, 90, 45, 5, '#d5e8d4', 'black', 1.5)
    svg += text(375, 75, "73 Aligned", "small bold")
    svg += text(375, 90, "Programs", "small")

    svg += rounded_rect(330, 110, 90, 45, 5, '#f8cecc', 'black', 1.5)
    svg += text(375, 130, "Response", "small bold")
    svg += text(375, 145, "Analysis", "small")

    # Legend
    svg += rounded_rect(30, 190, 12, 12, 2, '#2ecc71')
    svg += text(50, 200, "Responder", "small", "start")
    svg += rounded_rect(130, 190, 12, 12, 2, '#e74c3c')
    svg += text(150, 200, "Progressor", "small", "start")

    svg += svg_footer()
    return svg


# =============================================================================
# Figure 6 Schematics
# =============================================================================

def figure6_panel_a():
    """Figure 6A: Downstream Analysis Workflow."""
    width, height = 350, 250
    svg = svg_header(width, height)

    # Title
    svg += text(175, 20, "CITEgeist → Downstream Analysis", "label bold")

    # CITEgeist outputs (left column)
    outputs = [("Proportions", 50), ("GEX Layers", 100), ("Programs", 150)]
    for label, y in outputs:
        svg += rounded_rect(15, y, 80, 35, 5, COLORS['primary'], '#2980b9', 1.5)
        svg += text(55, y + 22, label, "small bold white")

    # Arrows
    for y in [67, 117, 167]:
        svg += arrow(100, y, 130, y)

    # Tools (middle)
    tools = [
        ("scanpy", 45, COLORS['accent2'], "Clustering"),
        ("PyDESeq2", 90, COLORS['accent1'], "Diff. Expression"),
        ("GSEApy", 135, COLORS['highlight'], "Pathway Enrichment"),
        ("COMMOT", 180, '#e74c3c', "Cell Communication"),
    ]
    for label, y, color, desc in tools:
        svg += rounded_rect(135, y, 75, 30, 5, color, 'black', 1)
        svg += text(172, y + 20, label, "small bold white")
        svg += text(220, y + 20, desc, "small neutral", "start")

    # Arrows to results
    for y in [60, 105, 150, 195]:
        svg += f'  <line x1="215" y1="{y}" x2="280" y2="{y}" stroke="{COLORS["neutral"]}" stroke-width="1" stroke-dasharray="4,2"/>\n'

    # Results box
    svg += rounded_rect(285, 100, 55, 50, 5, COLORS['background'], COLORS['border'])
    svg += text(312, 120, "Analysis", "small bold")
    svg += text(312, 135, "Results", "small bold")

    svg += svg_footer()
    return svg


def figure6_panel_e():
    """Figure 6E: Experimental Validation Summary."""
    width, height = 500, 180
    svg = svg_header(width, height)

    # Title
    svg += text(250, 20, "Experimental Validation: Midkine Discovery", "label bold")

    # Validation points
    points = [
        ("CITEgeist Discovery:", "MDK identified as spatially variable in ER+ breast cancer"),
        ("ChIP-seq Validation:", "ER binding at MDK locus confirmed, E2-dependent regulation"),
        ("ELISA Confirmation:", "Secreted MDK levels correlate with ER signaling"),
        ("Mechanism:", "MDK expression modulated by ER-FOXA1 axis"),
    ]

    for i, (title, desc) in enumerate(points):
        y = 50 + i * 30
        svg += circle(25, y, 5, COLORS['primary'])
        svg += text(35, y + 4, title, "small bold", "start")
        svg += text(160, y + 4, desc, "small neutral", "start")

    # Bottom insight box
    svg += rounded_rect(50, 170, 400, 30, 5, COLORS['background'], COLORS['border'])
    svg += text(250, 190, "CITEgeist spatial analysis → testable hypothesis → experimental confirmation", "small italic neutral")

    svg += svg_footer()
    return svg


# =============================================================================
# Generate all schematics
# =============================================================================

def generate_all_schematics():
    """Generate all SVG schematic files."""
    schematics = {
        'figure1_panel_a_pipeline.svg': figure1_panel_a,
        'figure1_panel_b_spatial_stats.svg': figure1_panel_b,
        'figure1_panel_c_resolution.svg': figure1_panel_c,
        'figure2_panel_a_marker_interest.svg': figure2_panel_a,
        'figure2_panel_b_profile_discovery.svg': figure2_panel_b,
        'figure2_panel_c_xenium_demo.svg': figure2_panel_c,
        'figure3_panel_a_deconvolution.svg': figure3_panel_a,
        'figure4_panel_a_nmf.svg': figure4_panel_a,
        'figure5_panel_a_integration.svg': figure5_panel_a,
        'figure6_panel_a_workflow.svg': figure6_panel_a,
        'figure6_panel_e_validation.svg': figure6_panel_e,
    }

    for filename, func in schematics.items():
        svg_content = func()
        filepath = OUTPUT_DIR / filename
        with open(filepath, 'w') as f:
            f.write(svg_content)
        print(f"Generated: {filepath}")

    print(f"\nAll schematics saved to: {OUTPUT_DIR}")
    print("\nThese SVGs can be:")
    print("  1. Opened directly in Illustrator/Inkscape for editing")
    print("  2. Combined with matplotlib data panels")
    print("  3. Converted to PNG using Inkscape CLI or online tools")


if __name__ == "__main__":
    generate_all_schematics()
