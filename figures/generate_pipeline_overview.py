"""
Generate CITEgeist pipeline overview figure as SVG.

High-level horizontal flow: Module 1 -> 2 -> 3 -> 4 -> 5
Suitable for a paper graphical abstract or methods figure.
"""

import svgwrite


def create_pipeline_figure(output_path="figures/pipeline_overview.svg"):
    # --- Layout constants ---
    box_w = 230
    box_h = 165
    box_rx = 12
    gap = 55           # horizontal gap between boxes
    y_top = 10         # top margin
    pad_x = 15         # horizontal margin

    n_modules = 5
    canvas_w = n_modules * box_w + (n_modules - 1) * gap + 2 * pad_x
    label_area_h = 40      # space above boxes for arrow labels
    canvas_h = label_area_h + box_h + 10  # label area + boxes + bottom margin

    dwg = svgwrite.Drawing(output_path, size=(canvas_w, canvas_h))

    # --- Style constants ---
    arrow_color = "#333333"
    text_light = "#FFFFFF"
    text_dark = "#333333"
    font = "Helvetica, Arial, sans-serif"

    # Font sizes — maximized
    num_size = 34
    title_size = 17
    desc_size = 14
    label_size = 13

    # Module definitions: (number, title, description, color)
    modules = [
        ("1", "Marker Interest\nDetection", "Filter antibody markers\nby spatial signal", "#5B7FA5"),
        ("2", "Profile\nAssembly", "Discover cell type\nmarker profiles", "#7BA58A"),
        ("3", "Deconvolution", "Estimate proportions &\ndeconvolve expression", "#C4965A"),
        ("4", "Program\nDiscovery", "NMF-based spatial gene\nprograms per cell type", "#B07A8F"),
        ("5", "Cross-Sample\nIntegration", "Align programs\nacross patients", "#8A7EB5"),
    ]

    # Arrow labels (between modules)
    arrow_labels = [
        "Interesting\nmarkers",
        "Cell type\nprofiles",
        "Deconvolved\nlayers",
        "Programs +\nrelationships",
    ]

    x_start = pad_x
    box_top = label_area_h  # boxes start below the label area

    # --- Draw each module box ---
    for i, (num, title, desc, color) in enumerate(modules):
        x = x_start + i * (box_w + gap)
        y = box_top
        cx = x + box_w / 2

        # Drop shadow
        dwg.add(dwg.rect(
            insert=(x + 3, y + 3), size=(box_w, box_h),
            rx=box_rx, ry=box_rx, fill="#000000", opacity=0.08,
        ))

        # Main box
        dwg.add(dwg.rect(
            insert=(x, y), size=(box_w, box_h),
            rx=box_rx, ry=box_rx, fill=color, stroke=color, stroke_width=1,
        ))

        # Module number
        dwg.add(dwg.text(
            f"Module {num}", insert=(cx, y + 38),
            text_anchor="middle", font_size=num_size, font_weight="bold",
            font_family=font, fill=text_light, opacity=0.9,
        ))

        # Separator line
        dwg.add(dwg.line(
            start=(x + 20, y + 50), end=(x + box_w - 20, y + 50),
            stroke=text_light, stroke_width=0.8, opacity=0.4,
        ))

        # Title — vertically centered in middle zone
        title_lines = title.split("\n")
        n_title = len(title_lines)
        # Middle zone: y+54 to y+110 (56px). Center the title block in it.
        title_line_h = 22
        title_block_h = n_title * title_line_h
        title_top = y + 54 + (56 - title_block_h) / 2 + 14  # +14 for baseline
        for j, line in enumerate(title_lines):
            dwg.add(dwg.text(
                line, insert=(cx, title_top + j * title_line_h),
                text_anchor="middle", font_size=title_size, font_weight="bold",
                font_family=font, fill=text_light,
            ))

        # Description — fixed zone: y+115 to y+155
        desc_lines = desc.split("\n")
        desc_line_h = 18
        n_desc = len(desc_lines)
        desc_block_h = n_desc * desc_line_h
        desc_top = y + 115 + (40 - desc_block_h) / 2 + 12  # +12 for baseline
        for j, line in enumerate(desc_lines):
            dwg.add(dwg.text(
                line, insert=(cx, desc_top + j * desc_line_h),
                text_anchor="middle", font_size=desc_size,
                font_family=font, fill=text_light, opacity=0.85,
            ))

    # --- Draw arrows between boxes ---
    for i in range(n_modules - 1):
        x1 = x_start + i * (box_w + gap) + box_w       # right edge of box
        x2 = x_start + (i + 1) * (box_w + gap)          # left edge of next box
        arrow_y = box_top + box_h / 2

        # Arrow shaft
        dwg.add(dwg.line(
            start=(x1 + 4, arrow_y), end=(x2 - 14, arrow_y),
            stroke=arrow_color, stroke_width=2,
        ))

        # Arrowhead
        hs = 9
        dwg.add(dwg.polygon(
            [(x2 - 4, arrow_y),
             (x2 - 4 - hs, arrow_y - hs / 2),
             (x2 - 4 - hs, arrow_y + hs / 2)],
            fill=arrow_color,
        ))

        # Arrow label — positioned above the boxes in the label area
        label_lines = arrow_labels[i].split("\n")
        label_cx = (x1 + x2) / 2
        label_line_h = 16
        n_label = len(label_lines)
        # Place labels in the label_area, bottom-aligned just above box_top
        label_bottom_baseline = box_top - 6
        label_first_baseline = label_bottom_baseline - (n_label - 1) * label_line_h
        for j, line in enumerate(label_lines):
            dwg.add(dwg.text(
                line, insert=(label_cx, label_first_baseline + j * label_line_h),
                text_anchor="middle", font_size=label_size, font_style="italic",
                font_family=font, fill=text_dark,
            ))

    dwg.save()
    print(f"SVG saved to: {output_path}")
    print(f"Canvas: {canvas_w} x {canvas_h}")
    return output_path


if __name__ == "__main__":
    create_pipeline_figure()
