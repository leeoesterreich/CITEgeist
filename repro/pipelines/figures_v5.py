"""Figure-set definitions for v5 manuscript reproducibility."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FigureTask:
    """Single figure generation task."""

    name: str
    script: Path
    expected_outputs: tuple[Path, ...]


def v5_main_figure_tasks(repo_root: Path) -> list[FigureTask]:
    """Figure 1-5 tasks in v5 manuscript set."""
    out = repo_root / "manuscript" / "figures" / "output"
    scripts = repo_root / "manuscript" / "figures"
    return [
        FigureTask("figure1", scripts / "generate_figure1.py", (out / "figure1_pipeline_overview.pdf",)),
        FigureTask("figure2", scripts / "generate_figure2.py", (out / "figure2_profile_discovery.pdf",)),
        FigureTask("figure3", scripts / "generate_figure3.py", (out / "figure3_benchmarking.pdf",)),
        FigureTask("figure4", scripts / "generate_figure4.py", (out / "figure4_midkine_esr1.pdf",)),
        FigureTask("figure5", scripts / "generate_figure5.py", (out / "figure5_full_pipeline.pdf",)),
    ]


def v5_supp_figure_tasks(repo_root: Path) -> list[FigureTask]:
    """Active supplementary figures in current v5 flow."""
    out = repo_root / "manuscript" / "figures" / "output"
    scripts = repo_root / "manuscript" / "figures"
    return [
        FigureTask(
            "supp_figure2",
            scripts / "generate_supp_figure2.py",
            (out / "supp_figure2_de_pathway.pdf",),
        )
    ]

