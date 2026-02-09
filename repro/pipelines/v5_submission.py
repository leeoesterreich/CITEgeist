"""Target graph for v5 manuscript reproducibility."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from repro.pipelines.figures_v5 import v5_main_figure_tasks, v5_supp_figure_tasks
from repro.pipelines.modules_v5 import module_commands


@dataclass(frozen=True)
class Target:
    """Executable reproducibility target."""

    name: str
    description: str
    commands: tuple[list[str], ...]


def build_targets(repo_root: Path) -> dict[str, Target]:
    """Build canonical target map."""
    figure_commands = [["python", str(task.script)] for task in v5_main_figure_tasks(repo_root)]
    supp_commands = [["python", str(task.script)] for task in v5_supp_figure_tasks(repo_root)]

    targets = {
        "quickstart": Target(
            name="quickstart",
            description="Validate environment and print canonical v5 commands.",
            commands=tuple(),
        ),
        "v5_figures": Target(
            name="v5_figures",
            description="Generate v5 main and active supplementary figures.",
            commands=tuple(figure_commands + supp_commands),
        ),
        "v5_full": Target(
            name="v5_full",
            description="Run modules and figures required for v5 submission package.",
            commands=tuple(module_commands(repo_root) + figure_commands + supp_commands),
        ),
    }
    return targets

