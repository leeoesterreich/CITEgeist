"""Module rerun command templates for v5 reproducibility."""

from __future__ import annotations

from pathlib import Path


def module_commands(repo_root: Path) -> list[list[str]]:
    """
    Canonical module rerun templates.

    These commands are intentionally explicit and parameterized by env/config
    instead of hardcoded institutional paths.
    """
    examples = repo_root / "examples"
    return [
        ["python", str(examples / "run_module12_discovery.py"), "--sample", "HCC22-088-P1-S1"],
        ["python", str(examples / "run_module3_unified.py"), "--sample", "HCC22-088-P1-S1"],
        ["python", str(examples / "run_module4_discovery.py"), "--sample", "HCC22-088-P1-S1"],
        ["python", str(examples / "run_module4b_bivariate.py"), "--sample", "HCC22-088-P1-S1"],
        ["python", str(examples / "run_module5_integration.py")],
    ]

