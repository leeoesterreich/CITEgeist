"""Shared helpers for reproducibility CLI commands."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any


REQUIRED_ENV_KEYS = ("CITEGEIST_DATA_ROOT", "CITEGEIST_OUTPUT_ROOT", "CITEGEIST_LICENSE_FILE")


def repo_root() -> Path:
    """Return repository root from this module location."""
    return Path(__file__).resolve().parents[2]


def load_yaml_or_json(path: Path) -> dict[str, Any]:
    """Load YAML/JSON config with a safe fallback when PyYAML is unavailable."""
    text = path.read_text()
    try:
        import yaml  # type: ignore

        data = yaml.safe_load(text)
        if not isinstance(data, dict):
            raise ValueError(f"Config at {path} must be a mapping")
        return data
    except ModuleNotFoundError:
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                "PyYAML is not installed and config is not valid JSON. "
                f"Install pyyaml or convert {path} to JSON."
            ) from exc
        if not isinstance(data, dict):
            raise ValueError(f"Config at {path} must be a mapping")
        return data


def load_paths_config(config_path: str | None) -> dict[str, str]:
    """
    Resolve canonical path config with precedence:
    CLI config file > environment variables.
    """
    out: dict[str, str] = {}
    if config_path:
        raw = load_yaml_or_json(Path(config_path))
        paths = raw.get("paths", raw)
        if not isinstance(paths, dict):
            raise ValueError("Config must define a top-level 'paths' mapping or be a mapping itself.")
        for key in REQUIRED_ENV_KEYS:
            value = paths.get(key)
            if value:
                out[key] = str(value)

    for key in REQUIRED_ENV_KEYS:
        out.setdefault(key, os.getenv(key, ""))
    out.setdefault("CITEGEIST_TMPDIR", os.getenv("CITEGEIST_TMPDIR", "/tmp"))
    return out


def validate_required_paths(paths_cfg: dict[str, str]) -> list[str]:
    """Return user-facing validation errors for missing/unavailable paths."""
    errors: list[str] = []
    for key in REQUIRED_ENV_KEYS:
        value = paths_cfg.get(key, "")
        if not value:
            errors.append(f"{key} is required but not set.")
            continue
        p = Path(value)
        if key == "CITEGEIST_LICENSE_FILE":
            if not p.is_file():
                errors.append(f"{key} does not point to an existing file: {value}")
        else:
            if not p.exists():
                errors.append(f"{key} does not exist: {value}")
    return errors


def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Attach shared CLI arguments."""
    parser.add_argument(
        "--config",
        default=None,
        help="Path to YAML/JSON configuration file (see repro/config/example.paths.yaml).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and checks without executing workflows.",
    )


def run_subprocess(cmd: list[str], dry_run: bool = False, cwd: Path | None = None) -> int:
    """Run a subprocess command with dry-run support."""
    printable = " ".join(cmd)
    print(f"[repro] {printable}")
    if dry_run:
        return 0

    import subprocess

    proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=False)
    return proc.returncode


def die(message: str, code: int = 2) -> int:
    """Print a user-facing fatal message and return an exit code."""
    print(f"[repro][error] {message}", file=sys.stderr)
    return code

