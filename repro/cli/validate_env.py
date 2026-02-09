"""Validate canonical reproducibility environment configuration."""

from __future__ import annotations

import argparse
from pathlib import Path

from repro.cli.common import add_common_args, die, load_paths_config, load_yaml_or_json, repo_root, validate_required_paths


def validate_manifest(manifest_path: Path) -> list[str]:
    """Lightweight manifest schema validation."""
    errors: list[str] = []
    if not manifest_path.exists():
        return [f"Manifest not found: {manifest_path}"]

    manifest = load_yaml_or_json(manifest_path)
    required_top = ("version", "submission_scope", "datasets", "targets", "figures", "vignettes")
    for key in required_top:
        if key not in manifest:
            errors.append(f"Manifest missing key: {key}")
    if manifest.get("submission_scope") != "v5":
        errors.append("submission_scope must be 'v5'")
    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate CITEgeist reproducibility environment and manifest.")
    add_common_args(parser)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    paths_cfg = load_paths_config(args.config)
    errors = validate_required_paths(paths_cfg)
    errors.extend(validate_manifest(repo_root() / "repro" / "MANIFEST.yaml"))

    if errors:
        for err in errors:
            print(f"[repro][invalid] {err}")
        return die("Environment validation failed.", code=1)

    print("[repro] Environment and manifest validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

