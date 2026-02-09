"""Top-level orchestrator for reproducibility targets."""

from __future__ import annotations

import argparse

from repro.cli.common import add_common_args, die, load_paths_config, repo_root, run_subprocess, validate_required_paths
from repro.pipelines.v5_submission import build_targets


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Canonical reproducibility orchestrator for manuscript v5.")
    add_common_args(parser)
    parser.add_argument(
        "--target",
        choices=("quickstart", "v5_figures", "v5_full"),
        required=True,
        help="Reproducibility target.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = repo_root()
    paths_cfg = load_paths_config(args.config)
    errors = validate_required_paths(paths_cfg)
    if errors:
        for err in errors:
            print(f"[repro][invalid] {err}")
        return die("Path validation failed.", code=1)

    targets = build_targets(root)
    target = targets[args.target]
    print(f"[repro] Target: {target.name}")
    print(f"[repro] Description: {target.description}")
    if not target.commands:
        print("[repro] No execution commands for this target (metadata/validation target).")
        return 0

    for cmd in target.commands:
        rc = run_subprocess(cmd, dry_run=args.dry_run, cwd=root)
        if rc != 0:
            return die(f"Target '{target.name}' failed.", code=rc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

