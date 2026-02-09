"""Run canonical v5 figure generation sets."""

from __future__ import annotations

import argparse
from pathlib import Path

from repro.cli.common import add_common_args, die, load_paths_config, repo_root, run_subprocess, validate_required_paths
from repro.pipelines.figures_v5 import v5_main_figure_tasks, v5_supp_figure_tasks


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run manuscript figure scripts for canonical v5 sets.")
    add_common_args(parser)
    parser.add_argument(
        "--set",
        choices=("v5_main", "v5_supp", "v5_all"),
        default="v5_main",
        help="Figure set to run.",
    )
    return parser


def _tasks_for_set(name: str, root: Path):
    if name == "v5_main":
        return v5_main_figure_tasks(root)
    if name == "v5_supp":
        return v5_supp_figure_tasks(root)
    return v5_main_figure_tasks(root) + v5_supp_figure_tasks(root)


def main() -> int:
    args = build_parser().parse_args()
    root = repo_root()
    paths_cfg = load_paths_config(args.config)
    errors = validate_required_paths(paths_cfg)
    if errors:
        for err in errors:
            print(f"[repro][invalid] {err}")
        return die("Path validation failed for figure run.", code=1)

    tasks = _tasks_for_set(args.set, root)
    for task in tasks:
        rc = run_subprocess(["python", str(task.script)], dry_run=args.dry_run, cwd=root)
        if rc != 0:
            return die(f"Figure task failed: {task.name}", code=rc)

    missing = [str(out) for task in tasks for out in task.expected_outputs if not out.exists()]
    if missing and not args.dry_run:
        for path in missing:
            print(f"[repro][missing-output] {path}")
        return die("One or more expected figure outputs are missing.", code=1)

    print(f"[repro] Completed figure set: {args.set}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

