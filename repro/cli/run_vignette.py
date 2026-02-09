"""Run or inspect canonical vignette notebooks."""

from __future__ import annotations

import argparse
from pathlib import Path

from repro.cli.common import add_common_args, die, load_paths_config, repo_root, run_subprocess, validate_required_paths


VIGNETTE_PATHS = {
    "1": "repro/vignettes/vignette_1_biopsy_heterogeneity.ipynb",
    "2": "repro/vignettes/vignette_2_surgical_d538g.ipynb",
    "3": "repro/vignettes/vignette_3_responder_macrophages.ipynb",
    "4": "repro/vignettes/vignette_4_interoperability_commot.ipynb",
    "5": "repro/vignettes/vignette_5_figure4_midkine_esr1.ipynb",
    "6": "repro/vignettes/vignette_6_figure5_full_pipeline.ipynb",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Inspect or execute canonical vignettes.")
    add_common_args(parser)
    parser.add_argument("--id", choices=("1", "2", "3", "4", "5", "6"), required=True, help="Vignette identifier.")
    parser.add_argument("--execute", action="store_true", help="Execute notebook with jupyter nbconvert.")
    parser.add_argument("--out", default=None, help="Output executed notebook path (optional).")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = repo_root()
    paths_cfg = load_paths_config(args.config)
    errors = validate_required_paths(paths_cfg)
    if errors:
        for err in errors:
            print(f"[repro][invalid] {err}")
        return die("Path validation failed for vignette run.", code=1)

    notebook = root / VIGNETTE_PATHS[args.id]
    if not notebook.exists():
        return die(f"Notebook not found: {notebook}", code=1)

    print(f"[repro] Vignette {args.id}: {notebook}")
    if not args.execute:
        return 0

    output_path = Path(args.out) if args.out else notebook.with_name(f"{notebook.stem}.executed.ipynb")
    cmd = [
        "jupyter",
        "nbconvert",
        "--to",
        "notebook",
        "--execute",
        "--inplace" if args.out is None else "--output",
    ]
    if args.out is None:
        cmd.append(str(notebook))
    else:
        cmd.extend([str(output_path), str(notebook)])
    rc = run_subprocess(cmd, dry_run=args.dry_run, cwd=root)
    if rc != 0:
        return die("Notebook execution failed.", code=rc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
