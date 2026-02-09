# Reproducibility Hub (v5)

This directory is the canonical manuscript reproducibility entrypoint for the **v5 submission set**.

## Scope

- Main figures: Figure 1-5
- Active supplementary figure: Supplementary Figure S2
- Vignettes: 1-6 (canonical wrappers under `python -m repro.cli.run_vignette`)
  - `5` is an extension of `2` for Figure 4 analysis packaging.
  - `6` demonstrates Figure 5 full-pipeline response analysis packaging.

## Quickstart

1. Configure paths via env vars or `repro/config/example.paths.yaml`.
2. Validate setup:

```bash
python -m repro.cli.validate_env --config repro/config/example.paths.yaml
```

3. Run canonical figure package:

```bash
python -m repro.cli.run_figures --set v5_all --config repro/config/example.paths.yaml
```

4. Run a vignette wrapper:

```bash
python -m repro.cli.run_vignette --id 2 --config repro/config/example.paths.yaml
```

## Orchestrator Targets

```bash
python -m repro.cli.repro --target quickstart --config repro/config/example.paths.yaml
python -m repro.cli.repro --target v5_figures --config repro/config/example.paths.yaml
python -m repro.cli.repro --target v5_full --config repro/config/example.paths.yaml
```

## Notes

- The compatibility phase is active: legacy scripts remain available while this hub becomes canonical.
- Legacy mapping reference: `repro/runbooks/legacy_mapping.md`.
- No canonical script should depend on hardcoded institutional absolute paths.
