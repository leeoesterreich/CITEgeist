# Reviewer Quickstart (v5)

## 1. Configure environment

Set required values:

- `CITEGEIST_DATA_ROOT`
- `CITEGEIST_OUTPUT_ROOT`
- `CITEGEIST_LICENSE_FILE`

Or edit `repro/config/example.paths.yaml`.

## 2. Validate

```bash
python -m repro.cli.validate_env --config repro/config/example.paths.yaml
```

## 3. Reproduce manuscript figures

```bash
python -m repro.cli.run_figures --set v5_all --config repro/config/example.paths.yaml
```

## 4. Inspect a vignette

```bash
python -m repro.cli.run_vignette --id 1 --config repro/config/example.paths.yaml
python -m repro.cli.run_vignette --id 5 --config repro/config/example.paths.yaml
python -m repro.cli.run_vignette --id 6 --config repro/config/example.paths.yaml
```
