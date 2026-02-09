# Reproducibility Guide (Canonical v5)

The canonical entrypoint for manuscript reproducibility is now:

- `repro/README.md`
- `repro/runbooks/reviewer_quickstart.md`
- CLI: `python -m repro.cli.*`

## For Reviewers

### 1. Download the code and data

1. Download the code from Figshare: [https://figshare.com/s/34e456fd7786e5211acc](https://figshare.com/s/34e456fd7786e5211acc)
2. Download data from GEO (see manuscript Data Availability).

### 2. Prepare data

Use the publication prep script (correct filename):

```bash
python publication_prep_scripts/delete_all_but_essential.py --folder GSE289326/processed_files/
```

### 3. Configure paths and validate

```bash
python -m repro.cli.validate_env --config repro/config/example.paths.yaml
```

Required path keys:

- `CITEGEIST_DATA_ROOT`
- `CITEGEIST_OUTPUT_ROOT`
- `CITEGEIST_LICENSE_FILE`

### 4. Generate manuscript figure package (v5)

```bash
python -m repro.cli.run_figures --set v5_all --config repro/config/example.paths.yaml
```

### 5. Run vignettes

```bash
python -m repro.cli.run_vignette --id 1 --config repro/config/example.paths.yaml
python -m repro.cli.run_vignette --id 2 --config repro/config/example.paths.yaml
python -m repro.cli.run_vignette --id 3 --config repro/config/example.paths.yaml
python -m repro.cli.run_vignette --id 4 --config repro/config/example.paths.yaml
```

## Notes

- Scope is the v5 manuscript submission package.
- Legacy scripts remain available during the compatibility phase, but `repro/` is canonical.
