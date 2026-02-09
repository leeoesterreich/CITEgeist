# Legacy-to-Canonical Mapping (Compatibility Phase)

Use this mapping during the transition window.

- `manuscript/figures/sbatch_generate_all.sh`
  canonical: `python -m repro.cli.run_figures --set v5_all --config repro/config/example.paths.yaml`
- `manuscript/figures/sbatch_generate_supp_figures.sh`
  canonical: `python -m repro.cli.run_figures --set v5_supp --config repro/config/example.paths.yaml`
- `examples/sbatch_sample.sh`
  canonical: `python -m repro.cli.repro --target v5_full --config repro/config/example.paths.yaml`

Legacy scripts remain callable, but `repro/` is the authoritative interface.

