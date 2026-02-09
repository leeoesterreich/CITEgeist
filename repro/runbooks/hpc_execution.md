# HPC Execution Notes

Use the same canonical `repro` CLI targets on HPC login/compute nodes.

## Recommended flow

1. Validate on login node:

```bash
python -m repro.cli.validate_env --config repro/config/example.paths.yaml
```

2. Submit cluster wrapper jobs that call canonical commands, e.g.:

```bash
python -m repro.cli.run_figures --set v5_all --config repro/config/example.paths.yaml
```

3. Keep all path configuration in env vars or the config file, not in script source.

