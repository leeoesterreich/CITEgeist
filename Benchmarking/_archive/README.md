# Archived Benchmark Results

This directory contains benchmark results from before the cleanup on 2026-02-20.

## Why Archived

These results were generated during exploratory development with inconsistent:
- Ground truth datasets (RNA-based vs protein-gated)
- Evaluation metrics
- Method configurations

## Canonical Results

The canonical benchmark results are now in:
- `Benchmarking/canonical_results.json`

## Contents

- `2026-02-20_pre_cleanup/` - Snapshot before cleanup (not committed, for local reference only)

## Regenerating Results

To regenerate canonical results:
```bash
python Benchmarking/scripts/consolidate_results.py
```
