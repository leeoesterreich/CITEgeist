from pathlib import Path

from repro.cli.validate_env import validate_manifest


def test_manifest_schema_minimal() -> None:
    path = Path("repro/MANIFEST.yaml")
    errors = validate_manifest(path)
    assert errors == []

