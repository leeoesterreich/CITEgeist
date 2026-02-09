import subprocess


def test_validate_help_smoke() -> None:
    proc = subprocess.run(["python", "-m", "repro.cli.validate_env", "--help"], check=False, capture_output=True, text=True)
    assert proc.returncode == 0
    assert "Validate CITEgeist reproducibility environment" in proc.stdout

