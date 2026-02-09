from repro.cli.common import validate_required_paths


def test_missing_paths_report_errors() -> None:
    cfg = {
        "CITEGEIST_DATA_ROOT": "",
        "CITEGEIST_OUTPUT_ROOT": "",
        "CITEGEIST_LICENSE_FILE": "",
    }
    errors = validate_required_paths(cfg)
    assert len(errors) == 3

