"""Smoke test to verify biotooler package can be imported."""


def test_import_biotooler():
    """Test that biotooler package can be imported successfully."""
    import biotooler

    assert biotooler.__version__ == "0.1.0"
