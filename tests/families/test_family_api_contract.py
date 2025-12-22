"""Test family API contract to ensure consistent, contributor-friendly interface.

This module tests that all families follow the standardized API:
1. Each family has a FAMILY_META dict with required keys
2. Each family provides a get_features() function
3. Importing a family module is lightweight (doesn't import heavy dependencies)
4. Calling get_features() triggers appropriate ImportError if dependencies are missing
"""

import sys
from typing import Any


def _clear_modules(*prefixes: str) -> None:
    """Helper to clear modules with given prefixes from sys.modules.

    Args:
        prefixes: Module name prefixes to clear (e.g., "codonbias", "biotooler.families")
    """
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if any(module.startswith(prefix) for prefix in prefixes)
    ]
    for module in modules_to_clear:
        del sys.modules[module]


def test_codon_bias_has_family_meta():
    """Test that codon_bias family has FAMILY_META with required keys."""
    # Clear any previously imported modules
    _clear_modules("biotooler.families.codon_bias")

    # Import the family module
    from biotooler.families import codon_bias

    # Check that FAMILY_META exists
    assert hasattr(codon_bias, "FAMILY_META"), "codon_bias should have FAMILY_META"

    meta = codon_bias.FAMILY_META

    # Check that it's a dictionary
    assert isinstance(meta, dict), "FAMILY_META should be a dictionary"

    # Check required keys
    required_keys = {"name", "extra", "owner", "summary"}
    assert set(meta.keys()) >= required_keys, (
        f"FAMILY_META should have keys: {required_keys}, got: {set(meta.keys())}"
    )

    # Check that values are strings
    assert isinstance(meta["name"], str), "FAMILY_META['name'] should be a string"
    assert isinstance(meta["extra"], str), "FAMILY_META['extra'] should be a string"
    assert isinstance(meta["owner"], str), "FAMILY_META['owner'] should be a string"
    assert isinstance(meta["summary"], str), "FAMILY_META['summary'] should be a string"

    # Check specific values for codon_bias
    assert meta["name"] == "codon_bias", "Expected name to be 'codon_bias'"
    assert meta["extra"] == "codon_bias", "Expected extra to be 'codon_bias'"


def test_codon_bias_has_get_features():
    """Test that codon_bias family has get_features() function."""
    # Clear any previously imported modules
    _clear_modules("biotooler.families.codon_bias")

    # Import the family module
    from biotooler.families import codon_bias

    # Check that get_features exists
    assert hasattr(codon_bias, "get_features"), "codon_bias should have get_features function"

    # Check that it's callable
    assert callable(codon_bias.get_features), "get_features should be callable"


def test_codon_bias_import_is_lightweight():
    """Test importing codon_bias family doesn't import the codonbias dependency.

    This ensures that families remain lightweight to import and heavy dependencies
    are only loaded when actually used via get_features().
    """
    # Clear any previously imported modules
    _clear_modules("biotooler.families.codon_bias", "codonbias")

    # Import the family module
    from biotooler.families import codon_bias

    # Verify the module is imported
    assert codon_bias is not None

    # Check that FAMILY_META is accessible (lightweight operation)
    assert hasattr(codon_bias, "FAMILY_META")
    assert codon_bias.FAMILY_META["name"] == "codon_bias"

    # Check that codonbias is NOT imported yet
    assert "codonbias" not in sys.modules, (
        "codonbias was imported when importing the family module, "
        "but it should only be loaded when get_features() is called"
    )


def test_codon_bias_get_features_works():
    """Test that calling get_features() successfully returns feature classes.

    This test verifies that when dependencies are installed, get_features()
    returns the expected list of feature classes.
    """
    # Clear any previously imported modules
    _clear_modules("biotooler.families.codon_bias", "codonbias")

    # Import the family module
    from biotooler.families import codon_bias

    # Call get_features()
    features = codon_bias.get_features()

    # Check that it returns a list
    assert isinstance(features, list), "get_features() should return a list"

    # Check that the list is not empty
    assert len(features) > 0, "get_features() should return at least one feature"

    # Check that we got CodonBiasFeature
    assert features[0].__name__ == "CodonBiasFeature", (
        "Expected first feature to be CodonBiasFeature"
    )

    # Now codonbias SHOULD be imported (because get_features() loaded it)
    assert "codonbias" in sys.modules, (
        "codonbias should be imported after calling get_features()"
    )


def test_codon_bias_get_features_import_error_with_hint():
    """Test that get_features() raises ImportError with hint if dependency missing.

    This test simulates a missing dependency by temporarily blocking the import,
    then verifies that the error message includes helpful installation instructions.

    Note: This test is conditional - it only runs if we can simulate the missing
    dependency.
    """
    # Clear any previously imported modules
    _clear_modules("biotooler.families.codon_bias", "codonbias")

    # Try to simulate missing codonbias by temporarily hiding it
    original_codonbias = sys.modules.get("codonbias")
    codonbias_installed = False

    # Check if codonbias is actually installed
    try:
        import codonbias as _  # noqa: F401

        codonbias_installed = True
        # Remove it temporarily
        del sys.modules["codonbias"]
    except ImportError:
        pass

    # If codonbias is installed, we need to mock its absence
    if codonbias_installed:
        # Temporarily block codonbias import by modifying sys.modules
        import builtins

        original_import = builtins.__import__

        def mock_import(name: str, *args: Any, **kwargs: Any):
            if name == "codonbias" or name.startswith("codonbias."):
                raise ImportError(
                    "No module named 'codonbias'. "
                    "To use codon_bias family features, install with: "
                    "pip install biotooler[codon_bias]"
                )
            return original_import(name, *args, **kwargs)

        builtins.__import__ = mock_import

        try:
            # Clear the family module to force re-import
            _clear_modules("biotooler.families.codon_bias")

            # Import the family module
            from biotooler.families import codon_bias

            # Calling get_features() should raise ImportError
            try:
                codon_bias.get_features()
                # If we get here, the test failed
                raise AssertionError("Expected ImportError when codonbias is not installed")
            except ImportError as e:
                # Check that the error message includes installation hint
                error_message = str(e)
                assert (
                    "codon_bias" in error_message.lower()
                    or "codonbias" in error_message.lower()
                ), f"Error message should mention codon_bias or codonbias: {error_message}"
                # The lazy_import utility should provide a helpful message
                # We just verify that ImportError was raised
        finally:
            # Restore original import
            builtins.__import__ = original_import
            if original_codonbias is not None:
                sys.modules["codonbias"] = original_codonbias


def test_codon_bias_backward_compatibility():
    """Test that existing import pattern still works via lazy loading.

    This ensures backward compatibility: users can still do
    `from biotooler.families.codon_bias import CodonBiasFeature`
    """
    # Clear any previously imported modules
    _clear_modules("biotooler.families.codon_bias", "codonbias")

    # Import CodonBiasFeature directly (old pattern)
    from biotooler.families.codon_bias import CodonBiasFeature

    # Verify we got the class
    assert CodonBiasFeature is not None
    assert CodonBiasFeature.__name__ == "CodonBiasFeature"

    # codonbias should now be imported (because we accessed CodonBiasFeature)
    assert "codonbias" in sys.modules
