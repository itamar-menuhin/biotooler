"""Tests for lazy import behavior of chimera family.

Verifies that importing biotooler.families.chimera does not eagerly load
pyChimera, and that get_features() raises appropriate ImportError when
pyChimera is not installed.
"""

import sys

import pytest


def test_import_chimera_does_not_load_pychimera():
    """Test that importing biotooler.families.chimera doesn't import pyChimera.

    This ensures the family module itself is lightweight and only loads
    the heavy pyChimera dependency when features are actually instantiated.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.chimera") or module.startswith("pyChimera")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the chimera family module
    import biotooler.families.chimera

    # Verify pyChimera is NOT loaded
    assert "pyChimera" not in sys.modules, (
        "pyChimera was imported when importing biotooler.families.chimera, "
        "but it should only be loaded when features are instantiated"
    )

    # Verify we can access FAMILY_META without loading pyChimera
    assert biotooler.families.chimera.FAMILY_META["name"] == "chimera"
    assert biotooler.families.chimera.FAMILY_META["extra"] == "chimera"

    # pyChimera still should not be imported
    assert "pyChimera" not in sys.modules


def test_get_features_raises_import_error_without_pychimera():
    """Test that get_features() raises ImportError with install instructions.

    When pyChimera is not installed, calling get_features() should raise
    an ImportError with a helpful message about how to install it.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.chimera") or module.startswith("pyChimera")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the chimera family module
    from biotooler.families.chimera import get_features

    # Calling get_features() should raise ImportError with install hint
    with pytest.raises(ImportError) as exc_info:
        get_features()

    # Check the error message contains the required information
    error_message = str(exc_info.value)
    assert "pyChimera" in error_message or "pychimera" in error_message.lower(), (
        "Error message should mention pyChimera"
    )
    assert 'pip install "biotooler[chimera]"' in error_message, (
        "Error message should include installation instructions"
    )


def test_feature_instantiation_raises_import_error_without_pychimera():
    """Test that instantiating ChimeraFeature raises ImportError without pyChimera.

    When pyChimera is not installed, attempting to create a ChimeraFeature
    instance should raise an ImportError with installation instructions.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.chimera") or module.startswith("pyChimera")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import ChimeraFeature (should work without loading pyChimera)
    from biotooler.families.chimera.feature import ChimeraFeature

    # Verify pyChimera is still not loaded after importing the class
    assert "pyChimera" not in sys.modules

    # Attempting to instantiate should raise ImportError
    with pytest.raises(ImportError) as exc_info:
        ChimeraFeature()

    # Check the error message
    error_message = str(exc_info.value)
    assert "pyChimera" in error_message or "pychimera" in error_message.lower()
    assert 'pip install "biotooler[chimera]"' in error_message
