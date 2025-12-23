"""Tests for lazy import behavior of viennarna family.

Verifies that importing biotooler.families.viennarna does not eagerly load
ViennaRNA (RNA module), and that get_features() raises appropriate ImportError when
ViennaRNA is not installed.
"""

import sys

import pytest


def test_import_viennarna_does_not_load_rna():
    """Test that importing biotooler.families.viennarna doesn't import RNA.

    This ensures the family module itself is lightweight and only loads
    the heavy ViennaRNA dependency when features are actually accessed.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if (
            module.startswith("biotooler.families.viennarna")
            or module == "RNA"
            or module.startswith("RNA.")
        )
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the viennarna family module
    import biotooler.families.viennarna

    # Verify RNA is NOT loaded
    assert "RNA" not in sys.modules, (
        "RNA was imported when importing biotooler.families.viennarna, "
        "but it should only be loaded when features are accessed"
    )

    # Verify we can access FAMILY_META without loading RNA
    assert biotooler.families.viennarna.FAMILY_META["name"] == "viennarna"
    assert biotooler.families.viennarna.FAMILY_META["extra"] == "viennarna"

    # RNA still should not be imported
    assert "RNA" not in sys.modules


def test_get_features_raises_import_error_without_viennarna():
    """Test that get_features() raises ImportError with install instructions.

    When ViennaRNA is not installed, calling get_features() should raise
    an ImportError with a helpful message about how to install it.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if (
            module.startswith("biotooler.families.viennarna")
            or module == "RNA"
            or module.startswith("RNA.")
        )
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the viennarna family module
    from biotooler.families.viennarna import get_features

    # Calling get_features() should trigger lazy_import which raises ImportError
    with pytest.raises(ImportError) as exc_info:
        get_features()

    # Check the error message contains the required information
    error_message = str(exc_info.value)
    assert "RNA" in error_message or "viennarna" in error_message.lower(), (
        "Error message should mention RNA or viennarna"
    )
    assert 'pip install "biotooler[viennarna]"' in error_message, (
        "Error message should include installation instructions"
    )


def test_require_viennarna_raises_import_error_without_viennarna():
    """Test that require_viennarna() raises ImportError without ViennaRNA.

    When ViennaRNA is not installed, calling require_viennarna() should raise
    an ImportError with installation instructions.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if (
            module.startswith("biotooler.families.viennarna")
            or module == "RNA"
            or module.startswith("RNA.")
        )
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import require_viennarna function
    from biotooler.families.viennarna.integration import require_viennarna

    # Calling require_viennarna should raise ImportError
    with pytest.raises(ImportError) as exc_info:
        require_viennarna()

    # Check the error message
    error_message = str(exc_info.value)
    assert "RNA" in error_message, "Error message should mention RNA module"
    assert 'pip install "biotooler[viennarna]"' in error_message, (
        "Error message should include installation instructions"
    )
