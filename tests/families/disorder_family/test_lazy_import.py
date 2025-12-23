"""Tests for lazy import behavior of disorder family.

Verifies that importing biotooler.families.disorder does not eagerly load
metapredict, and that get_features() raises appropriate ImportError when
metapredict is not installed.
"""

import sys

import pytest


def test_import_disorder_does_not_load_metapredict():
    """Test that importing biotooler.families.disorder doesn't import metapredict.

    This ensures the family module itself is lightweight and only loads
    the heavy metapredict dependency when features are actually accessed.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.disorder") or module.startswith("metapredict")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the disorder family module
    import biotooler.families.disorder

    # Verify metapredict is NOT loaded
    assert "metapredict" not in sys.modules, (
        "metapredict was imported when importing biotooler.families.disorder, "
        "but it should only be loaded when features are accessed"
    )

    # Verify we can access FAMILY_META without loading metapredict
    assert biotooler.families.disorder.FAMILY_META["name"] == "disorder"
    assert biotooler.families.disorder.FAMILY_META["extra"] == "disorder"

    # metapredict still should not be imported
    assert "metapredict" not in sys.modules


def test_get_features_raises_import_error_without_metapredict():
    """Test that get_features() raises ImportError with install instructions.

    When metapredict is not installed, calling get_features() should raise
    an ImportError with a helpful message about how to install it.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.disorder") or module.startswith("metapredict")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the disorder family module
    from biotooler.families.disorder import get_features

    # Calling get_features() triggers require_metapredict which should raise ImportError
    with pytest.raises(ImportError) as exc_info:
        get_features()

    # Check the error message contains the required information
    error_message = str(exc_info.value)
    assert "metapredict" in error_message.lower(), "Error message should mention metapredict"
    assert 'pip install "biotooler[disorder]"' in error_message, (
        "Error message should include installation instructions"
    )
