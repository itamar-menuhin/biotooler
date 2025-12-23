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


def test_get_features_returns_empty_list():
    """Test that get_features() returns empty list when no features implemented.

    Since the disorder family has no features implemented yet, get_features()
    should return an empty list without requiring metapredict to be installed.
    When features are added, they will use lazy_import at their module level.
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

    # Calling get_features() should return empty list without requiring metapredict
    features = get_features()
    assert features == [], "Expected empty list as no features are implemented yet"

    # metapredict should still not be imported
    assert "metapredict" not in sys.modules


def test_require_metapredict_raises_import_error():
    """Test that require_metapredict() raises ImportError with install instructions.

    The require_metapredict() helper from integration module should raise
    an ImportError with a helpful message about how to install metapredict.
    This will be used by future feature implementations.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.disorder") or module.startswith("metapredict")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the integration module
    from biotooler.families.disorder.integration import require_metapredict

    # Calling require_metapredict() should raise ImportError
    with pytest.raises(ImportError) as exc_info:
        require_metapredict()

    # Check the error message contains the required information
    error_message = str(exc_info.value)
    assert "metapredict" in error_message.lower(), "Error message should mention metapredict"
    assert 'pip install "biotooler[disorder]"' in error_message, (
        "Error message should include installation instructions"
    )
