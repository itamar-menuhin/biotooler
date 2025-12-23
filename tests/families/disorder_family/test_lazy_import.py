"""Tests for lazy import behavior of disorder family.

Verifies that importing biotooler.families.disorder does not eagerly load
metapredict or idrpred, and that get_features() raises appropriate ImportError when
dependencies are not installed.
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


def test_import_does_not_load_idrpred():
    """Test that importing biotooler.families.disorder doesn't import idrpred.

    This ensures the family module itself is lightweight and only loads
    the idrpred dependency when features are actually accessed.
    """
    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.disorder") or module.startswith("idrpred")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the disorder family module
    import biotooler.families.disorder

    # Verify idrpred is NOT loaded
    assert "idrpred" not in sys.modules, (
        "idrpred was imported when importing biotooler.families.disorder, "
        "but it should only be loaded when features are accessed"
    )

    # Verify we can access FAMILY_META without loading idrpred
    assert biotooler.families.disorder.FAMILY_META["name"] == "disorder"
    assert biotooler.families.disorder.FAMILY_META["extra"] == "disorder"

    # idrpred still should not be imported
    assert "idrpred" not in sys.modules


def test_get_features_returns_feature_list():
    """Test that get_features() returns list of feature classes.

    The disorder family now has features implemented, but get_features()
    should still not require metapredict to be installed at import time.
    metapredict is only loaded when features are actually instantiated and used.
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

    # Calling get_features() should return a list of feature classes
    features = get_features()
    assert isinstance(features, list), "Expected get_features() to return a list"
    assert len(features) > 0, "Expected at least one feature to be implemented"

    # Verify features are classes, not instances
    for feature in features:
        assert isinstance(feature, type), "Expected feature classes, not instances"

    # metapredict should still not be imported (lazy import happens in compute_vector)
    assert "metapredict" not in sys.modules


def test_require_metapredict_raises_import_error():
    """Test that require_metapredict() raises ImportError with install instructions.

    The require_metapredict() helper from integration module should raise
    an ImportError with a helpful message about how to install metapredict.
    This will be used by future feature implementations.

    This test only runs when metapredict is NOT installed.
    """
    # Check if metapredict is available
    try:
        import metapredict  # noqa: F401

        pytest.skip("Test only valid when metapredict is not installed")
    except ImportError:
        pass

    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if (module.startswith("biotooler.families.disorder") or module.startswith("metapredict"))
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


def test_missing_idrpred_cli_error_mentions_extra():
    """Test that require_idrpred_cli() raises ImportError with install instructions.

    The require_idrpred_cli() helper should check for the idrpred CLI tool
    and raise an ImportError with helpful message about how to install it.

    This test only runs when idrpred CLI is NOT available on PATH.
    """
    import shutil

    # Check if idrpred CLI is available
    if shutil.which("idrpred") is not None:
        pytest.skip("Test only valid when idrpred CLI is not on PATH")

    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.disorder")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the integration module
    from biotooler.families.disorder.integration import require_idrpred_cli

    # Calling require_idrpred_cli() should raise ImportError
    with pytest.raises(ImportError) as exc_info:
        require_idrpred_cli()

    # Check the error message contains the required information
    error_message = str(exc_info.value)
    assert "idrpred" in error_message.lower(), "Error message should mention idrpred"
    assert 'pip install "biotooler[disorder-idrpred]"' in error_message, (
        "Error message should include installation instructions with disorder-idrpred extra"
    )
    assert "PATH" in error_message, "Error message should mention PATH"


def test_missing_idrpred_pkg_error_mentions_extra():
    """Test that require_idrpred_pkg() raises ImportError with install instructions.

    The require_idrpred_pkg() helper should check for the idrpred package
    and raise an ImportError with helpful message about how to install it.

    This test only runs when idrpred package is NOT installed.
    """
    # Check if idrpred package is available
    try:
        import idrpred  # noqa: F401

        pytest.skip("Test only valid when idrpred package is not installed")
    except ImportError:
        pass

    # Clear any previously imported modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if (module.startswith("biotooler.families.disorder") or module.startswith("idrpred"))
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import the integration module
    from biotooler.families.disorder.integration import require_idrpred_pkg

    # Calling require_idrpred_pkg() should raise ImportError
    with pytest.raises(ImportError) as exc_info:
        require_idrpred_pkg()

    # Check the error message contains the required information
    error_message = str(exc_info.value)
    assert "idrpred" in error_message.lower(), "Error message should mention idrpred"
    assert 'pip install "biotooler[disorder-idrpred]"' in error_message, (
        "Error message should include installation instructions with disorder-idrpred extra"
    )
