"""Test to verify biotooler imports remain lightweight."""

import sys


def _clear_modules(*prefixes: str) -> None:
    """Helper to clear modules with given prefixes from sys.modules.

    Args:
        prefixes: Module name prefixes to clear (e.g., "codonbias", "biotooler")
    """
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if any(module.startswith(prefix) for prefix in prefixes)
    ]
    for module in modules_to_clear:
        del sys.modules[module]


def test_import_biotooler_fast():
    """Test that importing biotooler does not import heavy optional dependencies.

    This test ensures that optional dependencies like codonbias are not loaded
    when just importing the biotooler package. Heavy dependencies should only
    be imported when their specific features are explicitly used.
    """
    # Clear any previously imported modules that we want to test
    _clear_modules("biotooler", "codonbias")

    # Import biotooler and biotooler.features
    import biotooler
    import biotooler.features

    # Check that heavy dependencies are NOT imported
    # codonbias should not be loaded yet
    assert "codonbias" not in sys.modules, (
        "codonbias was imported when importing biotooler, "
        "but it should only be loaded when explicitly used"
    )

    # Verify we can access the package
    assert biotooler.__version__ == "0.1.0"

    # Verify that FeatureSet and IncrementalFeature are accessible
    # (these are lightweight, no heavy dependencies)
    from biotooler.features import FeatureSet, IncrementalFeature

    assert FeatureSet is not None
    assert IncrementalFeature is not None

    # codonbias still should not be imported
    assert "codonbias" not in sys.modules, (
        "codonbias was imported when accessing FeatureSet/IncrementalFeature, "
        "but it should only be loaded when explicitly used"
    )


def test_import_codon_bias_feature_loads_codonbias():
    """Test that importing CodonBiasFeature does load codonbias.

    This is the opposite test - when we explicitly use CodonBiasFeature,
    the codonbias dependency should be loaded.
    """
    # Clear any previously imported modules
    _clear_modules("codonbias", "biotooler.families.codon_bias")

    # Now import CodonBiasFeature
    from biotooler.families.codon_bias import CodonBiasFeature

    # Now codonbias SHOULD be imported
    assert "codonbias" in sys.modules, (
        "codonbias should be imported when CodonBiasFeature is explicitly imported"
    )
    assert CodonBiasFeature is not None


def test_lazy_import_through_features_module():
    """Test that lazy import through features module works correctly."""
    # Clear any previously imported modules
    _clear_modules("codonbias", "biotooler.families.codon_bias", "biotooler.features")

    # Import biotooler.features
    import biotooler.features

    # codonbias should NOT be loaded yet
    assert "codonbias" not in sys.modules

    # Now access CodonBiasFeature through lazy import
    CodonBiasFeature = biotooler.features.CodonBiasFeature

    # Now codonbias SHOULD be loaded
    assert "codonbias" in sys.modules
    assert CodonBiasFeature is not None
