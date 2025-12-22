"""Tests for the families registry module."""

import sys


def test_import_families_does_not_import_family_modules():
    """Test that importing biotooler.families does not import family modules.

    This ensures the registry is lightweight and doesn't trigger import of
    heavy dependencies like codonbias.
    """
    # Clear any previously imported family modules
    modules_to_clear = [
        module
        for module in list(sys.modules.keys())
        if module.startswith("biotooler.families.codon_bias") or module.startswith("codonbias")
    ]
    for module in modules_to_clear:
        del sys.modules[module]

    # Import biotooler.families - this should import registry but not family modules
    import biotooler.families

    # Verify that we can use the registry functions
    assert hasattr(biotooler.families, "list_families")
    assert hasattr(biotooler.families, "get_family_meta")
    assert hasattr(biotooler.families, "get_install_extras")

    # Verify that family modules are NOT imported
    assert "biotooler.families.codon_bias" not in sys.modules, (
        "biotooler.families.codon_bias was imported when importing biotooler.families, "
        "but registry should be lightweight and not import family modules"
    )

    # Verify that heavy dependencies are NOT imported
    assert "codonbias" not in sys.modules, (
        "codonbias was imported when importing biotooler.families, "
        "but registry should be lightweight and not import optional dependencies"
    )


def test_list_families_returns_existing_families():
    """Test that list_families returns entries for existing families.

    This verifies that at least the codon_bias family is present and properly
    formatted in the registry.
    """
    from biotooler.families import list_families

    families = list_families()

    # Should return at least 1 family (codon_bias)
    assert len(families) >= 1, "Expected at least one family in the registry"

    # Check that codon_bias is present
    family_names = [f["name"] for f in families]
    assert "codon_bias" in family_names, "codon_bias family should be in registry"

    # Verify structure of family entries
    for family in families:
        assert "name" in family, "Family entry must have 'name' field"
        assert "summary" in family, "Family entry must have 'summary' field"
        assert "extra" in family, "Family entry must have 'extra' field"
        assert "owner" in family, "Family entry must have 'owner' field"
        assert "heavy" in family, "Family entry must have 'heavy' field"

        # Verify types
        assert isinstance(family["name"], str), "name must be a string"
        assert isinstance(family["summary"], str), "summary must be a string"
        assert family["extra"] is None or isinstance(family["extra"], str), (
            "extra must be None or string"
        )
        assert isinstance(family["owner"], str), "owner must be a string"
        assert isinstance(family["heavy"], bool), "heavy must be a boolean"

    # Verify families are sorted by name
    assert families == sorted(families, key=lambda f: f["name"]), (
        "Families should be sorted by name"
    )


def test_get_family_meta_unknown_raises_clear_message():
    """Test that get_family_meta raises ValueError with clear message for unknown family."""
    from biotooler.families import get_family_meta

    # Try to get metadata for a non-existent family
    try:
        get_family_meta("nonexistent_family")
        raise AssertionError("Expected ValueError to be raised for unknown family")
    except ValueError as e:
        error_message = str(e)
        # Verify error message is clear and helpful
        assert "nonexistent_family" in error_message, (
            "Error message should mention the unknown family name"
        )
        assert "Available families" in error_message or "available families" in error_message, (
            "Error message should mention available families"
        )
        assert "codon_bias" in error_message, (
            "Error message should list available families (including codon_bias)"
        )


def test_get_family_meta_returns_correct_data():
    """Test that get_family_meta returns correct metadata for known families."""
    from biotooler.families import get_family_meta

    # Get metadata for codon_bias
    meta = get_family_meta("codon_bias")

    # Verify structure
    assert meta["name"] == "codon_bias"
    assert "summary" in meta
    assert "extra" in meta
    assert "owner" in meta
    assert "heavy" in meta

    # Verify types
    assert isinstance(meta["summary"], str)
    assert meta["extra"] is None or isinstance(meta["extra"], str)
    assert isinstance(meta["owner"], str)
    assert isinstance(meta["heavy"], bool)


def test_get_install_extras_returns_correct_mapping():
    """Test that get_install_extras returns correct extra to families mapping."""
    from biotooler.families import get_install_extras

    extras = get_install_extras()

    # Should be a dict
    assert isinstance(extras, dict)

    # codon_bias has an extra, so it should be in the result
    assert "codon_bias" in extras, "codon_bias extra should be in the mapping"
    assert "codon_bias" in extras["codon_bias"], (
        "codon_bias family should be in the codon_bias extra list"
    )

    # Verify structure: values should be lists of strings
    for extra_name, family_list in extras.items():
        assert isinstance(extra_name, str), "Extra name must be a string"
        assert isinstance(family_list, list), "Family list must be a list"
        for family_name in family_list:
            assert isinstance(family_name, str), "Family name must be a string"
