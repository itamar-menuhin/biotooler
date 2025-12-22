"""Tests for ChimeraFeature reference_set validation.

Verifies that ChimeraFeature correctly handles ReferenceSequenceSet parameter
and validates that it can provide the required sequence types.
"""

import pytest

from biotooler.core.reference_sequences import ReferenceSequenceSet


def test_chimera_feature_accepts_reference_set_with_proteins():
    """Test that ChimeraFeature accepts reference_set with protein sequences."""
    # This test requires pychimera to be installed, so we skip if not available
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a reference set with CDS and proteins
    cds = {
        "gene1": "ATGAAATAA",  # M K *
        "gene2": "ATGGCATAA",  # M A *
    }
    proteins = {
        "gene1": "MK",
        "gene2": "MA",
    }
    ref_set = ReferenceSequenceSet(cds, proteins)

    # Create feature with reference_set
    feature = ChimeraFeature(reference_set=ref_set, algorithm="cARS")

    # Verify it stored protein sequences
    assert feature._sequence_type == "protein"
    assert len(feature.reference_seqs) == 2
    assert "MK" in feature.reference_seqs
    assert "MA" in feature.reference_seqs


def test_chimera_feature_accepts_reference_set_with_cds_only():
    """Test that ChimeraFeature accepts reference_set with CDS only (derives proteins)."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a reference set with CDS only (proteins will be derived)
    cds = {
        "gene1": "ATGAAATAA",  # M K *
        "gene2": "ATGGCATAA",  # M A *
    }
    ref_set = ReferenceSequenceSet(cds)

    # Create feature with reference_set
    feature = ChimeraFeature(reference_set=ref_set, algorithm="cARS")

    # Should derive proteins from CDS
    assert feature._sequence_type == "protein"
    assert len(feature.reference_seqs) == 2


def test_chimera_feature_falls_back_to_cds_strings():
    """Test that ChimeraFeature falls back to CDS strings if protein derivation fails."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a reference set with CDS that has internal stops
    cds = {
        "gene1": "ATGTAATAA",  # M * * (internal stop)
    }
    ref_set = ReferenceSequenceSet(cds)

    # Create feature with reference_set - should fall back to CDS
    feature = ChimeraFeature(reference_set=ref_set, algorithm="cARS")

    # Should use CDS strings as fallback
    assert feature._sequence_type == "cds"
    assert len(feature.reference_seqs) == 1
    assert feature.reference_seqs[0] == "ATGTAATAA"


def test_chimera_feature_requires_reference_seqs_or_set():
    """Test that ChimeraFeature requires either reference_seqs or reference_set."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Should raise ValueError if neither is provided
    with pytest.raises(ValueError, match="Either reference_seqs or reference_set must be provided"):
        ChimeraFeature()


def test_chimera_feature_accepts_reference_seqs_list():
    """Test that ChimeraFeature still accepts legacy reference_seqs list."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature with legacy reference_seqs parameter
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Should store CDS sequences
    assert feature._sequence_type == "cds"
    assert len(feature.reference_seqs) == 2
    assert feature.reference_seqs == ref_seqs


def test_chimera_feature_prefers_protein_strings():
    """Test that ChimeraFeature prefers protein strings over CDS strings."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a reference set with both CDS and proteins
    cds = {
        "gene1": "ATGAAATAA",
    }
    proteins = {
        "gene1": "MK",
    }
    ref_set = ReferenceSequenceSet(cds, proteins)

    # Create feature with reference_set
    feature = ChimeraFeature(reference_set=ref_set, algorithm="cARS")

    # Should prefer proteins
    assert feature._sequence_type == "protein"
    assert feature.reference_seqs == ["MK"]
