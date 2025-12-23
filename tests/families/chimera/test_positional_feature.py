"""Tests for ChimeraFeature PositionalFeature protocol implementation.

Verifies that ChimeraFeature implements the PositionalFeature protocol correctly
and that windowing is properly gated to use positional vectors.
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features.aggregation import AggregationSpec, PositionSpace


def test_chimera_feature_has_position_space_property():
    """Test that ChimeraFeature has position_space property."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature with minimal setup
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Check position_space property exists and returns correct value
    assert hasattr(feature, "position_space")
    assert feature.position_space == PositionSpace.CODON


def test_chimera_feature_has_vector_keys_property():
    """Test that ChimeraFeature has vector_keys property."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature with cARS algorithm
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Check vector_keys property exists
    assert hasattr(feature, "vector_keys")
    assert isinstance(feature.vector_keys, dict)

    # Check that it returns the correct key
    assert "cARS_score" in feature.vector_keys
    assert isinstance(feature.vector_keys["cARS_score"], AggregationSpec)

    # Check that the aggregation function is mean
    assert feature.vector_keys["cARS_score"].aggregation_fn == np.mean


def test_chimera_feature_vector_keys_for_pscars():
    """Test that vector_keys returns correct key for PScARS algorithm."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature with PScARS algorithm
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="PScARS")

    # Check that it returns the correct key
    assert "PScARS_score" in feature.vector_keys
    assert isinstance(feature.vector_keys["PScARS_score"], AggregationSpec)


def test_chimera_feature_has_compute_vector_method():
    """Test that ChimeraFeature has compute_vector method."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Check compute_vector method exists
    assert hasattr(feature, "compute_vector")
    assert callable(feature.compute_vector)


def test_chimera_feature_compute_vector_works():
    """Test that compute_vector returns per-position values."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Create a test record
    record = SeqRecord(Seq("ATGAAATAA"), id="test")

    # compute_vector should return per-position values
    result = feature.compute_vector(record)

    # Check the result
    assert "cARS_score" in result
    assert isinstance(result["cARS_score"], np.ndarray)
    # 9 nucleotides = 3 codons
    assert result["cARS_score"].shape == (3,)


def test_chimera_feature_implements_positional_protocol():
    """Test that ChimeraFeature implements all required PositionalFeature protocol methods."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Check all required protocol properties and methods exist
    assert hasattr(feature, "position_space")
    assert hasattr(feature, "vector_keys")
    assert hasattr(feature, "compute_vector")

    # Verify they are the right types
    assert isinstance(feature.position_space, PositionSpace)
    assert isinstance(feature.vector_keys, dict)
    assert callable(feature.compute_vector)


def test_chimera_feature_backward_compatibility():
    """Test that ChimeraFeature maintains backward compatibility with __call__."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create feature
    ref_seqs = ["ATGAAATAA", "ATGGCATAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Check __call__ method still exists
    assert callable(feature)

    # Create a test record
    record = SeqRecord(Seq("ATGAAATAA"), id="test")

    # __call__ should now work and return a dict
    result = feature(record)
    assert isinstance(result, dict)
    assert "cARS_score" in result
