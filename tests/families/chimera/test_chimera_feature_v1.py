"""Tests for ChimeraFeature implementation (v1).

Tests scalar computation, vector computation, error handling, and
RNA/DNA normalization for the chimera feature family.
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.reference_sequences import ReferenceSequenceSet


def test_chimera_scalar_computation_deterministic():
    """Test deterministic scalar cARS computation."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError as e:
        pytest.skip(f"chimera not installed: {e}")

    # Create reference sequences with known patterns
    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA", "ATGCCCCCCCTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Target sequence that partially matches references
    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")
    result = feature(record)

    # Verify result structure
    assert "cARS_score" in result
    assert isinstance(result["cARS_score"], (float, np.floating))

    # Score should be positive (matches exist in reference)
    assert result["cARS_score"] >= 0.0


def test_chimera_missing_reference_set_error():
    """Test ChimeraFeature raises error when neither reference_seqs nor reference_set provided."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    # Should raise ValueError if neither is provided
    with pytest.raises(ValueError, match="Either reference_seqs or reference_set must be provided"):
        ChimeraFeature()


def test_chimera_scalar_with_reference_set():
    """Test scalar computation with ReferenceSequenceSet."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    # Create reference set
    cds = {
        "gene1": "ATGAAAAAATAA",
        "gene2": "ATGGGGGGGTAA",
        "gene3": "ATGCCCCCCCTAA",
    }
    ref_set = ReferenceSequenceSet(cds)

    # Create feature with reference_set
    feature = ChimeraFeature(reference_set=ref_set, algorithm="cARS")

    # Compute on target sequence
    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")
    result = feature(record)

    assert "cARS_score" in result
    assert isinstance(result["cARS_score"], (float, np.floating))
    assert result["cARS_score"] >= 0.0


def test_chimera_rna_normalization():
    """Test that RNA input (with U) is properly normalized to DNA (with T)."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # DNA sequence
    dna_record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="dna")
    dna_result = feature(dna_record)

    # RNA sequence (same as DNA but with U instead of T)
    rna_record = SeqRecord(Seq("AUGAAAAGGGUAA"), id="rna")
    rna_result = feature(rna_record)

    # Results should be identical after normalization
    assert dna_result["cARS_score"] == rna_result["cARS_score"]


def test_chimera_orf_nt_span():
    """Test subsequence extraction with orf_nt_span parameter."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Full sequence
    full_record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="full")
    full_result = feature(full_record)

    # Longer sequence with the same sequence embedded
    longer_record = SeqRecord(Seq("NNNNNNATGAAAAGGGTAANNNNNN"), id="longer")
    orf_result = feature(longer_record, orf_nt_span=(6, 19))

    # Results should be identical
    assert full_result["cARS_score"] == orf_result["cARS_score"]


def test_chimera_vector_computation():
    """Test that compute_vector returns per-position values."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA", "ATGCCCCCCCTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Compute vector
    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")
    vector_result = feature.compute_vector(record)

    # Check result structure
    assert "cARS_score" in vector_result
    assert isinstance(vector_result["cARS_score"], np.ndarray)

    # Vector should have one value per codon (13 nt / 3 = 4 codons)
    expected_length = len(str(record.seq)) // 3
    assert vector_result["cARS_score"].shape == (expected_length,)

    # All values should be non-negative
    assert np.all(vector_result["cARS_score"] >= 0.0)


def test_chimera_vector_scalar_consistency():
    """Test that scalar score equals mean of vector values."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA", "ATGCCCCCCCTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")

    # Compute both scalar and vector
    scalar_result = feature(record)
    vector_result = feature.compute_vector(record)

    # Scalar should be the mean of vector values
    expected_scalar = np.mean(vector_result["cARS_score"])
    assert np.isclose(scalar_result["cARS_score"], expected_scalar)


def test_chimera_pscars_algorithm():
    """Test position-specific cARS algorithm."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA"]
    win_params = {"size": 40, "center": 0, "by_start": True, "by_stop": True}
    feature = ChimeraFeature(
        reference_seqs=ref_seqs,
        algorithm="PScARS",
        win_params=win_params,
        max_len=40,
        max_pos=0.5,
    )

    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")
    result = feature(record)

    # Check result structure
    assert "PScARS_score" in result
    assert isinstance(result["PScARS_score"], (float, np.floating))
    assert result["PScARS_score"] >= 0.0


def test_chimera_suffix_array_caching():
    """Test that suffix array is cached between calls."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Initially no suffix array
    assert feature._suffix_array is None

    # First call should build suffix array
    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")
    feature(record)

    # Now suffix array should be cached
    assert feature._suffix_array is not None
    cached_sa = feature._suffix_array

    # Second call should reuse cached suffix array
    feature(record)
    assert feature._suffix_array is cached_sa  # Same object


def test_chimera_empty_sequence():
    """Test handling of empty or very short sequences."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cARS")

    # Very short sequence (less than one codon after conversion)
    short_record = SeqRecord(Seq("AT"), id="short")
    result = feature(short_record)

    # Should return 0 for empty codon sequence
    assert result["cARS_score"] == 0.0


def test_chimera_unsupported_algorithm():
    """Test that unsupported algorithms raise NotImplementedError."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    ref_seqs = ["ATGAAAAAATAA", "ATGGGGGGGTAA"]
    feature = ChimeraFeature(reference_seqs=ref_seqs, algorithm="cMap")

    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")

    # cMap is not yet implemented
    with pytest.raises(NotImplementedError, match="cMap.*not yet implemented"):
        feature(record)


def test_chimera_vector_with_reference_set():
    """Test vector computation with ReferenceSequenceSet."""
    try:
        from biotooler.families.chimera.feature import ChimeraFeature
    except ImportError:
        pytest.skip("chimera not installed")

    # Create reference set with proteins
    cds = {
        "gene1": "ATGAAAAAATAA",
        "gene2": "ATGGGGGGGTAA",
    }
    proteins = {
        "gene1": "MKK",
        "gene2": "MGG",
    }
    ref_set = ReferenceSequenceSet(cds, proteins)

    feature = ChimeraFeature(reference_set=ref_set, algorithm="cARS")

    record = SeqRecord(Seq("ATGAAAAGGGTAA"), id="test")
    vector_result = feature.compute_vector(record)

    assert "cARS_score" in vector_result
    assert isinstance(vector_result["cARS_score"], np.ndarray)
    assert len(vector_result["cARS_score"]) > 0
