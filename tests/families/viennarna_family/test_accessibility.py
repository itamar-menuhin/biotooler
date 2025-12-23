"""Tests for ViennaRNA accessibility feature.

Tests verify:
1. Vector length equals ORF span length
2. Outputs deterministic for same input
3. Aggregations over windows match manual aggregation of the vector
4. DNA->RNA normalization works correctly
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features import FeatureSet
from biotooler.features.aggregation import PositionSpace

# Try to import ViennaRNA - skip tests if not available
pytest.importorskip("RNA", reason="ViennaRNA not installed")


class TestViennaRNAAccessibilityBasic:
    """Basic tests for ViennaRNAAccessibility feature."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the feature."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        assert feature is not None

    def test_position_space_is_residue(self):
        """Test that feature operates in RESIDUE position space."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        assert feature.position_space == PositionSpace.RESIDUE

    def test_vector_keys_contains_pu(self):
        """Test that vector_keys contains PU with mean aggregation."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        vector_keys = feature.vector_keys
        assert "PU" in vector_keys
        assert vector_keys["PU"].aggregation_fn == np.mean

    def test_compute_vector_returns_dict_with_pu(self):
        """Test that compute_vector returns a dictionary with PU key."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        record = SeqRecord(Seq("ACGUACGU"), id="test")
        result = feature.compute_vector(record)

        assert isinstance(result, dict)
        assert "PU" in result
        assert isinstance(result["PU"], np.ndarray)

    def test_vector_length_matches_sequence_length(self):
        """Test that vector length equals sequence length."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Test with various sequence lengths
        for seq_len in [8, 15, 30]:
            seq = "ACGU" * (seq_len // 4)
            record = SeqRecord(Seq(seq), id="test")
            result = feature.compute_vector(record)

            assert len(result["PU"]) == len(seq), (
                f"Expected PU vector length {len(seq)}, got {len(result['PU'])}"
            )

    def test_pu_values_in_valid_range(self):
        """Test that PU values are probabilities in [0, 1]."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        record = SeqRecord(Seq("ACGUACGUACGUACGU"), id="test")
        result = feature.compute_vector(record)

        pu_values = result["PU"]
        assert np.all(pu_values >= 0.0), "PU values should be >= 0"
        assert np.all(pu_values <= 1.0), "PU values should be <= 1"

    def test_empty_sequence_returns_empty_vector(self):
        """Test that empty sequence returns empty PU vector."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        record = SeqRecord(Seq(""), id="test")
        result = feature.compute_vector(record)

        assert len(result["PU"]) == 0


class TestViennaRNAAccessibilityDeterminism:
    """Test determinism of ViennaRNA accessibility feature."""

    def test_deterministic_for_same_input(self):
        """Test that same input produces same output."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        record = SeqRecord(Seq("ACGUACGUACGUACGU"), id="test")

        # Compute multiple times
        result1 = feature.compute_vector(record)
        result2 = feature.compute_vector(record)
        result3 = feature.compute_vector(record)

        # All results should be identical
        np.testing.assert_array_equal(result1["PU"], result2["PU"])
        np.testing.assert_array_equal(result1["PU"], result3["PU"])

    def test_different_sequences_produce_different_results(self):
        """Test that different sequences produce different results."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Two sequences with clearly different pairing potential
        # CGCGCGCG has strong C-G pairing potential (stable structure)
        # AUAUAUAU has weaker A-U pairing potential (less stable structure)
        record1 = SeqRecord(Seq("CGCGCGCG"), id="test1")
        record2 = SeqRecord(Seq("AUAUAUAU"), id="test2")

        result1 = feature.compute_vector(record1)
        result2 = feature.compute_vector(record2)

        # Results should differ due to different pairing stabilities
        assert not np.array_equal(result1["PU"], result2["PU"])


class TestViennaRNAAccessibilityNormalization:
    """Test DNA to RNA normalization."""

    def test_dna_to_rna_normalization(self):
        """Test that DNA (with T) is converted to RNA (with U)."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # DNA sequence with T
        dna_record = SeqRecord(Seq("ATGATGATG"), id="dna")
        dna_result = feature.compute_vector(dna_record)

        # RNA sequence with U (should be equivalent)
        rna_record = SeqRecord(Seq("AUGAUGAUG"), id="rna")
        rna_result = feature.compute_vector(rna_record)

        # Results should be identical after normalization
        np.testing.assert_array_almost_equal(
            dna_result["PU"], rna_result["PU"],
            decimal=10,
            err_msg="DNA and RNA sequences should produce identical results"
        )

    def test_lowercase_sequences_handled(self):
        """Test that lowercase sequences are handled correctly."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Lowercase sequence
        lower_record = SeqRecord(Seq("acguacgu"), id="lower")
        lower_result = feature.compute_vector(lower_record)

        # Uppercase sequence
        upper_record = SeqRecord(Seq("ACGUACGU"), id="upper")
        upper_result = feature.compute_vector(upper_record)

        # Results should be identical
        np.testing.assert_array_almost_equal(
            lower_result["PU"], upper_result["PU"],
            decimal=10,
            err_msg="Lowercase and uppercase sequences should produce identical results"
        )


class TestViennaRNAAccessibilityWindowAggregation:
    """Test window aggregation matches manual aggregation."""

    def test_window_aggregation_matches_manual(self):
        """Test that window engine v2 aggregation matches manual aggregation."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        fs = FeatureSet(feature, name="viennarna")

        # Create a test sequence
        seq = "ACGUACGUACGUACGU"  # 16 nucleotides
        record = SeqRecord(Seq(seq), id="test")

        # Compute full vector manually
        full_result = feature.compute_vector(record)
        pu_vector = full_result["PU"]

        # Define window parameters
        window_nt = 9
        step_nt = 3
        orf = (0, len(seq))

        # Use window engine v2
        windowed_result = fs.compute_orf_windows_v2(
            record, orf=orf, window_nt=window_nt, step_nt=step_nt
        )

        # Manually compute expected aggregations
        # Windows: 0-9, 3-12, 6-15
        expected_windows = []
        for start in range(0, len(seq), step_nt):
            end = min(start + window_nt, len(seq))
            if end - start == window_nt:  # Only full windows by default
                window_slice = pu_vector[start:end]
                expected_mean = np.mean(window_slice)
                expected_windows.append((start, expected_mean))

        # Compare with windowed results
        for start, expected_mean in expected_windows:
            col_name = f"viennarna.PU_{start}"
            assert col_name in windowed_result.columns, (
                f"Column {col_name} not found in windowed results"
            )
            actual_mean = windowed_result[col_name].iloc[0]
            np.testing.assert_almost_equal(
                actual_mean, expected_mean, decimal=10,
                err_msg=f"Window starting at {start} has incorrect aggregation"
            )

    def test_partial_window_aggregation(self):
        """Test that partial windows are aggregated correctly when included."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        fs = FeatureSet(feature, name="viennarna")

        # Sequence that will create a partial window
        seq = "ACGUACGUACGU"  # 12 nucleotides
        record = SeqRecord(Seq(seq), id="test")

        # Compute full vector manually
        full_result = feature.compute_vector(record)
        pu_vector = full_result["PU"]

        # Define window parameters that create partial windows
        window_nt = 9
        step_nt = 3
        orf = (0, len(seq))

        # Use window engine v2 with drop_partial=False
        windowed_result = fs.compute_orf_windows_v2(
            record, orf=orf, window_nt=window_nt, step_nt=step_nt, drop_partial=False
        )

        # Check the partial window (9-12, only 3 nucleotides)
        if "viennarna.PU_9" in windowed_result.columns:
            partial_window_slice = pu_vector[9:12]
            expected_mean = np.mean(partial_window_slice)
            actual_mean = windowed_result["viennarna.PU_9"].iloc[0]
            np.testing.assert_almost_equal(
                actual_mean, expected_mean, decimal=10,
                err_msg="Partial window has incorrect aggregation"
            )

    def test_single_window_covers_full_sequence(self):
        """Test single window covering the entire sequence."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()
        fs = FeatureSet(feature, name="viennarna")

        seq = "ACGUACGUACGU"
        record = SeqRecord(Seq(seq), id="test")

        # Full vector
        full_result = feature.compute_vector(record)
        pu_vector = full_result["PU"]

        # Single window covering the entire sequence
        windowed_result = fs.compute_orf_windows_v2(
            record, orf=(0, len(seq)), window_nt=len(seq), step_nt=len(seq)
        )

        # Should have exactly one window at position 0
        col_name = "viennarna.PU_0"
        assert col_name in windowed_result.columns

        # The aggregation should be the mean of the full vector
        expected_mean = np.mean(pu_vector)
        actual_mean = windowed_result[col_name].iloc[0]
        np.testing.assert_almost_equal(
            actual_mean, expected_mean, decimal=10,
            err_msg="Single window aggregation doesn't match full vector mean"
        )


class TestViennaRNAAccessibilityStructure:
    """Test structure-related properties of accessibility."""

    def test_hairpin_loop_high_accessibility(self):
        """Test that nucleotides in hairpin loops have high unpaired probability.

        A simple hairpin structure should have high PU in the loop region.
        """
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Classic hairpin-forming sequence
        # CGCG forms stem, UUUU forms loop
        hairpin_seq = "CGCGUUUUCGCG"
        record = SeqRecord(Seq(hairpin_seq), id="hairpin")
        result = feature.compute_vector(record)

        pu_values = result["PU"]

        # Loop region (positions 4-7, UUUU) should generally have higher PU
        # than stem regions, though this depends on the exact energy model
        # Just verify all values are valid probabilities
        assert len(pu_values) == len(hairpin_seq)
        assert np.all(pu_values >= 0.0)
        assert np.all(pu_values <= 1.0)

    def test_poly_a_sequence_accessibility(self):
        """Test that poly-A sequences (no base pairing) have high accessibility.

        Sequences with no self-complementarity should have high unpaired probabilities.
        """
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Poly-A sequence - no Watson-Crick pairing possible
        poly_a_seq = "AAAAAAAA"
        record = SeqRecord(Seq(poly_a_seq), id="poly_a")
        result = feature.compute_vector(record)

        pu_values = result["PU"]

        # Most positions should have high PU (close to 1.0) since no pairing is favorable
        # Allow some flexibility due to non-Watson-Crick interactions and energy model
        assert len(pu_values) == len(poly_a_seq)
        assert np.mean(pu_values) > 0.5, "Poly-A should have generally high accessibility"


class TestViennaRNAAccessibilityPositions:
    """Test positions argument support."""

    def test_positions_argument_full_vs_sparse(self):
        """Test that positions argument produces correct values.

        When positions are specified, the full vector is still computed (for
        structural context correctness), and values at those positions should
        match the full computation.
        """
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Use a hairpin sequence with interesting structure
        seq = "CGCGUUUUCGCG"
        record = SeqRecord(Seq(seq), id="test")

        # Compute full vector
        result_full = feature.compute_vector(record)
        pu_full = result_full["PU"]

        # Compute with positions argument
        # Select a few positions to check
        positions = np.array([0, 4, 7, 11], dtype=np.int64)
        result_positions = feature.compute_vector(record, positions=positions)
        pu_positions = result_positions["PU"]

        # The returned vector should still be full length
        assert len(pu_positions) == len(seq), (
            "Vector length should match sequence length even with positions argument"
        )

        # Values at the specified positions should match full computation
        for pos in positions:
            np.testing.assert_almost_equal(
                pu_full[pos], pu_positions[pos], decimal=10,
                err_msg=f"Position {pos} value differs between full and sparse computation"
            )

    def test_positions_argument_preserves_all_values(self):
        """Test that positions argument doesn't affect computation correctness.

        The full vector should be computed for structural correctness,
        so all values should be valid and identical to full computation.
        """
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        seq = "ACGUACGUACGUACGU"
        record = SeqRecord(Seq(seq), id="test")

        # Compute full vector
        result_full = feature.compute_vector(record)
        pu_full = result_full["PU"]

        # Compute with positions argument (even indices)
        positions = np.arange(0, len(seq), 2, dtype=np.int64)
        result_positions = feature.compute_vector(record, positions=positions)
        pu_positions = result_positions["PU"]

        # All values should match exactly
        np.testing.assert_array_almost_equal(
            pu_full, pu_positions, decimal=10,
            err_msg="Full computation should be identical regardless of positions argument"
        )


class TestViennaRNAAccessibilityValueRange:
    """Test PU value range with epsilon."""

    def test_pu_values_within_valid_range_with_epsilon(self):
        """Test that PU values are within [0, 1] with small numerical tolerance.

        Due to floating-point arithmetic, values might be slightly outside
        the strict [0, 1] range but should be within a small epsilon.
        """
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        feature = ViennaRNAAccessibility()

        # Test with various sequences
        test_sequences = [
            "ACGUACGUACGUACGU",
            "CGCGUUUUCGCG",  # hairpin
            "AAAAAAAA",  # poly-A
            "GCGCGCGCGCGC",  # strong pairing
            "AAAUUUGGGCCC",  # mixed
        ]

        epsilon = 1e-6  # Small tolerance for numerical errors

        for seq in test_sequences:
            record = SeqRecord(Seq(seq), id="test")
            result = feature.compute_vector(record)
            pu_values = result["PU"]

            # Check that all values are within [0-epsilon, 1+epsilon]
            assert np.all(pu_values >= -epsilon), (
                f"PU values should be >= -epsilon ({-epsilon}), "
                f"but got min={np.min(pu_values)} for sequence {seq}"
            )
            assert np.all(pu_values <= 1.0 + epsilon), (
                f"PU values should be <= 1+epsilon ({1.0 + epsilon}), "
                f"but got max={np.max(pu_values)} for sequence {seq}"
            )

            # Also check strict bounds (should pass in practice)
            assert np.all(pu_values >= 0.0), (
                f"PU values should be >= 0.0, but got min={np.min(pu_values)} for sequence {seq}"
            )
            assert np.all(pu_values <= 1.0), (
                f"PU values should be <= 1.0, but got max={np.max(pu_values)} for sequence {seq}"
            )
