"""Tests for derived IDRPred features computed from IDRPRED_IDR mask.

Tests verify:
1. No IDRPred CLI dependency (uses synthetic IDRPRED_IDR vectors)
2. Correct computation of all derived features
3. Edge cases (empty, all-0, all-1, single residue)
4. Various patterns for longest run and segment counting
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestIDRPredDerivedScalarsBasic:
    """Basic tests for IDRPredDerivedScalars feature."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the feature."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        assert feature is not None

    def test_returns_dict_with_all_keys(self):
        """Test that compute returns dictionary with all expected keys."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Synthetic IDRPRED_IDR vector
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 1, 1, 0, 0, 0, 1])
        result = feature(record)

        assert isinstance(result, dict)
        assert "IDRPRED_FRAC_IDR" in result
        assert "IDRPRED_LONGEST_IDR_LEN" in result
        assert "IDRPRED_NUM_IDR_SEGMENTS" in result

    def test_missing_idrpred_idr_raises_error(self):
        """Test that missing IDRPRED_IDR in annotations raises error."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # No IDRPRED_IDR in annotations

        with pytest.raises(ValueError, match="IDRPRED_IDR not found"):
            feature(record)


class TestIDRPredFractionEdgeCases:
    """Test IDRPRED_FRAC_IDR computation with edge cases."""

    def test_idrpred_fraction_empty_sequence(self):
        """Test with empty sequence."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 0.0

    def test_idrpred_fraction_all_zeros(self):
        """Test with all-zero mask (no IDRs)."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 0.0

    def test_idrpred_fraction_all_ones(self):
        """Test with all-one mask (entire sequence is IDR)."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 1, 1, 1])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 1.0

    def test_idrpred_fraction_single_residue_zero(self):
        """Test with single residue, no IDR."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("M"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 0.0

    def test_idrpred_fraction_single_residue_one(self):
        """Test with single residue, in IDR."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("M"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([1])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 1.0

    def test_idrpred_fraction_mixed_pattern(self):
        """Test with mixed IDR pattern."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # 4 out of 9 in IDRs
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 1, 1, 0, 0, 0, 1])
        result = feature(record)

        expected_frac = 4.0 / 9.0
        assert abs(result["IDRPRED_FRAC_IDR"] - expected_frac) < 1e-10

    def test_idrpred_fraction_half_and_half(self):
        """Test with exactly half in IDRs."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR" + "P"), id="test")
        # 5 out of 10 in IDRs
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 0.5


class TestIDRPredLongestRunCases:
    """Test IDRPRED_LONGEST_IDR_LEN computation."""

    def test_idrpred_longest_run_empty_sequence(self):
        """Test longest run with empty sequence."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 0

    def test_idrpred_longest_run_all_zeros(self):
        """Test longest run when no residues are in IDRs."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 0

    def test_idrpred_longest_run_all_ones(self):
        """Test longest run when entire sequence is IDR."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 1, 1, 1])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 5

    def test_idrpred_longest_run_single_residue_in_idr(self):
        """Test longest run with single residue in IDR."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 0, 0])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 1

    def test_idrpred_longest_run_single_contiguous_run(self):
        """Test longest run with a single contiguous run."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Single run of 3 residues in IDR
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 1, 1, 0, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 3

    def test_idrpred_longest_run_multiple_runs_return_longest(self):
        """Test longest run with multiple runs (should return longest)."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGRPQM"), id="test")
        # Two runs: length 2 and length 4
        record.annotations["IDRPRED_IDR"] = np.array([0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 4

    def test_idrpred_longest_run_at_start(self):
        """Test longest run at sequence start."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 2

    def test_idrpred_longest_run_at_end(self):
        """Test longest run at sequence end."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 1, 1])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 2

    def test_idrpred_longest_run_alternating(self):
        """Test longest run with alternating pattern."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Alternating: longest run is 1
        record.annotations["IDRPRED_IDR"] = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1])
        result = feature(record)

        assert result["IDRPRED_LONGEST_IDR_LEN"] == 1


class TestIDRPredSegmentCountCases:
    """Test IDRPRED_NUM_IDR_SEGMENTS computation."""

    def test_idrpred_segment_count_empty_sequence(self):
        """Test segment count with empty sequence."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 0

    def test_idrpred_segment_count_all_zeros(self):
        """Test segment count when no residues are in IDRs."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 0

    def test_idrpred_segment_count_all_ones(self):
        """Test segment count when entire sequence is one IDR."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 1, 1, 1])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

    def test_idrpred_segment_count_single_residue_segment(self):
        """Test segment count with single residue in IDR."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 0, 0])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

    def test_idrpred_segment_count_single_contiguous_segment(self):
        """Test segment count with a single contiguous segment."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Single segment of 3 residues
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 1, 1, 0, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

    def test_idrpred_segment_count_two_segments(self):
        """Test segment count with two distinct segments."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Two segments: [2-4] and [8]
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 1, 1, 1, 0, 0, 0, 1])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 2

    def test_idrpred_segment_count_three_segments(self):
        """Test segment count with three distinct segments."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGRPQM"), id="test")
        # Three segments: [1-2], [5-6], [10]
        record.annotations["IDRPRED_IDR"] = np.array([0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 3

    def test_idrpred_segment_count_alternating(self):
        """Test segment count with alternating pattern."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Alternating: 5 segments of length 1 each
        record.annotations["IDRPRED_IDR"] = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1])
        result = feature(record)

        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 5

    def test_idrpred_segment_count_at_boundaries(self):
        """Test segment count with segments at sequence boundaries."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")

        # Segment at start
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 0, 0, 0])
        result = feature(record)
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

        # Segment at end
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 1, 1])
        result = feature(record)
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

        # Segments at both boundaries
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 0, 1, 1])
        result = feature(record)
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 2


class TestIDRPredDerivedScalarsEdgeCases:
    """Test edge cases for derived IDRPred features."""

    def test_all_metrics_empty_sequence(self):
        """Test all metrics with empty sequence."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 0.0
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 0
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 0

    def test_all_metrics_all_zeros(self):
        """Test all metrics with all-zero mask."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([0, 0, 0, 0, 0])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 0.0
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 0
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 0

    def test_all_metrics_all_ones(self):
        """Test all metrics with all-one mask."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 1, 1, 1])
        result = feature(record)

        assert result["IDRPRED_FRAC_IDR"] == 1.0
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 5
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

    def test_idrpred_idr_as_list(self):
        """Test that IDRPRED_IDR can be provided as a list (converted to array)."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        # Provide as list instead of numpy array
        record.annotations["IDRPRED_IDR"] = [0, 0, 1, 1, 1]
        result = feature(record)

        # Should work (converted to array internally)
        assert result["IDRPRED_FRAC_IDR"] == 0.6  # 3/5
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 3
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 1

    def test_invalid_idrpred_idr_shape(self):
        """Test that multi-dimensional IDRPRED_IDR raises error."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        # 2D array (invalid)
        record.annotations["IDRPRED_IDR"] = np.array([[0, 1], [1, 0]])

        with pytest.raises(ValueError, match="must be 1-dimensional"):
            feature(record)

    def test_invalid_idrpred_idr_type(self):
        """Test that invalid IDRPRED_IDR type raises error."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        # Invalid type
        record.annotations["IDRPRED_IDR"] = "not an array"

        with pytest.raises(ValueError, match="could not be converted to array"):
            feature(record)


class TestIDRPredDerivedScalarsNoIDRPredCall:
    """Test that derived features never call IDRPred CLI."""

    def test_no_idrpred_call(self):
        """Test that computing derived features does not call IDRPred CLI."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        # This test ensures that IDRPredDerivedScalars only uses IDRPRED_IDR
        # from annotations and never calls the IDRPred CLI
        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")

        # Manually set IDRPRED_IDR (simulating it was computed earlier)
        synthetic_idrpred_idr = np.array([0, 0, 1, 1, 1, 0, 0, 0, 1])
        record.annotations["IDRPRED_IDR"] = synthetic_idrpred_idr

        # Compute derived features
        result = feature(record)

        # Verify results were computed from our synthetic data
        assert "IDRPRED_FRAC_IDR" in result
        assert "IDRPRED_LONGEST_IDR_LEN" in result
        assert "IDRPRED_NUM_IDR_SEGMENTS" in result

        # Verify fraction matches our synthetic data (4/9 = 0.444...)
        expected_frac = float(np.mean(synthetic_idrpred_idr))
        assert abs(result["IDRPRED_FRAC_IDR"] - expected_frac) < 1e-10


class TestIDRPredDerivedScalarsIntegration:
    """Integration tests for IDRPredDerivedScalars with complex patterns."""

    def test_complex_pattern_1(self):
        """Test with complex pattern: multiple segments of varying lengths."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("M" * 20), id="test")
        # Pattern: segments of length 1, 2, 3, 4 with gaps
        record.annotations["IDRPRED_IDR"] = np.array(
            [0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        )
        result = feature(record)

        # 10 out of 20 in IDRs
        assert result["IDRPRED_FRAC_IDR"] == 0.5
        # Longest segment is 4
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 4
        # Four segments
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 4

    def test_complex_pattern_2(self):
        """Test with complex pattern: many small segments."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("M" * 15), id="test")
        # Pattern: 5 segments of length 2 each with single gaps
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0])
        result = feature(record)

        # 10 out of 15 in IDRs
        assert abs(result["IDRPRED_FRAC_IDR"] - 10.0 / 15.0) < 1e-10
        # Longest segment is 2
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 2
        # Five segments
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 5

    def test_complex_pattern_3(self):
        """Test with complex pattern: large gap in the middle."""
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        feature = IDRPredDerivedScalars()
        record = SeqRecord(Seq("M" * 15), id="test")
        # Pattern: two segments with large gap
        record.annotations["IDRPRED_IDR"] = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        result = feature(record)

        # 10 out of 15 in IDRs
        assert abs(result["IDRPRED_FRAC_IDR"] - 10.0 / 15.0) < 1e-10
        # Longest segment is 5
        assert result["IDRPRED_LONGEST_IDR_LEN"] == 5
        # Two segments
        assert result["IDRPRED_NUM_IDR_SEGMENTS"] == 2


class TestIDRPredDerivedScalarsExport:
    """Test that the derived feature is properly exported."""

    def test_feature_in_get_features(self):
        """Test that IDRPredDerivedScalars is in get_features() list."""
        from biotooler.families.disorder import get_features
        from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

        features = get_features()
        assert IDRPredDerivedScalars in features

    def test_get_features_returns_multiple_classes(self):
        """Test that get_features() returns multiple feature classes."""
        from biotooler.families.disorder import get_features

        features = get_features()
        assert isinstance(features, list)
        assert len(features) >= 2  # At least predictor and derived features
        # Check that items are classes, not instances
        for feature_class in features:
            assert isinstance(feature_class, type)
