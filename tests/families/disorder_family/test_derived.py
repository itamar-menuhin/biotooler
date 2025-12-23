"""Tests for derived disorder features computed from DISORDER_P.

Tests verify:
1. No metapredict dependency (uses synthetic DISORDER_P vectors)
2. Correct computation of all derived features
3. Edge cases (empty, all ordered, all disordered, single residue)
4. Threshold behavior
5. Longest IDR computation with various patterns
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestDisorderDerivedScalarsBasic:
    """Basic tests for DisorderDerivedScalars feature."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the feature."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        assert feature is not None

    def test_default_threshold(self):
        """Test that default threshold is 0.5."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        assert feature.threshold == 0.5

    def test_custom_threshold(self):
        """Test that custom threshold can be set."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.7)
        assert feature.threshold == 0.7

    def test_threshold_validation(self):
        """Test that threshold must be in range [0, 1]."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        # Valid thresholds
        DisorderDerivedScalars(threshold=0.0)
        DisorderDerivedScalars(threshold=0.5)
        DisorderDerivedScalars(threshold=1.0)

        # Invalid thresholds
        with pytest.raises(ValueError, match="threshold must be in range"):
            DisorderDerivedScalars(threshold=-0.1)
        with pytest.raises(ValueError, match="threshold must be in range"):
            DisorderDerivedScalars(threshold=1.1)

    def test_returns_dict_with_all_keys(self):
        """Test that compute returns dictionary with all expected keys."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Synthetic DISORDER_P vector
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.6, 0.7, 0.8, 0.3, 0.2, 0.1, 0.9])
        result = feature(record)

        assert isinstance(result, dict)
        assert "DISORDER_FRAC" in result
        assert "DISORDER_LONGEST_IDR" in result
        assert "DISORDER_MEAN" in result
        assert "DISORDER_P95" in result

    def test_missing_disorder_p_raises_error(self):
        """Test that missing DISORDER_P in annotations raises error."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # No DISORDER_P in annotations

        with pytest.raises(ValueError, match="DISORDER_P not found"):
            feature(record)


class TestDisorderDerivedScalarsComputation:
    """Test computation of individual derived features."""

    def test_disorder_frac_all_disordered(self):
        """Test DISORDER_FRAC when all residues are disordered."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")
        # All values above threshold
        record.annotations["DISORDER_P"] = np.array([0.6, 0.7, 0.8, 0.9, 1.0])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 1.0

    def test_disorder_frac_all_ordered(self):
        """Test DISORDER_FRAC when all residues are ordered."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")
        # All values below threshold
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.3, 0.4, 0.0])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 0.0

    def test_disorder_frac_mixed(self):
        """Test DISORDER_FRAC with mixed disorder."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # 4 out of 9 above threshold (0.5)
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.6, 0.7, 0.8, 0.3, 0.2, 0.1, 0.9])
        result = feature(record)

        expected_frac = 4.0 / 9.0
        assert abs(result["DISORDER_FRAC"] - expected_frac) < 1e-10

    def test_disorder_frac_at_threshold(self):
        """Test DISORDER_FRAC with values exactly at threshold (should be included)."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")
        # 3 values at or above threshold (>= 0.5)
        record.annotations["DISORDER_P"] = np.array([0.4, 0.5, 0.5, 0.6, 0.3])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 0.6  # 3/5

    def test_disorder_frac_different_threshold(self):
        """Test DISORDER_FRAC with different threshold."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.7)
        record = SeqRecord(Seq("MKALV"), id="test")
        # 2 out of 5 above 0.7
        record.annotations["DISORDER_P"] = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 0.6  # 3/5 (0.7, 0.8, 0.9)

    def test_disorder_mean(self):
        """Test DISORDER_MEAN computation."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.2, 0.4, 0.6, 0.8, 1.0])
        result = feature(record)

        expected_mean = 0.6  # (0.2 + 0.4 + 0.6 + 0.8 + 1.0) / 5
        assert abs(result["DISORDER_MEAN"] - expected_mean) < 1e-10

    def test_disorder_p95(self):
        """Test DISORDER_P95 computation."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("M" * 20), id="test")
        # Create array where 95th percentile should be clear
        disorder_p = np.linspace(0.0, 1.0, 20)
        record.annotations["DISORDER_P"] = disorder_p
        result = feature(record)

        # 95th percentile of 20 values from 0 to 1 should be around 0.95
        expected_p95 = np.percentile(disorder_p, 95)
        assert abs(result["DISORDER_P95"] - expected_p95) < 1e-10


class TestDisorderLongestIDR:
    """Test DISORDER_LONGEST_IDR computation."""

    def test_longest_idr_single_run(self):
        """Test longest IDR with a single contiguous run."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Single run of 3 disordered residues
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.6, 0.7, 0.8, 0.3, 0.2, 0.1, 0.4])
        result = feature(record)

        assert result["DISORDER_LONGEST_IDR"] == 3

    def test_longest_idr_multiple_runs(self):
        """Test longest IDR with multiple runs (should return longest)."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALVSWGRPQM"), id="test")
        # Two runs: length 2 and length 4
        record.annotations["DISORDER_P"] = np.array(
            [0.1, 0.6, 0.7, 0.3, 0.2, 0.6, 0.7, 0.8, 0.9, 0.3, 0.1, 0.2]
        )
        result = feature(record)

        assert result["DISORDER_LONGEST_IDR"] == 4

    def test_longest_idr_all_disordered(self):
        """Test longest IDR when entire sequence is disordered."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.6, 0.7, 0.8, 0.9, 1.0])
        result = feature(record)

        assert result["DISORDER_LONGEST_IDR"] == 5

    def test_longest_idr_no_disordered(self):
        """Test longest IDR when no residues are disordered."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.3, 0.4, 0.0])
        result = feature(record)

        assert result["DISORDER_LONGEST_IDR"] == 0

    def test_longest_idr_single_residue_above(self):
        """Test longest IDR with single residue above threshold."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.6, 0.3, 0.2])
        result = feature(record)

        assert result["DISORDER_LONGEST_IDR"] == 1

    def test_longest_idr_at_boundaries(self):
        """Test longest IDR at sequence boundaries."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALV"), id="test")

        # Run at start
        record.annotations["DISORDER_P"] = np.array([0.6, 0.7, 0.3, 0.2, 0.1])
        result = feature(record)
        assert result["DISORDER_LONGEST_IDR"] == 2

        # Run at end
        record.annotations["DISORDER_P"] = np.array([0.1, 0.2, 0.3, 0.6, 0.7])
        result = feature(record)
        assert result["DISORDER_LONGEST_IDR"] == 2

    def test_longest_idr_alternating(self):
        """Test longest IDR with alternating pattern."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        # Alternating: longest run is 1
        record.annotations["DISORDER_P"] = np.array([0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.2, 0.5])
        result = feature(record)

        assert result["DISORDER_LONGEST_IDR"] == 1


class TestDisorderDerivedScalarsEdgeCases:
    """Test edge cases for derived disorder features."""

    def test_empty_sequence(self):
        """Test with empty sequence."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["DISORDER_P"] = np.array([])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 0.0
        assert result["DISORDER_LONGEST_IDR"] == 0
        assert np.isnan(result["DISORDER_MEAN"])
        assert np.isnan(result["DISORDER_P95"])

    def test_single_residue_above_threshold(self):
        """Test with single residue above threshold."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("M"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.8])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 1.0
        assert result["DISORDER_LONGEST_IDR"] == 1
        assert result["DISORDER_MEAN"] == 0.8
        assert result["DISORDER_P95"] == 0.8

    def test_single_residue_below_threshold(self):
        """Test with single residue below threshold."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars(threshold=0.5)
        record = SeqRecord(Seq("M"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.2])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 0.0
        assert result["DISORDER_LONGEST_IDR"] == 0
        assert result["DISORDER_MEAN"] == 0.2
        assert result["DISORDER_P95"] == 0.2

    def test_all_zeros(self):
        """Test with all zero disorder probabilities."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["DISORDER_P"] = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 0.0
        assert result["DISORDER_LONGEST_IDR"] == 0
        assert result["DISORDER_MEAN"] == 0.0
        assert result["DISORDER_P95"] == 0.0

    def test_all_ones(self):
        """Test with all maximum disorder probabilities."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        record.annotations["DISORDER_P"] = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        result = feature(record)

        assert result["DISORDER_FRAC"] == 1.0
        assert result["DISORDER_LONGEST_IDR"] == 5
        assert result["DISORDER_MEAN"] == 1.0
        assert result["DISORDER_P95"] == 1.0

    def test_disorder_p_as_list(self):
        """Test that DISORDER_P can be provided as a list (converted to array)."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        # Provide as list instead of numpy array
        record.annotations["DISORDER_P"] = [0.1, 0.2, 0.6, 0.7, 0.8]
        result = feature(record)

        # Should work (converted to array internally)
        assert result["DISORDER_FRAC"] == 0.6  # 3/5
        assert result["DISORDER_LONGEST_IDR"] == 3

    def test_invalid_disorder_p_shape(self):
        """Test that multi-dimensional DISORDER_P raises error."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        # 2D array (invalid)
        record.annotations["DISORDER_P"] = np.array([[0.1, 0.2], [0.3, 0.4]])

        with pytest.raises(ValueError, match="must be 1-dimensional"):
            feature(record)

    def test_invalid_disorder_p_type(self):
        """Test that invalid DISORDER_P type raises error."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALV"), id="test")
        # Invalid type
        record.annotations["DISORDER_P"] = "not an array"

        with pytest.raises(ValueError, match="could not be converted to array"):
            feature(record)


class TestDisorderDerivedScalarsExport:
    """Test that the derived feature is properly exported."""

    def test_feature_in_get_features(self):
        """Test that DisorderDerivedScalars is in get_features() list."""
        from biotooler.families.disorder import get_features
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        features = get_features()
        assert DisorderDerivedScalars in features

    def test_get_features_returns_multiple_classes(self):
        """Test that get_features() returns both predictor and derived features."""
        from biotooler.families.disorder import get_features

        features = get_features()
        assert isinstance(features, list)
        assert len(features) >= 2  # At least DisorderProfileMetapredict and DisorderDerivedScalars
        # Check that items are classes, not instances
        for feature_class in features:
            assert isinstance(feature_class, type)

    def test_can_import_derived_without_metapredict(self):
        """Test that derived features can be imported without metapredict installed."""
        # This test verifies that importing derived features doesn't trigger
        # metapredict import, which is important for users who only want
        # to compute derived features from existing DISORDER_P data
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        # Should succeed even if metapredict is not installed
        feature = DisorderDerivedScalars()
        assert feature is not None


class TestDisorderDerivedScalarsNoMetapredictCall:
    """Test that derived features never call metapredict."""

    def test_no_metapredict_call(self):
        """Test that computing derived features does not call metapredict."""
        from biotooler.families.disorder.derived import DisorderDerivedScalars

        # This test ensures that DisorderDerivedScalars only uses DISORDER_P
        # from annotations and never calls a predictor
        feature = DisorderDerivedScalars()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")

        # Manually set DISORDER_P (simulating it was computed earlier)
        synthetic_disorder_p = np.array([0.1, 0.2, 0.6, 0.7, 0.8, 0.3, 0.2, 0.1, 0.9])
        record.annotations["DISORDER_P"] = synthetic_disorder_p

        # Compute derived features
        result = feature(record)

        # Verify results were computed from our synthetic data
        assert "DISORDER_FRAC" in result
        assert "DISORDER_LONGEST_IDR" in result
        assert "DISORDER_MEAN" in result
        assert "DISORDER_P95" in result

        # Verify mean matches our synthetic data
        expected_mean = float(np.mean(synthetic_disorder_p))
        assert abs(result["DISORDER_MEAN"] - expected_mean) < 1e-10
