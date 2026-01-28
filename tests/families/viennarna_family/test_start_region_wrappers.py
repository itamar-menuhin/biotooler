"""Tests for ViennaRNA start region wrappers.

Tests verify that wrappers equal manually calling underlying features and aggregating.
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Try to import ViennaRNA - skip tests if not available
pytest.importorskip("RNA", reason="ViennaRNA not installed")


class TestWindowMFEStartRegion:
    """Tests for WindowMFEStartRegion wrapper."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the wrapper."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        wrapper = WindowMFEStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        assert wrapper is not None
        assert wrapper.region_start == 0
        assert wrapper.region_end == 30
        assert wrapper.window_size == 15
        assert wrapper.step == 10

    def test_parameter_validation(self):
        """Test that __init__ validates parameters."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        # Negative region_start
        with pytest.raises(ValueError, match="region_start must be non-negative"):
            WindowMFEStartRegion(region_start=-1, region_end=30, window_size=15, step=10)

        # region_end <= region_start
        with pytest.raises(ValueError, match="region_end must be greater than region_start"):
            WindowMFEStartRegion(region_start=30, region_end=30, window_size=15, step=10)
        with pytest.raises(ValueError, match="region_end must be greater than region_start"):
            WindowMFEStartRegion(region_start=30, region_end=20, window_size=15, step=10)

        # Non-positive window_size
        with pytest.raises(ValueError, match="window_size must be positive"):
            WindowMFEStartRegion(region_start=0, region_end=30, window_size=0, step=10)
        with pytest.raises(ValueError, match="window_size must be positive"):
            WindowMFEStartRegion(region_start=0, region_end=30, window_size=-5, step=10)

        # Non-positive step
        with pytest.raises(ValueError, match="step must be positive"):
            WindowMFEStartRegion(region_start=0, region_end=30, window_size=15, step=0)
        with pytest.raises(ValueError, match="step must be positive"):
            WindowMFEStartRegion(region_start=0, region_end=30, window_size=15, step=-5)

    def test_wrapper_equals_manual_underlying_feature(self):
        """Test that wrapper results equal manually calling WindowMFEFeature."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        # Create test sequence
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 40 bp
        record = SeqRecord(Seq(sequence), id="test")

        # Define region parameters
        region_start = 0
        region_end = 30
        window_size = 15
        step = 10

        # Compute using wrapper
        wrapper = WindowMFEStartRegion(
            region_start=region_start,
            region_end=region_end,
            window_size=window_size,
            step=step,
        )
        wrapper_result = wrapper(record)

        # Compute manually: generate window starts and call underlying feature
        manual_window_starts = list(range(region_start, region_end, step))
        manual_feature = WindowMFEFeature(
            window_starts=manual_window_starts, window_size=window_size
        )
        manual_result = manual_feature(record)

        # Results should be identical
        assert wrapper_result == manual_result

    def test_wrapper_generates_correct_window_starts(self):
        """Test that wrapper generates the correct window starts."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        sequence = "A" * 100
        record = SeqRecord(Seq(sequence), id="test")

        # Test case 1: step=10, region 0-30 -> windows at 0, 10, 20
        wrapper = WindowMFEStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        result = wrapper(record)
        assert "MFE_0" in result
        assert "MFE_10" in result
        assert "MFE_20" in result
        assert "MFE_30" not in result  # region_end is exclusive
        assert len(result) == 3

        # Test case 2: step=5, region 10-25 -> windows at 10, 15, 20
        wrapper = WindowMFEStartRegion(region_start=10, region_end=25, window_size=10, step=5)
        result = wrapper(record)
        assert "MFE_10" in result
        assert "MFE_15" in result
        assert "MFE_20" in result
        assert "MFE_25" not in result  # region_end is exclusive
        assert len(result) == 3

    def test_wrapper_handles_window_beyond_sequence(self):
        """Test that wrapper correctly handles windows extending beyond sequence."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        # Sequence is only 30 bp
        sequence = "ACGTACGTACGTACGTACGTACGTACGTAC"
        record = SeqRecord(Seq(sequence), id="test")

        # Region that would generate windows extending beyond sequence
        wrapper = WindowMFEStartRegion(region_start=0, region_end=40, window_size=20, step=15)
        result = wrapper(record)

        # Window at 0 (0-20) fits
        assert "MFE_0" in result
        # Window at 15 (15-35) extends beyond 30, should be skipped
        assert "MFE_15" not in result
        # Window at 30 is at the end, should be skipped
        assert "MFE_30" not in result

    def test_empty_region_returns_empty_dict(self):
        """Test that region with no windows returns empty dict."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
        record = SeqRecord(Seq(sequence), id="test")

        # When region_start=10, region_end=11, step=100 -> range(10, 11, 100) = [10]
        # This produces one window at position 10
        wrapper = WindowMFEStartRegion(region_start=10, region_end=11, window_size=10, step=100)
        result = wrapper(record)

        # One window at position 10 should be included (10 < 11)
        assert "MFE_10" in result

    def test_wrapper_no_hardcoded_defaults(self):
        """Test that wrapper requires all parameters (no hardcoded defaults)."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        # All parameters are required - this should work
        wrapper = WindowMFEStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        assert wrapper.region_start == 0
        assert wrapper.region_end == 30
        assert wrapper.window_size == 15
        assert wrapper.step == 10

    def test_wrapper_stable_naming(self):
        """Test that wrapper produces stable key names."""
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
        record = SeqRecord(Seq(sequence), id="test")

        wrapper = WindowMFEStartRegion(region_start=0, region_end=20, window_size=10, step=5)
        result = wrapper(record)

        # Check that keys follow the stable pattern MFE_<start>
        expected_keys = {"MFE_0", "MFE_5", "MFE_10", "MFE_15"}
        assert set(result.keys()) == expected_keys

    def test_wrapper_from_init_export(self):
        """Test that wrapper can be imported from __init__."""
        from biotooler.families.viennarna import WindowMFEStartRegion

        wrapper = WindowMFEStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        assert wrapper.region_start == 0

    def test_get_features_includes_wrapper(self):
        """Test that get_features() includes WindowMFEStartRegion."""
        from biotooler.families.viennarna import get_features
        from biotooler.families.viennarna.start_region_mfe import WindowMFEStartRegion

        features = get_features()
        assert WindowMFEStartRegion in features


class TestAccessibilityStartRegion:
    """Tests for AccessibilityStartRegion wrapper."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the wrapper."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        wrapper = AccessibilityStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        assert wrapper is not None
        assert wrapper.region_start == 0
        assert wrapper.region_end == 30
        assert wrapper.window_size == 15
        assert wrapper.step == 10

    def test_parameter_validation(self):
        """Test that __init__ validates parameters."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        # Negative region_start
        with pytest.raises(ValueError, match="region_start must be non-negative"):
            AccessibilityStartRegion(region_start=-1, region_end=30, window_size=15, step=10)

        # region_end <= region_start
        with pytest.raises(ValueError, match="region_end must be greater than region_start"):
            AccessibilityStartRegion(region_start=30, region_end=30, window_size=15, step=10)

        # Non-positive window_size
        with pytest.raises(ValueError, match="window_size must be positive"):
            AccessibilityStartRegion(region_start=0, region_end=30, window_size=0, step=10)

        # Non-positive step
        with pytest.raises(ValueError, match="step must be positive"):
            AccessibilityStartRegion(region_start=0, region_end=30, window_size=15, step=0)

    def test_wrapper_equals_manual_underlying_feature(self):
        """Test that wrapper results equal manually calling ViennaRNAAccessibility + aggregation."""
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        # Create test sequence
        sequence = "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"  # 40 nt
        record = SeqRecord(Seq(sequence), id="test")

        # Define region parameters
        region_start = 0
        region_end = 30
        window_size = 15
        step = 10

        # Compute using wrapper
        wrapper = AccessibilityStartRegion(
            region_start=region_start,
            region_end=region_end,
            window_size=window_size,
            step=step,
        )
        wrapper_result = wrapper(record)

        # Compute manually: call underlying feature and aggregate
        manual_feature = ViennaRNAAccessibility()
        vectors = manual_feature.compute_vector(record)
        pu_vector = vectors["PU"]

        manual_result = {}
        manual_window_starts = list(range(region_start, region_end, step))
        for window_start in manual_window_starts:
            window_end = min(window_start + window_size, len(pu_vector))
            if window_start < len(pu_vector):
                window_pu = pu_vector[window_start:window_end]
                if len(window_pu) > 0:
                    manual_result[f"PU_{window_start}"] = float(np.mean(window_pu))

        # Results should be identical
        assert set(wrapper_result.keys()) == set(manual_result.keys())
        for key in wrapper_result:
            np.testing.assert_almost_equal(
                wrapper_result[key],
                manual_result[key],
                decimal=10,
                err_msg=f"Mismatch for {key}",
            )

    def test_wrapper_generates_correct_window_starts(self):
        """Test that wrapper generates the correct window starts."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        sequence = "A" * 100
        record = SeqRecord(Seq(sequence), id="test")

        # Test case 1: step=10, region 0-30 -> windows at 0, 10, 20
        wrapper = AccessibilityStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        result = wrapper(record)
        assert "PU_0" in result
        assert "PU_10" in result
        assert "PU_20" in result
        assert "PU_30" not in result  # region_end is exclusive
        assert len(result) == 3

        # Test case 2: step=5, region 10-25 -> windows at 10, 15, 20
        wrapper = AccessibilityStartRegion(region_start=10, region_end=25, window_size=10, step=5)
        result = wrapper(record)
        assert "PU_10" in result
        assert "PU_15" in result
        assert "PU_20" in result
        assert "PU_25" not in result  # region_end is exclusive
        assert len(result) == 3

    def test_wrapper_handles_window_beyond_sequence(self):
        """Test that wrapper correctly handles windows extending beyond sequence."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        # Sequence is only 30 nt
        sequence = "ACGUACGUACGUACGUACGUACGUACGUAC"
        record = SeqRecord(Seq(sequence), id="test")

        # Region that would generate windows extending beyond sequence
        wrapper = AccessibilityStartRegion(region_start=0, region_end=40, window_size=20, step=15)
        result = wrapper(record)

        # Window at 0 (0-20) fits
        assert "PU_0" in result
        # Window at 15 (15-35) extends beyond 30, but should be clipped to 15-30
        assert "PU_15" in result
        # Window at 30 is at the end, should be skipped (empty window)
        assert "PU_30" not in result

    def test_empty_region_returns_empty_dict(self):
        """Test that region with no windows returns empty dict."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        sequence = "ACGUACGUACGUACGUACGUACGUACGUACGU"
        record = SeqRecord(Seq(sequence), id="test")

        # Region with region_start >= len(sequence) produces no windows
        # Window starts at position beyond sequence should be skipped
        wrapper = AccessibilityStartRegion(region_start=50, region_end=60, window_size=10, step=5)
        result = wrapper(record)

        assert result == {}

    def test_wrapper_no_recomputation(self):
        """Test that wrapper computes PU vector only once (no recomputation)."""
        from unittest.mock import patch

        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        sequence = "ACGUACGUACGUACGUACGUACGUACGUACGU"
        record = SeqRecord(Seq(sequence), id="test")

        # Patch compute_vector to count calls
        with patch(
            "biotooler.families.viennarna.accessibility.ViennaRNAAccessibility.compute_vector"
        ) as mock_compute:
            # Set up mock to return a valid PU vector
            mock_compute.return_value = {"PU": np.ones(len(sequence))}

            wrapper = AccessibilityStartRegion(
                region_start=0, region_end=30, window_size=10, step=5
            )
            wrapper(record)

            # compute_vector should be called exactly once, not once per window
            assert mock_compute.call_count == 1

    def test_wrapper_no_hardcoded_defaults(self):
        """Test that wrapper requires all parameters (no hardcoded defaults)."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        # All parameters are required - this should work
        wrapper = AccessibilityStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        assert wrapper.region_start == 0
        assert wrapper.region_end == 30
        assert wrapper.window_size == 15
        assert wrapper.step == 10

    def test_wrapper_stable_naming(self):
        """Test that wrapper produces stable key names."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        sequence = "ACGUACGUACGUACGUACGUACGUACGUACGU"
        record = SeqRecord(Seq(sequence), id="test")

        wrapper = AccessibilityStartRegion(region_start=0, region_end=20, window_size=10, step=5)
        result = wrapper(record)

        # Check that keys follow the stable pattern PU_<start>
        expected_keys = {"PU_0", "PU_5", "PU_10", "PU_15"}
        assert set(result.keys()) == expected_keys

    def test_wrapper_from_init_export(self):
        """Test that wrapper can be imported from __init__."""
        from biotooler.families.viennarna import AccessibilityStartRegion

        wrapper = AccessibilityStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        assert wrapper.region_start == 0

    def test_get_features_includes_wrapper(self):
        """Test that get_features() includes AccessibilityStartRegion."""
        from biotooler.families.viennarna import get_features
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        features = get_features()
        assert AccessibilityStartRegion in features

    def test_wrapper_values_in_valid_range(self):
        """Test that aggregated PU values are valid probabilities in [0, 1]."""
        from biotooler.families.viennarna.start_region_accessibility import (
            AccessibilityStartRegion,
        )

        sequence = "ACGUACGUACGUACGUACGUACGUACGUACGU"
        record = SeqRecord(Seq(sequence), id="test")

        wrapper = AccessibilityStartRegion(region_start=0, region_end=30, window_size=15, step=10)
        result = wrapper(record)

        # All values should be probabilities
        for key, value in result.items():
            assert 0.0 <= value <= 1.0, f"{key} has invalid value {value}"
