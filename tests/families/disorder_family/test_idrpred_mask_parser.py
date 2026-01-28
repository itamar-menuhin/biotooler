"""Tests for IDRPred TSV parsing and mask conversion.

These tests verify the robustness of TSV parsing and region-to-mask conversion
without requiring the idrpred CLI to be installed.
"""

import numpy as np


class TestParseTsvRegions:
    """Test parse_tsv_regions function."""

    def test_parse_basic_tsv(self):
        """Test parsing basic TSV with IDR regions."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "1\t10\n20\t30\n"
        regions = parse_tsv_regions(tsv_output)

        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_empty_tsv(self):
        """Test parsing empty TSV output (no IDRs predicted)."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = ""
        regions = parse_tsv_regions(tsv_output)

        assert len(regions) == 0
        assert regions == []

    def test_parse_tsv_with_empty_lines(self):
        """Test parsing TSV with empty lines."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "1\t10\n\n20\t30\n\n"
        regions = parse_tsv_regions(tsv_output)

        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_tsv_with_comment_lines(self):
        """Test parsing TSV with comment/header lines."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "# IDRPred output\n# Start\tEnd\n1\t10\n20\t30\n"
        regions = parse_tsv_regions(tsv_output)

        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_tsv_with_extra_columns(self):
        """Test parsing TSV with extra columns beyond start/end."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        # IDRPred might output additional columns (e.g., confidence, method)
        tsv_output = "1\t10\t0.95\tConsensus\n20\t30\t0.87\tVSL2b\n"
        regions = parse_tsv_regions(tsv_output)

        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_tsv_ignores_invalid_lines(self):
        """Test that parser ignores lines with invalid data."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "1\t10\ninvalid\tdata\n20\t30\nabc\txyz\n"
        regions = parse_tsv_regions(tsv_output)

        # Only valid lines should be parsed
        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_tsv_ignores_single_column_lines(self):
        """Test that parser ignores lines with only one column."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "1\t10\n20\n30\t40\n"
        regions = parse_tsv_regions(tsv_output)

        # Single column lines should be ignored
        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (30, 40)

    def test_parse_tsv_ignores_zero_or_negative_coordinates(self):
        """Test that parser ignores regions with zero or negative coordinates."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "0\t10\n-5\t5\n1\t10\n20\t30\n"
        regions = parse_tsv_regions(tsv_output)

        # Only regions with positive coordinates should be kept
        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_tsv_ignores_inverted_regions(self):
        """Test that parser ignores regions where start > end."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "10\t1\n1\t10\n30\t20\n20\t30\n"
        regions = parse_tsv_regions(tsv_output)

        # Only valid regions (start <= end) should be kept
        assert len(regions) == 2
        assert regions[0] == (1, 10)
        assert regions[1] == (20, 30)

    def test_parse_tsv_determinism(self):
        """Test that parsing is deterministic."""
        from biotooler.families.disorder.idrpred_feature import parse_tsv_regions

        tsv_output = "1\t10\n20\t30\n40\t50\n"

        # Parse multiple times
        regions1 = parse_tsv_regions(tsv_output)
        regions2 = parse_tsv_regions(tsv_output)
        regions3 = parse_tsv_regions(tsv_output)

        # Should be identical
        assert regions1 == regions2 == regions3


class TestRegionsToMask:
    """Test regions_to_mask function."""

    def test_empty_regions_returns_all_zeros(self):
        """Test that empty regions list returns all zeros."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        mask = regions_to_mask([], seq_length=10)

        assert len(mask) == 10
        assert np.all(mask == 0.0)
        assert mask.dtype == np.float64

    def test_single_region(self):
        """Test converting a single region to mask."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Region [1, 5] in 1-based inclusive coordinates
        # Should mark positions 0-4 (0-based) as IDR
        mask = regions_to_mask([(1, 5)], seq_length=10)

        expected = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_multiple_non_overlapping_regions(self):
        """Test converting multiple non-overlapping regions."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Regions [1, 3] and [6, 8] in 1-based inclusive
        mask = regions_to_mask([(1, 3), (6, 8)], seq_length=10)

        expected = np.array([1, 1, 1, 0, 0, 1, 1, 1, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_overlapping_regions(self):
        """Test that overlapping regions are merged."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Overlapping regions [1, 5] and [3, 7]
        mask = regions_to_mask([(1, 5), (3, 7)], seq_length=10)

        # Should merge into one contiguous region [1, 7]
        expected = np.array([1, 1, 1, 1, 1, 1, 1, 0, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_adjacent_regions(self):
        """Test adjacent regions."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Adjacent regions [1, 3] and [4, 6]
        mask = regions_to_mask([(1, 3), (4, 6)], seq_length=10)

        expected = np.array([1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_coordinate_clamping_below_range(self):
        """Test that coordinates below 1 are clamped."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Region with start < 1 should be clamped to 1
        mask = regions_to_mask([(-5, 3)], seq_length=10)

        # Should be clamped to [1, 3]
        expected = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_coordinate_clamping_above_range(self):
        """Test that coordinates above seq_length are clamped."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Region with end > seq_length should be clamped
        mask = regions_to_mask([(8, 15)], seq_length=10)

        # Should be clamped to [8, 10]
        expected = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1, 1], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_coordinate_clamping_both_sides(self):
        """Test clamping on both sides."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Region spanning beyond both ends
        mask = regions_to_mask([(-5, 15)], seq_length=10)

        # Should be clamped to [1, 10] = entire sequence
        expected = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_entire_sequence_as_idr(self):
        """Test marking entire sequence as IDR."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        mask = regions_to_mask([(1, 10)], seq_length=10)

        expected = np.ones(10, dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_single_residue_region(self):
        """Test single residue region."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        # Single residue at position 5 (1-based)
        mask = regions_to_mask([(5, 5)], seq_length=10)

        expected = np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_mask_dtype_is_float64(self):
        """Test that mask has float64 dtype."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        mask = regions_to_mask([(1, 5)], seq_length=10)

        assert mask.dtype == np.float64

    def test_mask_shape_matches_sequence_length(self):
        """Test that mask shape matches sequence length."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        for seq_length in [5, 10, 50, 100]:
            mask = regions_to_mask([(1, 3)], seq_length=seq_length)
            assert mask.shape == (seq_length,)

    def test_determinism(self):
        """Test that mask generation is deterministic."""
        from biotooler.families.disorder.idrpred_feature import regions_to_mask

        regions = [(1, 5), (10, 15), (20, 25)]

        # Generate mask multiple times
        mask1 = regions_to_mask(regions, seq_length=30)
        mask2 = regions_to_mask(regions, seq_length=30)
        mask3 = regions_to_mask(regions, seq_length=30)

        # Should be identical
        np.testing.assert_array_equal(mask1, mask2)
        np.testing.assert_array_equal(mask2, mask3)


class TestIntegration:
    """Test integration of parsing and mask conversion."""

    def test_parse_and_convert_basic(self):
        """Test end-to-end: parse TSV and convert to mask."""
        from biotooler.families.disorder.idrpred_feature import (
            parse_tsv_regions,
            regions_to_mask,
        )

        tsv_output = "1\t5\n10\t12\n"
        regions = parse_tsv_regions(tsv_output)
        mask = regions_to_mask(regions, seq_length=15)

        expected = np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_parse_and_convert_no_idrs(self):
        """Test when no IDRs are predicted."""
        from biotooler.families.disorder.idrpred_feature import (
            parse_tsv_regions,
            regions_to_mask,
        )

        tsv_output = ""
        regions = parse_tsv_regions(tsv_output)
        mask = regions_to_mask(regions, seq_length=10)

        expected = np.zeros(10, dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)

    def test_parse_and_convert_with_comments(self):
        """Test with TSV containing comments."""
        from biotooler.families.disorder.idrpred_feature import (
            parse_tsv_regions,
            regions_to_mask,
        )

        tsv_output = "# Header\n# Start\tEnd\n1\t3\n7\t9\n"
        regions = parse_tsv_regions(tsv_output)
        mask = regions_to_mask(regions, seq_length=10)

        expected = np.array([1, 1, 1, 0, 0, 0, 1, 1, 1, 0], dtype=np.float64)
        np.testing.assert_array_equal(mask, expected)
