"""Tests for sliding window utilities."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.windowing import iter_windows


class TestIterWindows:
    """Tests for iter_windows function."""

    def test_exact_fit_no_overlap(self):
        """Test windows with exact fit and no overlap."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=4))

        assert len(windows) == 2
        assert str(windows[0].seq) == "ACGT"
        assert str(windows[1].seq) == "ACGT"

    def test_exact_fit_with_overlap(self):
        """Test windows with overlap."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=2))

        assert len(windows) == 3
        assert str(windows[0].seq) == "ACGT"
        assert str(windows[1].seq) == "GTAC"
        assert str(windows[2].seq) == "ACGT"

    def test_drop_partial_true(self):
        """Test that partial windows are dropped when drop_partial=True."""
        record = SeqRecord(Seq("ACGTACG"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=2, drop_partial=True))

        # Should get windows at positions 0 and 2, but not 4 (partial)
        assert len(windows) == 2
        assert str(windows[0].seq) == "ACGT"
        assert str(windows[1].seq) == "GTAC"

    def test_drop_partial_false(self):
        """Test that partial windows are included when drop_partial=False."""
        record = SeqRecord(Seq("ACGTACG"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=2, drop_partial=False))

        # Should get windows at positions 0, 2, 4, and 6
        assert len(windows) == 4
        assert str(windows[0].seq) == "ACGT"  # 0-4
        assert str(windows[1].seq) == "GTAC"  # 2-6
        assert str(windows[2].seq) == "ACG"   # 4-7 (partial)
        assert str(windows[3].seq) == "G"     # 6-7 (partial)

    def test_window_metadata_parent_id(self):
        """Test that parent_id is set correctly in window annotations."""
        record = SeqRecord(Seq("ACGTACGT"), id="parent_seq")
        windows = list(iter_windows(record, window_size=4, step=4))

        for window in windows:
            assert window.annotations["parent_id"] == "parent_seq"

    def test_window_metadata_positions(self):
        """Test that start and end positions are set correctly."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=2))

        assert windows[0].annotations["start"] == 0
        assert windows[0].annotations["end"] == 4

        assert windows[1].annotations["start"] == 2
        assert windows[1].annotations["end"] == 6

        assert windows[2].annotations["start"] == 4
        assert windows[2].annotations["end"] == 8

    def test_window_metadata_index(self):
        """Test that window_index is set correctly."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=2))

        assert windows[0].annotations["window_index"] == 0
        assert windows[1].annotations["window_index"] == 1
        assert windows[2].annotations["window_index"] == 2

    def test_window_id_format(self):
        """Test that window IDs are formatted correctly."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=4))

        assert windows[0].id == "seq1_window_0"
        assert windows[1].id == "seq1_window_1"

    def test_window_description_format(self):
        """Test that window descriptions are formatted correctly."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=4))

        assert windows[0].description == "Window 0 from seq1 [0:4]"
        assert windows[1].description == "Window 1 from seq1 [4:8]"

    def test_parent_annotations_copied(self):
        """Test that parent annotations are copied to windows."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["custom_field"] = "custom_value"

        windows = list(iter_windows(record, window_size=4, step=4))

        for window in windows:
            assert window.annotations["molecule_type"] == "DNA"
            assert window.annotations["custom_field"] == "custom_value"

    def test_window_size_larger_than_sequence(self):
        """Test behavior when window_size is larger than sequence."""
        record = SeqRecord(Seq("ACG"), id="seq1")
        windows = list(iter_windows(record, window_size=10, step=5, drop_partial=True))

        # Should get no windows since the only window would be partial
        assert len(windows) == 0

    def test_window_size_larger_than_sequence_drop_partial_false(self):
        """Test that single partial window is returned when drop_partial=False."""
        record = SeqRecord(Seq("ACG"), id="seq1")
        windows = list(iter_windows(record, window_size=10, step=5, drop_partial=False))

        assert len(windows) == 1
        assert str(windows[0].seq) == "ACG"

    def test_step_larger_than_window_size(self):
        """Test with step size larger than window size (gaps between windows)."""
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=2, step=4))

        assert len(windows) == 3
        assert str(windows[0].seq) == "AC"  # 0:2
        assert str(windows[1].seq) == "AC"  # 4:6
        assert str(windows[2].seq) == "AC"  # 8:10

    def test_window_size_one(self):
        """Test with window_size of 1."""
        record = SeqRecord(Seq("ACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=1, step=1))

        assert len(windows) == 4
        assert str(windows[0].seq) == "A"
        assert str(windows[1].seq) == "C"
        assert str(windows[2].seq) == "G"
        assert str(windows[3].seq) == "T"

    def test_step_one_with_larger_window(self):
        """Test with step=1 and larger window (maximum overlap)."""
        record = SeqRecord(Seq("ACGT"), id="seq1")
        windows = list(iter_windows(record, window_size=3, step=1))

        assert len(windows) == 2
        assert str(windows[0].seq) == "ACG"
        assert str(windows[1].seq) == "CGT"

    def test_invalid_window_size_zero(self):
        """Test that window_size of 0 raises ValueError."""
        record = SeqRecord(Seq("ACGT"), id="seq1")
        with pytest.raises(ValueError, match="window_size must be positive"):
            list(iter_windows(record, window_size=0, step=1))

    def test_invalid_window_size_negative(self):
        """Test that negative window_size raises ValueError."""
        record = SeqRecord(Seq("ACGT"), id="seq1")
        with pytest.raises(ValueError, match="window_size must be positive"):
            list(iter_windows(record, window_size=-1, step=1))

    def test_invalid_step_zero(self):
        """Test that step of 0 raises ValueError."""
        record = SeqRecord(Seq("ACGT"), id="seq1")
        with pytest.raises(ValueError, match="step must be positive"):
            list(iter_windows(record, window_size=2, step=0))

    def test_invalid_step_negative(self):
        """Test that negative step raises ValueError."""
        record = SeqRecord(Seq("ACGT"), id="seq1")
        with pytest.raises(ValueError, match="step must be positive"):
            list(iter_windows(record, window_size=2, step=-1))

    def test_empty_sequence(self):
        """Test with empty sequence."""
        record = SeqRecord(Seq(""), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=2))

        assert len(windows) == 0

    def test_sequence_normalization(self):
        """Test that sequence normalization from get_seq_str is used."""
        # Create record with lowercase and whitespace
        record = SeqRecord(Seq("ac gt"), id="seq1")
        windows = list(iter_windows(record, window_size=4, step=4))

        # Should be normalized to uppercase with whitespace removed
        assert len(windows) == 1
        assert str(windows[0].seq) == "ACGT"

    def test_caching_efficiency(self):
        """Test that get_seq_str caching is used (sequence extracted once)."""
        record = SeqRecord(Seq("acgtacgt"), id="seq1")

        # First call to iter_windows
        list(iter_windows(record, window_size=4, step=2))

        # Check that cache was created
        assert "_biotooler_seq_str" in record.annotations
        assert record.annotations["_biotooler_seq_str"] == "ACGTACGT"

        # Modify the cached value
        record.annotations["_biotooler_seq_str"] = "MODIFIED!"

        # Second call should use cached (modified) value
        windows2 = list(iter_windows(record, window_size=4, step=4))

        # Should use the modified cached value
        assert str(windows2[0].seq) == "MODI"

    def test_large_sequence_with_small_windows(self):
        """Test with a large sequence and small windows."""
        # Create a long sequence
        seq = "ACGT" * 100  # 400 characters
        record = SeqRecord(Seq(seq), id="seq1")
        windows = list(iter_windows(record, window_size=10, step=10))

        assert len(windows) == 40
        assert str(windows[0].seq) == "ACGTACGTAC"
        # Last window starts at position 390 (39*10)
        assert str(windows[-1].seq) == "GTACGTACGT"

    def test_returns_iterator(self):
        """Test that iter_windows returns an iterator."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        result = iter_windows(record, window_size=4, step=2)

        # Should be an iterator, not a list
        assert hasattr(result, "__iter__")
        assert hasattr(result, "__next__")

    def test_iterator_can_be_consumed_once(self):
        """Test that iterator can only be consumed once."""
        record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        iterator = iter_windows(record, window_size=4, step=2)

        # First consumption
        windows1 = list(iterator)
        assert len(windows1) == 3

        # Second consumption should yield nothing
        windows2 = list(iterator)
        assert len(windows2) == 0
