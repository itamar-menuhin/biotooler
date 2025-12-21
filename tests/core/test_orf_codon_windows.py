"""Tests for iter_orf_codon_windows function."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.windowing import iter_orf_codon_windows


class TestIterOrfCodonWindows:
    """Tests for iter_orf_codon_windows function."""

    def test_basic_codon_windowing(self):
        """Test basic ORF windowing with codon-aligned step."""
        # 15 nt = 5 codons
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        # Window size 9nt (3 codons), step 3nt (1 codon)
        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=3)
        )

        # Should get windows at positions 0, 3, 6
        assert len(windows) == 3
        assert str(windows[0].seq) == "ATGAAACCC"
        assert str(windows[1].seq) == "AAACCCGGG"
        assert str(windows[2].seq) == "CCCGGGTTT"

    def test_window_start_annotations(self):
        """Test that window_start is relative to ORF start."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=3)
        )

        assert windows[0].annotations["window_start"] == 0
        assert windows[1].annotations["window_start"] == 3
        assert windows[2].annotations["window_start"] == 6

    def test_window_end_annotations(self):
        """Test that window_end is relative to ORF start."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=3)
        )

        assert windows[0].annotations["window_end"] == 9
        assert windows[1].annotations["window_end"] == 12
        assert windows[2].annotations["window_end"] == 15

    def test_window_index_annotations(self):
        """Test that window_index starts at 0 and increments."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=3)
        )

        assert windows[0].annotations["window_index"] == 0
        assert windows[1].annotations["window_index"] == 1
        assert windows[2].annotations["window_index"] == 2

    def test_orf_coordinates_in_annotations(self):
        """Test that ORF coordinates are included in window annotations."""
        record = SeqRecord(Seq("NNNNATGAAACCCGGGTTTNNNN"), id="seq1")

        windows = list(
            iter_orf_codon_windows(record, orf=(4, 19), window_nt=9, step_nt=3)
        )

        for window in windows:
            assert window.annotations["orf_start"] == 4
            assert window.annotations["orf_end"] == 19

    def test_absolute_positions_with_offset_orf(self):
        """Test absolute positions when ORF doesn't start at 0."""
        record = SeqRecord(Seq("NNNNATGAAACCCGGGTTTNNNN"), id="seq1")

        windows = list(
            iter_orf_codon_windows(record, orf=(4, 19), window_nt=9, step_nt=3)
        )

        # Window 0: starts at 4+0=4, ends at 4+9=13
        assert windows[0].annotations["start"] == 4
        assert windows[0].annotations["end"] == 13

        # Window 1: starts at 4+3=7, ends at 4+12=16
        assert windows[1].annotations["start"] == 7
        assert windows[1].annotations["end"] == 16

    def test_step_not_multiple_of_3_raises_error(self):
        """Test that step_nt not multiple of 3 raises ValueError."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(ValueError, match="step_nt must be a multiple of 3"):
            list(iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=4))

    def test_step_1_raises_error(self):
        """Test that step_nt=1 raises ValueError."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(ValueError, match="step_nt must be a multiple of 3"):
            list(iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=1))

    def test_step_6_accepted(self):
        """Test that step_nt=6 (multiple of 3) is accepted."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        # Should work without error
        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=6)
        )
        assert len(windows) == 2
        assert windows[0].annotations["window_start"] == 0
        assert windows[1].annotations["window_start"] == 6

    def test_drop_partial_true(self):
        """Test that partial windows are dropped when drop_partial=True."""
        # 14 nt, windows of 9 nt with step 3 nt
        record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

        windows = list(
            iter_orf_codon_windows(
                record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=True
            )
        )

        # Should only get windows at 0 and 3 (full size)
        assert len(windows) == 2
        assert len(windows[0].seq) == 9
        assert len(windows[1].seq) == 9

    def test_drop_partial_false(self):
        """Test that partial windows are included when drop_partial=False."""
        record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

        windows = list(
            iter_orf_codon_windows(
                record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=False
            )
        )

        # Should get windows at 0, 3, 6, 9, 12
        assert len(windows) == 5
        assert len(windows[0].seq) == 9  # Full
        assert len(windows[1].seq) == 9  # Full
        assert len(windows[2].seq) == 8  # Partial
        assert len(windows[3].seq) == 5  # Partial
        assert len(windows[4].seq) == 2  # Partial

    def test_window_size_larger_than_orf(self):
        """Test behavior when window is larger than ORF."""
        record = SeqRecord(Seq("ATGAAA"), id="seq1")

        # Window size 12 > ORF size 6
        windows = list(
            iter_orf_codon_windows(
                record, orf=(0, 6), window_nt=12, step_nt=3, drop_partial=True
            )
        )

        # No full windows fit
        assert len(windows) == 0

    def test_window_size_larger_than_orf_drop_partial_false(self):
        """Test partial window when window_nt > ORF and drop_partial=False."""
        record = SeqRecord(Seq("ATGAAA"), id="seq1")

        windows = list(
            iter_orf_codon_windows(
                record, orf=(0, 6), window_nt=12, step_nt=3, drop_partial=False
            )
        )

        # With step=3, we get windows at 0 and 3 (both partial)
        assert len(windows) == 2
        assert str(windows[0].seq) == "ATGAAA"  # Window at 0, size 6
        assert str(windows[1].seq) == "AAA"  # Window at 3, size 3

    def test_parent_annotations_copied(self):
        """Test that parent annotations are copied to windows."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
        record.annotations["molecule_type"] = "DNA"
        record.annotations["custom"] = "value"

        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=3)
        )

        for window in windows:
            assert window.annotations["molecule_type"] == "DNA"
            assert window.annotations["custom"] == "value"

    def test_window_id_format(self):
        """Test window ID format includes window_index."""
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test_seq")

        windows = list(
            iter_orf_codon_windows(record, orf=(0, 15), window_nt=9, step_nt=3)
        )

        assert windows[0].id == "test_seq_orf_window_0"
        assert windows[1].id == "test_seq_orf_window_1"
        assert windows[2].id == "test_seq_orf_window_2"

    def test_empty_orf(self):
        """Test behavior with empty ORF."""
        record = SeqRecord(Seq("ATGAAACCC"), id="seq1")

        windows = list(
            iter_orf_codon_windows(record, orf=(3, 3), window_nt=9, step_nt=3)
        )

        assert len(windows) == 0

    def test_single_codon_orf(self):
        """Test with ORF containing only one codon."""
        record = SeqRecord(Seq("ATGNNNNNN"), id="seq1")

        # Window size 3, step 3
        windows = list(
            iter_orf_codon_windows(record, orf=(0, 3), window_nt=3, step_nt=3)
        )

        assert len(windows) == 1
        assert str(windows[0].seq) == "ATG"
        assert windows[0].annotations["window_start"] == 0

    def test_negative_window_size_raises_error(self):
        """Test that negative window_nt raises ValueError."""
        record = SeqRecord(Seq("ATGAAACCC"), id="seq1")

        with pytest.raises(ValueError, match="window_nt must be positive"):
            list(iter_orf_codon_windows(record, orf=(0, 9), window_nt=-3, step_nt=3))

    def test_zero_window_size_raises_error(self):
        """Test that zero window_nt raises ValueError."""
        record = SeqRecord(Seq("ATGAAACCC"), id="seq1")

        with pytest.raises(ValueError, match="window_nt must be positive"):
            list(iter_orf_codon_windows(record, orf=(0, 9), window_nt=0, step_nt=3))

    def test_negative_step_raises_error(self):
        """Test that negative step_nt raises ValueError."""
        record = SeqRecord(Seq("ATGAAACCC"), id="seq1")

        with pytest.raises(ValueError, match="step_nt must be positive"):
            list(iter_orf_codon_windows(record, orf=(0, 9), window_nt=3, step_nt=-3))

    def test_zero_step_raises_error(self):
        """Test that zero step_nt raises ValueError."""
        record = SeqRecord(Seq("ATGAAACCC"), id="seq1")

        with pytest.raises(ValueError, match="step_nt must be positive"):
            list(iter_orf_codon_windows(record, orf=(0, 9), window_nt=3, step_nt=0))

    def test_returns_iterator(self):
        """Test that iter_orf_codon_windows returns an iterator."""
        record = SeqRecord(Seq("ATGAAACCC"), id="seq1")

        result = iter_orf_codon_windows(record, orf=(0, 9), window_nt=6, step_nt=3)

        # Should be an iterator
        assert hasattr(result, "__iter__")
        assert hasattr(result, "__next__")
