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
        assert str(windows[2].seq) == "ACG"  # 4-7 (partial)
        assert str(windows[3].seq) == "G"  # 6-7 (partial)

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

    def test_iter_windows_protein_step_1_ok(self):
        """Test that iter_windows works with protein sequences using step=1."""
        # Create a protein sequence
        protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein1")
        protein_record.annotations["molecule_type"] = "protein"

        # Extract overlapping windows with step=1
        windows = list(iter_windows(protein_record, window_size=5, step=1))

        # Should get 5 windows: positions 0, 1, 2, 3, 4
        assert len(windows) == 5

        # Verify each window
        assert str(windows[0].seq) == "MKALV"
        assert windows[0].annotations["start"] == 0
        assert windows[0].annotations["end"] == 5

        assert str(windows[1].seq) == "KALVS"
        assert windows[1].annotations["start"] == 1
        assert windows[1].annotations["end"] == 6

        assert str(windows[2].seq) == "ALVSW"
        assert windows[2].annotations["start"] == 2
        assert windows[2].annotations["end"] == 7

        assert str(windows[3].seq) == "LVSWG"
        assert windows[3].annotations["start"] == 3
        assert windows[3].annotations["end"] == 8

        assert str(windows[4].seq) == "VSWGR"
        assert windows[4].annotations["start"] == 4
        assert windows[4].annotations["end"] == 9

        # Verify molecule_type is copied to windows
        for window in windows:
            assert window.annotations["molecule_type"] == "protein"

    def test_iter_windows_step_not_multiple_of_3_ok_for_dna(self):
        """Test that iter_windows accepts step not multiple of 3 for DNA sequences."""
        # Create a DNA sequence (20 bases)
        dna_record = SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="dna1")
        dna_record.annotations["molecule_type"] = "DNA"

        # Use window_size=10, step=4 (not multiple of 3)
        # This should NOT raise an error in generic windowing
        windows = list(iter_windows(dna_record, window_size=10, step=4))

        # Positions: 0, 4, 8, 12, 16
        # 0-10 (len=10), 4-14 (len=10), 8-18 (len=10), 12-22 exceeds length
        # With drop_partial=True, we keep full windows only
        assert len(windows) == 3

        assert str(windows[0].seq) == "ACGTACGTAC"
        assert windows[0].annotations["start"] == 0
        assert windows[0].annotations["end"] == 10

        assert str(windows[1].seq) == "ACGTACGTAC"
        assert windows[1].annotations["start"] == 4
        assert windows[1].annotations["end"] == 14

        assert str(windows[2].seq) == "ACGTACGTAC"
        assert windows[2].annotations["start"] == 8
        assert windows[2].annotations["end"] == 18

        # Verify molecule_type is copied to windows
        for window in windows:
            assert window.annotations["molecule_type"] == "DNA"


class TestComputeWindowIndices:
    """Tests for compute_window_indices helper function."""

    def test_codon_space_basic(self):
        """Test basic codon-space window indexing."""
        from biotooler.core.windowing import compute_window_indices

        # 15 nt = 5 codons, window size 9 nt, step 3 nt
        indices = compute_window_indices(
            15, window_size=9, step=3, position_space="codon"
        )

        assert len(indices) == 3
        assert indices[0] == (0, 9)
        assert indices[1] == (3, 12)
        assert indices[2] == (6, 15)

    def test_codon_space_with_start_offset(self):
        """Test codon-space windowing with start offset."""
        from biotooler.core.windowing import compute_window_indices

        indices = compute_window_indices(
            15, window_size=9, step=3, position_space="codon", start_offset=3
        )

        assert len(indices) == 2
        assert indices[0] == (3, 12)
        assert indices[1] == (6, 15)

    def test_codon_space_step_not_multiple_of_3_error(self):
        """Test that step not multiple of 3 raises error for codon space."""
        from biotooler.core.windowing import compute_window_indices

        with pytest.raises(ValueError, match="step must be a multiple of 3"):
            compute_window_indices(15, window_size=9, step=1, position_space="codon")

        with pytest.raises(ValueError, match="step must be a multiple of 3"):
            compute_window_indices(15, window_size=9, step=4, position_space="codon")

    def test_codon_space_start_offset_not_multiple_of_3_error(self):
        """Test that start_offset not multiple of 3 raises error for codon space."""
        from biotooler.core.windowing import compute_window_indices

        with pytest.raises(ValueError, match="start_offset must be a multiple of 3"):
            compute_window_indices(
                15, window_size=9, step=3, position_space="codon", start_offset=1
            )

    def test_residue_space_step_1_allowed(self):
        """Test that step=1 is allowed for residue space."""
        from biotooler.core.windowing import compute_window_indices

        # 10 residues, window size 5, step 1
        indices = compute_window_indices(
            10, window_size=5, step=1, position_space="residue"
        )

        # Should get windows at positions 0, 1, 2, 3, 4, 5
        assert len(indices) == 6
        assert indices[0] == (0, 5)
        assert indices[1] == (1, 6)
        assert indices[2] == (2, 7)
        assert indices[3] == (3, 8)
        assert indices[4] == (4, 9)
        assert indices[5] == (5, 10)

    def test_residue_space_any_step_allowed(self):
        """Test that any positive step is allowed for residue space."""
        from biotooler.core.windowing import compute_window_indices

        # Step 2 (not multiple of 3) should work for residue space
        indices = compute_window_indices(
            10, window_size=4, step=2, position_space="residue"
        )

        assert len(indices) == 4
        assert indices[0] == (0, 4)
        assert indices[1] == (2, 6)
        assert indices[2] == (4, 8)
        assert indices[3] == (6, 10)

    def test_drop_partial_true(self):
        """Test that partial windows are dropped when drop_partial=True."""
        from biotooler.core.windowing import compute_window_indices

        # 14 nt, window 9 nt, step 3 nt
        indices = compute_window_indices(
            14, window_size=9, step=3, position_space="codon", drop_partial=True
        )

        # Should only get full windows at 0 and 3
        assert len(indices) == 2
        assert indices[0] == (0, 9)
        assert indices[1] == (3, 12)

    def test_drop_partial_false(self):
        """Test that partial windows are included when drop_partial=False."""
        from biotooler.core.windowing import compute_window_indices

        # 14 nt, window 9 nt, step 3 nt
        indices = compute_window_indices(
            14, window_size=9, step=3, position_space="codon", drop_partial=False
        )

        # Should get windows at 0, 3, 6, 9, 12
        assert len(indices) == 5
        assert indices[0] == (0, 9)  # Full
        assert indices[1] == (3, 12)  # Full
        assert indices[2] == (6, 14)  # Partial (8 nt)
        assert indices[3] == (9, 14)  # Partial (5 nt)
        assert indices[4] == (12, 14)  # Partial (2 nt)

    def test_invalid_window_size(self):
        """Test that invalid window_size raises ValueError."""
        from biotooler.core.windowing import compute_window_indices

        with pytest.raises(ValueError, match="window_size must be positive"):
            compute_window_indices(15, window_size=0, step=3, position_space="codon")

        with pytest.raises(ValueError, match="window_size must be positive"):
            compute_window_indices(15, window_size=-1, step=3, position_space="codon")

    def test_invalid_step(self):
        """Test that invalid step raises ValueError."""
        from biotooler.core.windowing import compute_window_indices

        with pytest.raises(ValueError, match="step must be positive"):
            compute_window_indices(15, window_size=9, step=0, position_space="codon")

        with pytest.raises(ValueError, match="step must be positive"):
            compute_window_indices(15, window_size=9, step=-3, position_space="codon")

    def test_invalid_position_space(self):
        """Test that invalid position_space raises ValueError."""
        from biotooler.core.windowing import compute_window_indices

        with pytest.raises(ValueError, match="position_space must be"):
            compute_window_indices(
                15, window_size=9, step=3, position_space="invalid"
            )

    def test_empty_result_when_window_too_large(self):
        """Test empty result when window size exceeds sequence length."""
        from biotooler.core.windowing import compute_window_indices

        indices = compute_window_indices(
            6, window_size=9, step=3, position_space="codon", drop_partial=True
        )

        assert len(indices) == 0

    def test_single_window_exact_fit(self):
        """Test single window when it exactly fits the sequence."""
        from biotooler.core.windowing import compute_window_indices

        indices = compute_window_indices(
            9, window_size=9, step=3, position_space="codon"
        )

        assert len(indices) == 1
        assert indices[0] == (0, 9)


class TestComputeWindowIndexArrays:
    """Tests for compute_window_index_arrays helper function."""

    def test_codon_output_space(self):
        """Test conversion to codon index space."""
        from biotooler.core.windowing import compute_window_index_arrays

        # 15 nt = 5 codons, window 9 nt (3 codons), step 3 nt (1 codon)
        arrays = compute_window_index_arrays(
            15, window_size=9, step=3, position_space="codon", output_space="codon"
        )

        assert len(arrays) == 3

        # Window 0: nt 0-9 -> codons 0-2 (indices 0, 1, 2)
        assert arrays[0].tolist() == [0, 1, 2]

        # Window 1: nt 3-12 -> codons 1-3 (indices 1, 2, 3)
        assert arrays[1].tolist() == [1, 2, 3]

        # Window 2: nt 6-15 -> codons 2-4 (indices 2, 3, 4)
        assert arrays[2].tolist() == [2, 3, 4]

    def test_residue_output_space_step_1(self):
        """Test residue space with step=1."""
        from biotooler.core.windowing import compute_window_index_arrays

        # 10 residues, window 3, step 1
        arrays = compute_window_index_arrays(
            10, window_size=3, step=1, position_space="residue", output_space="residue"
        )

        assert len(arrays) == 8

        # First few windows
        assert arrays[0].tolist() == [0, 1, 2]
        assert arrays[1].tolist() == [1, 2, 3]
        assert arrays[2].tolist() == [2, 3, 4]

        # Last window
        assert arrays[7].tolist() == [7, 8, 9]

    def test_residue_output_space_default(self):
        """Test that output_space=None keeps residue indices."""
        from biotooler.core.windowing import compute_window_index_arrays

        arrays = compute_window_index_arrays(
            10, window_size=3, step=2, position_space="residue"
        )

        assert len(arrays) == 4
        assert arrays[0].tolist() == [0, 1, 2]
        assert arrays[1].tolist() == [2, 3, 4]
        assert arrays[2].tolist() == [4, 5, 6]
        assert arrays[3].tolist() == [6, 7, 8]

    def test_with_start_offset(self):
        """Test with start_offset in codon space."""
        from biotooler.core.windowing import compute_window_index_arrays

        arrays = compute_window_index_arrays(
            15,
            window_size=9,
            step=3,
            position_space="codon",
            output_space="codon",
            start_offset=3,
        )

        assert len(arrays) == 2

        # Window at nt 3-12 -> codons 1-3
        assert arrays[0].tolist() == [1, 2, 3]

        # Window at nt 6-15 -> codons 2-4
        assert arrays[1].tolist() == [2, 3, 4]

    def test_invalid_output_space(self):
        """Test that invalid output_space raises ValueError."""
        from biotooler.core.windowing import compute_window_index_arrays

        with pytest.raises(ValueError, match="output_space must be"):
            compute_window_index_arrays(
                15,
                window_size=9,
                step=3,
                position_space="codon",
                output_space="invalid",
            )

    def test_drop_partial_affects_array_count(self):
        """Test that drop_partial affects the number of arrays."""
        from biotooler.core.windowing import compute_window_index_arrays

        arrays_with_drop = compute_window_index_arrays(
            14,
            window_size=9,
            step=3,
            position_space="codon",
            output_space="codon",
            drop_partial=True,
        )

        arrays_without_drop = compute_window_index_arrays(
            14,
            window_size=9,
            step=3,
            position_space="codon",
            output_space="codon",
            drop_partial=False,
        )

        assert len(arrays_with_drop) == 2  # Only full windows
        assert len(arrays_without_drop) == 5  # Including partial windows

    def test_array_dtypes(self):
        """Test that returned arrays are numpy arrays."""
        import numpy as np

        from biotooler.core.windowing import compute_window_index_arrays

        arrays = compute_window_index_arrays(
            15, window_size=9, step=3, position_space="codon", output_space="codon"
        )

        for arr in arrays:
            assert isinstance(arr, np.ndarray)
            assert arr.dtype == np.int_  # Should be integer type
