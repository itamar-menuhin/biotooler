"""Tests for ORF store utilities."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_store import OrfSpan, attach_orf, get_orf, select_orf_by_index
from biotooler.core.record import coerce_record


class TestSelectOrfByIndex:
    """Tests for select_orf_by_index function."""

    def test_select_first_orf(self):
        """Test selecting the first ORF from candidates."""
        candidates = [(0, 9), (9, 18), (18, 27)]
        result = select_orf_by_index(candidates, 0)
        assert result == (0, 9)

    def test_select_middle_orf(self):
        """Test selecting a middle ORF from candidates."""
        candidates = [(0, 9), (9, 18), (18, 27)]
        result = select_orf_by_index(candidates, 1)
        assert result == (9, 18)

    def test_select_last_orf(self):
        """Test selecting the last ORF from candidates."""
        candidates = [(0, 9), (9, 18), (18, 27)]
        result = select_orf_by_index(candidates, 2)
        assert result == (18, 27)

    def test_single_candidate(self):
        """Test selecting from a single candidate."""
        candidates = [(0, 9)]
        result = select_orf_by_index(candidates, 0)
        assert result == (0, 9)

    def test_negative_index_raises_error(self):
        """Test that negative indices raise clear error."""
        candidates = [(0, 9), (9, 18)]
        with pytest.raises(IndexError, match=r"ORF index -1 is out of range.*valid indices: 0-1"):
            select_orf_by_index(candidates, -1)

    def test_index_too_large_raises_error(self):
        """Test that index >= len(candidates) raises clear error."""
        candidates = [(0, 9), (9, 18)]
        with pytest.raises(IndexError, match=r"ORF index 2 is out of range.*valid indices: 0-1"):
            select_orf_by_index(candidates, 2)

    def test_way_too_large_index_raises_error(self):
        """Test that very large index raises clear error."""
        candidates = [(0, 9), (9, 18)]
        with pytest.raises(IndexError, match=r"ORF index 100 is out of range.*valid indices: 0-1"):
            select_orf_by_index(candidates, 100)

    def test_empty_candidates_raises_error(self):
        """Test that empty candidate list raises clear error."""
        candidates = []
        with pytest.raises(
            IndexError, match=r"ORF index 0 is out of range: no candidates available"
        ):
            select_orf_by_index(candidates, 0)

    def test_empty_candidates_with_nonzero_index_raises_error(self):
        """Test that empty candidate list with non-zero index raises clear error."""
        candidates = []
        with pytest.raises(
            IndexError, match=r"ORF index 5 is out of range: no candidates available"
        ):
            select_orf_by_index(candidates, 5)

    def test_type_annotation_orfspan(self):
        """Test that OrfSpan type alias is defined and usable."""
        # This is mainly for type checkers, but we can verify it exists
        assert OrfSpan is not None
        # Verify it's a tuple type hint
        orf: OrfSpan = (0, 9)
        assert orf == (0, 9)


class TestAttachOrfAndGetOrf:
    """Tests for attach_orf and get_orf functions."""

    def test_attach_and_get_roundtrip(self):
        """Test that attaching and getting an ORF works correctly."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        orf = (0, 9)

        # Attach the ORF
        result = attach_orf(record, orf)

        # Verify the record is returned
        assert result is record

        # Retrieve the ORF
        retrieved_orf = get_orf(record)
        assert retrieved_orf == orf

    def test_attach_with_custom_key(self):
        """Test attaching ORF with a custom key."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        orf = (3, 12)

        attach_orf(record, orf, key="my.custom.key")

        # Should not be found with default key
        with pytest.raises(KeyError):
            get_orf(record)

        # Should be found with custom key
        retrieved_orf = get_orf(record, key="my.custom.key")
        assert retrieved_orf == orf

    def test_get_orf_missing_raises_clear_error(self):
        """Test that get_orf raises clear error when ORF is missing."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test_record")

        with pytest.raises(
            KeyError,
            match=(
                r"ORF not found in record 'test_record': no annotation with key "
                r"'biotooler\.orf'\. Use attach_orf\(\) to store an ORF first\."
            ),
        ):
            get_orf(record)

    def test_get_orf_missing_custom_key_raises_clear_error(self):
        """Test that get_orf with custom key raises clear error when missing."""
        record = SeqRecord(Seq("ATGAAATAA"), id="my_seq")

        with pytest.raises(
            KeyError,
            match=r"ORF not found in record 'my_seq': no annotation with key 'custom\.key'\. "
                  r"Use attach_orf\(\) to store an ORF first\."
        ):
            get_orf(record, key="custom.key")

    def test_attach_multiple_orfs_with_different_keys(self):
        """Test attaching multiple ORFs with different keys."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        orf1 = (0, 9)
        orf2 = (3, 12)

        attach_orf(record, orf1, key="orf1")
        attach_orf(record, orf2, key="orf2")

        assert get_orf(record, key="orf1") == orf1
        assert get_orf(record, key="orf2") == orf2

    def test_attach_overwrites_existing_orf(self):
        """Test that attaching a new ORF overwrites the existing one."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        orf1 = (0, 9)
        orf2 = (3, 12)

        attach_orf(record, orf1)
        attach_orf(record, orf2)

        # Should get the second ORF
        assert get_orf(record) == orf2

    def test_attach_with_dna_record(self):
        """Test attaching ORF to DNA SeqRecord works."""
        record = coerce_record("ATGAAATAA", "DNA", id="dna_test")
        orf = (0, 9)

        attach_orf(record, orf)
        assert get_orf(record) == orf

    def test_attach_with_rna_record(self):
        """Test attaching ORF to RNA SeqRecord works."""
        record = coerce_record("AUGUGUUAA", "RNA", id="rna_test")
        orf = (0, 9)

        attach_orf(record, orf)
        assert get_orf(record) == orf

    def test_attach_to_protein_record_raises_error(self):
        """Test that attaching ORF to protein SeqRecord raises ValueError."""
        record = coerce_record("MKLVLS", "protein", id="protein_test")
        orf = (0, 6)

        with pytest.raises(
            ValueError,
            match=r"ORF operations are only supported for DNA/RNA sequences, not protein sequences"
        ):
            attach_orf(record, orf)

    def test_attach_without_molecule_type_annotation(self):
        """Test attaching ORF to SeqRecord without molecule_type annotation works."""
        # SeqRecord without molecule_type should work fine
        record = SeqRecord(Seq("ATGAAATAA"), id="no_mol_type")
        orf = (0, 9)

        attach_orf(record, orf)
        assert get_orf(record) == orf

    def test_orf_coordinates_are_preserved(self):
        """Test that ORF coordinates are preserved exactly."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        # Test various coordinate ranges
        test_cases = [
            (0, 9),
            (5, 15),
            (100, 200),
            (0, 0),  # Empty span
        ]

        for orf in test_cases:
            attach_orf(record, orf)
            retrieved = get_orf(record)
            assert retrieved == orf

    def test_attach_returns_same_record_instance(self):
        """Test that attach_orf returns the same record instance."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        orf = (0, 9)

        result = attach_orf(record, orf)

        # Should be the exact same object
        assert result is record
        # And the annotation should be on the original record
        assert "biotooler.orf" in record.annotations
