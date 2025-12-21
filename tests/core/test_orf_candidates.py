"""Tests for ORF candidate finding functionality."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_candidates import find_orf_candidates


class TestFindOrfCandidates:
    """Tests for find_orf_candidates function."""

    def test_simple_single_orf(self):
        """Test a simple sequence with one ORF."""
        result = find_orf_candidates("ATGAAATAA")
        assert result == [(0, 9)]

    def test_multiple_starts_and_stops(self):
        """Test sequence with multiple start and stop codons across frames."""
        # ATG AAA TAG ATG CCC TAG TAA
        # 0   3   6   9   12  15  18
        # Start at 0: stops at 6 (0->9), 15 (0->18), 18 (0->21)
        # Start at 9: stops at 15 (9->18), 18 (9->21)
        seq = "ATGAAATAGATGCCCTAGTAA"
        result = find_orf_candidates(seq)
        expected = [(0, 9), (0, 18), (0, 21), (9, 18), (9, 21)]
        assert result == expected

    def test_complex_multiple_frames(self):
        """Test complex sequence with starts and stops in different frames."""
        # Create a sequence with codons in different frames
        # Position: 0   3   6   9   12  15  18  21
        # Frame 0:  ATG GGG TAA ATG CCC TAG AAA ATG
        # Frame 1:   T  GGG GTA AAT GCC CTA GAA AAT G
        # Frame 2:    TG GGG TAA ATG CCC TAG AAA ATG
        seq = "ATGGGGTAAATGCCCTAGAAAATG"
        result = find_orf_candidates(seq)
        # Start at 0 can reach stop at 6 (0->9) and 15 (0->18)
        # Start at 9 can reach stop at 15 (9->18)
        # Start at 21 has no in-frame stop
        expected = [(0, 9), (0, 18), (9, 18)]
        assert result == expected

    def test_rna_with_u_conversion(self):
        """Test RNA sequence with U is properly converted to T."""
        # RNA sequence: AUG UGU UAA
        rna_seq = "AUGUGUUAA"
        result = find_orf_candidates(rna_seq)
        assert result == [(0, 9)]

    def test_rna_with_aug_start_and_u_stops(self):
        """Test RNA with AUG start and U-containing stops."""
        # Should work because U->T conversion happens
        result = find_orf_candidates("AUGAAAUAA", start="AUG", stops=("UAA", "UAG", "UGA"))
        assert result == [(0, 9)]

    def test_no_candidates_empty_list(self):
        """Test sequence with no valid ORF candidates returns empty list."""
        result = find_orf_candidates("AAACCCGGG")
        assert result == []

    def test_no_stop_codons(self):
        """Test sequence with start but no stop codons."""
        result = find_orf_candidates("ATGAAACCCGGG")
        assert result == []

    def test_no_start_codons(self):
        """Test sequence with stop but no start codons."""
        result = find_orf_candidates("AAACCCGGGTAA")
        assert result == []

    def test_out_of_frame_stop_not_included(self):
        """Test that out-of-frame stops are not included."""
        # ATG AA TAA (stop at position 5, not in frame with start at 0)
        # Position 5 - 0 = 5, which is not divisible by 3
        seq = "ATGAATAA"
        result = find_orf_candidates(seq)
        assert result == []

    def test_in_frame_validation(self):
        """Test that only in-frame pairs are returned."""
        # ATG GGG TAA (stop at position 6, in frame)
        # ATG G TAA at position 4 would be out of frame
        seq = "ATGGGGTAA"
        result = find_orf_candidates(seq)
        # Position 6 - 0 = 6, divisible by 3, so in frame
        assert result == [(0, 9)]

    def test_stop_before_start_not_included(self):
        """Test that stop codons before start are not included."""
        seq = "TAAATG"
        result = find_orf_candidates(seq)
        assert result == []

    def test_accepts_seq_object(self):
        """Test that function accepts Bio.Seq objects."""
        seq = Seq("ATGAAATAA")
        result = find_orf_candidates(seq)
        assert result == [(0, 9)]

    def test_accepts_seqrecord_object(self):
        """Test that function accepts Bio.SeqRecord objects."""
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        result = find_orf_candidates(record)
        assert result == [(0, 9)]

    def test_custom_start_codon(self):
        """Test using custom start codon."""
        seq = "TTGAAATAA"  # TTG as alternative start
        result = find_orf_candidates(seq, start="TTG")
        assert result == [(0, 9)]

    def test_custom_stop_codons(self):
        """Test using custom stop codons."""
        seq = "ATGAAATGA"
        result = find_orf_candidates(seq, stops=("TGA",))
        assert result == [(0, 9)]

    def test_lowercase_sequence(self):
        """Test that lowercase sequences are handled correctly."""
        result = find_orf_candidates("atgaaataa")
        assert result == [(0, 9)]

    def test_empty_sequence(self):
        """Test empty sequence returns empty list."""
        result = find_orf_candidates("")
        assert result == []

    def test_sequence_too_short_for_orf(self):
        """Test sequence shorter than minimum ORF."""
        result = find_orf_candidates("ATG")
        assert result == []

    def test_sorted_by_start_then_end(self):
        """Test that results are sorted by start position, then end position."""
        # Create a scenario where one start has multiple stops
        # ATG AAA TAG ATG CCC TAG TAA
        seq = "ATGAAATAGATGCCCTAGTAA"
        result = find_orf_candidates(seq)
        # Should be sorted: (0,9), (0,18), (0,21), (9,18), (9,21)
        assert result == [(0, 9), (0, 18), (0, 21), (9, 18), (9, 21)]
        # Verify sorting
        assert result == sorted(result)

    def test_adjacent_start_codons(self):
        """Test handling of adjacent start codons."""
        # ATG ATG TAA
        seq = "ATGATGTAA"
        result = find_orf_candidates(seq)
        # Start at 0 -> stop at 6 (0->9)
        # Start at 3 -> stop at 6 (3->9)
        assert result == [(0, 9), (3, 9)]

    def test_adjacent_stop_codons(self):
        """Test handling of adjacent stop codons."""
        # ATG TAA TAG
        # Position: 0   3   6
        # TAA at position 3, TAG at position 6
        seq = "ATGTAATAG"
        result = find_orf_candidates(seq)
        # Start at 0 can reach stop at 3 (0->6) and stop at 6 (0->9)
        assert result == [(0, 6), (0, 9)]

    def test_overlapping_orfs(self):
        """Test overlapping ORFs are all captured."""
        # ATG AAA ATG CCC TAA
        # Position: 0   3   6   9   12
        # Start at 0 -> stop at 12 (0->15)
        # Start at 6 -> stop at 12 (6->15)
        seq = "ATGAAAATGCCCTAA"
        result = find_orf_candidates(seq)
        assert result == [(0, 15), (6, 15)]

    def test_long_sequence_multiple_orfs(self):
        """Test longer sequence with multiple distinct ORFs."""
        # Create multiple distinct ORFs
        # ATG AAA TAA GGG ATG CCC TAG
        # Position: 0   3   6   9   12  15  18
        seq = "ATGAAATAAGGGATGCCCTAG"
        result = find_orf_candidates(seq)
        # Start at 0 can reach stop at 6 (0->9) and stop at 18 (0->21)
        # Start at 12 can reach stop at 18 (12->21)
        assert result == [(0, 9), (0, 21), (12, 21)]

    def test_all_stop_codons_default(self):
        """Test that all three default stop codons work."""
        # Test TAA
        assert find_orf_candidates("ATGAAATAA") == [(0, 9)]
        # Test TAG
        assert find_orf_candidates("ATGAAATAG") == [(0, 9)]
        # Test TGA
        assert find_orf_candidates("ATGAAATGA") == [(0, 9)]

    def test_return_includes_stop_codon(self):
        """Test that returned end position includes the stop codon (stop_idx + 3)."""
        # ATG AAA TAA
        # Positions: 0-2 (ATG), 3-5 (AAA), 6-8 (TAA)
        # Stop codon starts at 6, so result should be (0, 9) = (0, 6+3)
        result = find_orf_candidates("ATGAAATAA")
        assert result == [(0, 9)]
        # Verify the span includes all 9 bases
        seq = "ATGAAATAA"
        start, end = result[0]
        assert seq[start:end] == "ATGAAATAA"

    def test_rejects_protein_seqrecord(self):
        """Test that protein SeqRecord raises ValueError."""
        import pytest

        from biotooler.core.record import coerce_record

        # Create a protein SeqRecord
        protein_record = coerce_record("MKLVLS", "protein", id="test")

        # Should raise ValueError
        with pytest.raises(ValueError, match="DNA/RNA sequences"):
            find_orf_candidates(protein_record)

    def test_accepts_dna_seqrecord(self):
        """Test that DNA SeqRecord is accepted."""
        from biotooler.core.record import coerce_record

        # Create a DNA SeqRecord
        dna_record = coerce_record("ATGAAATAA", "DNA", id="test")
        result = find_orf_candidates(dna_record)
        assert result == [(0, 9)]

    def test_accepts_rna_seqrecord(self):
        """Test that RNA SeqRecord is accepted."""
        from biotooler.core.record import coerce_record

        # Create an RNA SeqRecord
        rna_record = coerce_record("AUGUGUUAA", "RNA", id="test")
        result = find_orf_candidates(rna_record)
        assert result == [(0, 9)]

    def test_accepts_seqrecord_without_molecule_type(self):
        """Test that SeqRecord without molecule_type annotation is accepted."""
        # SeqRecord without molecule_type should work fine
        record = SeqRecord(Seq("ATGAAATAA"), id="test")
        result = find_orf_candidates(record)
        assert result == [(0, 9)]
