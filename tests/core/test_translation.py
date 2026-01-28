"""Tests for translation utilities."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_store import attach_orf
from biotooler.core.translation import ensure_protein_record


class TestEnsureProteinRecordProteinInput:
    """Tests for ensure_protein_record with protein input."""

    def test_protein_passthrough(self):
        """Test that protein sequences pass through unchanged."""
        protein = SeqRecord(Seq("MKALV"), id="test", description="test protein")
        protein.annotations["molecule_type"] = "protein"

        result = ensure_protein_record(protein)

        assert str(result.seq) == "MKALV"
        assert result.id == "test"
        assert result.description == "test protein"
        assert result.annotations["molecule_type"] == "protein"

    def test_protein_with_stop(self):
        """Test protein sequence with stop codon."""
        protein = SeqRecord(Seq("MKALV*"), id="test")
        protein.annotations["molecule_type"] = "protein"

        result = ensure_protein_record(protein)

        assert str(result.seq) == "MKALV*"


class TestEnsureProteinRecordDNATranslation:
    """Tests for ensure_protein_record with DNA input."""

    def test_basic_dna_translation(self):
        """Test basic DNA to protein translation."""
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna)

        assert str(result.seq) == "MKALV"
        assert result.annotations["molecule_type"] == "protein"
        assert result.annotations["translation_performed"] is True
        assert result.annotations["translation_table"] == 1
        assert result.annotations["translation_source"] == "DNA"

    def test_dna_with_stop_codon(self):
        """Test DNA translation with stop codon."""
        # ATGAAAGCCCTGGTGTAA = MKALV* (TAA is stop)
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTGTAA"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        # By default, strip_terminal_stop=True
        result = ensure_protein_record(dna)
        assert str(result.seq) == "MKALV"

        # With strip_terminal_stop=False
        result = ensure_protein_record(dna, strip_terminal_stop=False)
        assert str(result.seq) == "MKALV*"

    def test_dna_internal_stop_error(self):
        """Test that internal stop codons raise error by default."""
        # ATGTAAAAAGCCCTG = M*KAL (TAA is internal stop)
        dna = SeqRecord(Seq("ATGTAAAAAGCCCTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        with pytest.raises(ValueError, match="Internal stop codon"):
            ensure_protein_record(dna)

    def test_dna_internal_stop_ignore(self):
        """Test that internal stops can be ignored."""
        # ATGTAAAAAGCCCTG = M*KAL (TAA is internal stop)
        dna = SeqRecord(Seq("ATGTAAAAAGCCCTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, on_internal_stop="ignore")
        assert str(result.seq) == "M*KAL"

    def test_dna_different_genetic_code(self):
        """Test translation with different genetic code table."""
        # TGA is normally stop, but codes for Trp (W) in table 2 (vertebrate mitochondrial)
        dna = SeqRecord(Seq("ATGTGAGCC"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        # Standard code (table 1): TGA = stop
        result1 = ensure_protein_record(dna, table=1, on_internal_stop="ignore")
        assert str(result1.seq) == "M*A"

        # Vertebrate mitochondrial code (table 2): TGA = Trp
        result2 = ensure_protein_record(dna, table=2)
        assert str(result2.seq) == "MWA"
        assert result2.annotations["translation_table"] == 2


class TestEnsureProteinRecordRNATranslation:
    """Tests for ensure_protein_record with RNA input."""

    def test_basic_rna_translation(self):
        """Test basic RNA to protein translation."""
        rna = SeqRecord(Seq("AUGAAAGCCCUGGUG"), id="test")
        rna.annotations["molecule_type"] = "RNA"

        result = ensure_protein_record(rna)

        assert str(result.seq) == "MKALV"
        assert result.annotations["molecule_type"] == "protein"
        assert result.annotations["translation_performed"] is True
        assert result.annotations["translation_source"] == "RNA"

    def test_rna_with_stop(self):
        """Test RNA translation with stop codon."""
        # AUGAAAGCCCUGGUGUAA = MKALV* (UAA is stop)
        rna = SeqRecord(Seq("AUGAAAGCCCUGGUGUAA"), id="test")
        rna.annotations["molecule_type"] = "RNA"

        result = ensure_protein_record(rna)
        assert str(result.seq) == "MKALV"

        result = ensure_protein_record(rna, strip_terminal_stop=False)
        assert str(result.seq) == "MKALV*"


class TestEnsureProteinRecordORFHandling:
    """Tests for ORF handling in ensure_protein_record."""

    def test_explicit_orf(self):
        """Test translation with explicit ORF."""
        # Full sequence: NNNNATGAAAGCCNNN
        # ORF region [4:13] = ATGAAAGCC (9 bases, 3 codons) -> MKA
        dna = SeqRecord(Seq("NNNNATGAAAGCCNNN"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf=(4, 13))

        assert str(result.seq) == "MKA"
        assert result.annotations["translation_region"] == (4, 13)
        assert result.annotations["translation_region_source"] == "explicit_orf"

    def test_attached_orf(self):
        """Test translation with attached ORF."""
        # ORF region [4:13] = ATGAAAGCC (9 bases, 3 codons) -> MKA
        dna = SeqRecord(Seq("NNNNATGAAAGCCNNN"), id="test")
        dna.annotations["molecule_type"] = "DNA"
        dna = attach_orf(dna, (4, 13))

        result = ensure_protein_record(dna)

        assert str(result.seq) == "MKA"
        assert result.annotations["translation_region"] == (4, 13)
        assert result.annotations["translation_region_source"] == "attached_orf"

    def test_explicit_orf_overrides_attached(self):
        """Test that explicit ORF overrides attached ORF."""
        # Attached ORF [4:13] = ATGAAAGCC -> MKA (9 bases, 3 codons)
        # Explicit ORF [16:22] = CTGGTG -> LV (6 bases, 2 codons)
        dna = SeqRecord(Seq("NNNNATGAAAGCCNNNCTGGTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"
        dna = attach_orf(dna, (4, 13))  # MKA

        result = ensure_protein_record(dna, orf=(16, 22))

        assert str(result.seq) == "LV"
        assert result.annotations["translation_region"] == (16, 22)
        assert result.annotations["translation_region_source"] == "explicit_orf"

    def test_use_orf_if_present_false(self):
        """Test that use_orf_if_present=False ignores attached ORF."""
        # Use a sequence with length divisible by 3 to avoid partial codon trimming
        # 18 bases = 6 codons
        dna = SeqRecord(Seq("NNNNATGAAAGCCNNNAA"), id="test")
        dna.annotations["molecule_type"] = "DNA"
        dna = attach_orf(dna, (4, 13))

        result = ensure_protein_record(dna, use_orf_if_present=False)

        # Should translate full sequence (with N's causing issues or being ambiguous)
        # Full sequence frame 0: NNNNATGAAAGCCNNNAA -> XMKAXN (N translates to X typically)
        assert result.annotations["translation_region"] == (0, len(dna.seq))
        assert result.annotations["translation_region_source"] == "full_sequence_frame0"

    def test_full_sequence_without_orf(self):
        """Test translation of full sequence without any ORF."""
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna)

        assert str(result.seq) == "MKALV"
        assert result.annotations["translation_region"] == (0, 15)
        assert result.annotations["translation_region_source"] == "full_sequence_frame0"


class TestEnsureProteinRecordEdgeCases:
    """Tests for edge cases and error handling."""

    def test_invalid_on_internal_stop(self):
        """Test that invalid on_internal_stop raises error."""
        dna = SeqRecord(Seq("ATGAAAGCC"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        with pytest.raises(ValueError, match="on_internal_stop must be"):
            ensure_protein_record(dna, on_internal_stop="invalid")

    def test_missing_molecule_type(self):
        """Test that missing molecule_type raises error."""
        record = SeqRecord(Seq("ATGAAAGCC"), id="test")

        with pytest.raises(KeyError, match="molecule_type"):
            ensure_protein_record(record)

    def test_invalid_molecule_type(self):
        """Test that invalid molecule_type raises error."""
        record = SeqRecord(Seq("ATGAAAGCC"), id="test")
        record.annotations["molecule_type"] = "invalid"

        with pytest.raises(ValueError, match="molecule_type must be"):
            ensure_protein_record(record)

    def test_empty_sequence(self):
        """Test translation of empty sequence."""
        dna = SeqRecord(Seq(""), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna)
        assert str(result.seq) == ""

    def test_short_sequence(self):
        """Test translation of sequence shorter than codon."""
        dna = SeqRecord(Seq("AT"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna)
        # 2 bases: incomplete codon, trimmed to 0 bases, results in empty protein
        assert len(result.seq) == 0

    def test_partial_codon_trimming(self):
        """Test that partial codons are trimmed during translation.

        Sequences with length not divisible by 3 have trailing bases trimmed
        to avoid BiopythonWarning. This is consistent with standard translation
        behavior where incomplete codons at the end are ignored.
        """
        # 7 bases: ATG GCT + 1 extra base -> should translate 6 bases -> MK
        dna7 = SeqRecord(Seq("ATGAAAA"), id="test7")
        dna7.annotations["molecule_type"] = "DNA"
        result7 = ensure_protein_record(dna7)
        assert str(result7.seq) == "MK"  # ATG + AAA = MK, trailing A trimmed

        # 8 bases: ATG GCT + 2 extra bases -> should translate 6 bases -> MK
        dna8 = SeqRecord(Seq("ATGAAAAA"), id="test8")
        dna8.annotations["molecule_type"] = "DNA"
        result8 = ensure_protein_record(dna8)
        assert str(result8.seq) == "MK"  # ATG + AAA = MK, trailing AA trimmed

        # 9 bases: ATG AAA GCC -> should translate all 9 bases -> MKA
        dna9 = SeqRecord(Seq("ATGAAAGCC"), id="test9")
        dna9.annotations["molecule_type"] = "DNA"
        result9 = ensure_protein_record(dna9)
        assert str(result9.seq) == "MKA"  # All 9 bases translated

    def test_preserves_id_and_description(self):
        """Test that ID and description are preserved."""
        dna = SeqRecord(Seq("ATGAAAGCC"), id="seq123", description="test sequence")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna)

        assert result.id == "seq123"
        assert result.description == "test sequence"
