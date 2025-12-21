"""Tests for core record utilities."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.errors import InvalidSequenceError
from biotooler.core.record import coerce_record, get_molecule_type


class TestCoerceRecord:
    """Tests for coerce_record function."""

    def test_normalize_whitespace(self):
        """Test that whitespace is stripped from sequences."""
        # Test with various whitespace
        record = coerce_record("AC GT\n\tAA", "DNA")
        assert str(record.seq) == "ACGTAA"

        record = coerce_record(" A C G T ", "DNA")
        assert str(record.seq) == "ACGT"

    def test_normalize_lowercase(self):
        """Test that sequences are converted to uppercase."""
        record = coerce_record("acgt", "DNA")
        assert str(record.seq) == "ACGT"

        record = coerce_record("AcGt", "DNA")
        assert str(record.seq) == "ACGT"

    def test_accepts_string(self):
        """Test that string sequences are accepted."""
        record = coerce_record("ACGT", "DNA")
        assert str(record.seq) == "ACGT"
        assert record.annotations["molecule_type"] == "DNA"

    def test_accepts_seq(self):
        """Test that Bio.Seq.Seq objects are accepted."""
        seq = Seq("ACGT")
        record = coerce_record(seq, "DNA")
        assert str(record.seq) == "ACGT"
        assert record.annotations["molecule_type"] == "DNA"

    def test_accepts_seqrecord(self):
        """Test that SeqRecord objects are accepted."""
        seq_record = SeqRecord(Seq("ACGT"), id="test_id", description="test description")
        record = coerce_record(seq_record, "DNA")
        assert str(record.seq) == "ACGT"
        assert record.annotations["molecule_type"] == "DNA"

    def test_seqrecord_preserves_id_when_not_overridden(self):
        """Test that existing SeqRecord id is preserved when not explicitly overridden."""
        seq_record = SeqRecord(Seq("ACGT"), id="original_id", description="original description")
        record = coerce_record(seq_record, "DNA")
        assert record.id == "original_id"
        assert record.description == "original description"

    def test_seqrecord_overrides_id_when_provided(self):
        """Test that explicit id/description override existing SeqRecord values."""
        seq_record = SeqRecord(Seq("ACGT"), id="original_id", description="original description")
        record = coerce_record(seq_record, "DNA", id="new_id", description="new description")
        assert record.id == "new_id"
        assert record.description == "new description"

    def test_sets_default_id_for_string(self):
        """Test that default id is set for string sequences."""
        record = coerce_record("ACGT", "DNA")
        assert record.id == "<unknown id>"
        assert record.description == "<unknown description>"

    def test_molecule_type_annotation_dna(self):
        """Test that DNA molecule_type is set correctly."""
        record = coerce_record("ACGT", "DNA")
        assert record.annotations["molecule_type"] == "DNA"

    def test_molecule_type_annotation_rna(self):
        """Test that RNA molecule_type is set correctly."""
        record = coerce_record("ACGU", "RNA")
        assert record.annotations["molecule_type"] == "RNA"

    def test_molecule_type_annotation_protein(self):
        """Test that protein molecule_type is set correctly."""
        record = coerce_record("ACDEFGHIKLMNPQRSTVWY", "protein")
        assert record.annotations["molecule_type"] == "protein"

    def test_invalid_molecule_type(self):
        """Test that invalid molecule_type raises ValueError."""
        with pytest.raises(ValueError, match="molecule_type must be"):
            coerce_record("ACGT", "invalid")

    def test_dna_valid_characters(self):
        """Test that standard DNA characters are accepted."""
        record = coerce_record("ACGT", "DNA")
        assert str(record.seq) == "ACGT"

    def test_rna_valid_characters(self):
        """Test that standard RNA characters are accepted."""
        record = coerce_record("ACGU", "RNA")
        assert str(record.seq) == "ACGU"

    def test_protein_valid_characters(self):
        """Test that standard protein characters are accepted."""
        record = coerce_record("ACDEFGHIKLMNPQRSTVWY", "protein")
        assert str(record.seq) == "ACDEFGHIKLMNPQRSTVWY"

    def test_dna_ambiguous_allowed(self):
        """Test that ambiguous DNA characters are accepted when allow_ambiguous=True."""
        record = coerce_record("ACGTN", "DNA", allow_ambiguous=True)
        assert str(record.seq) == "ACGTN"

        # Test more ambiguous codes
        record = coerce_record("NRYWSMKHBVD", "DNA", allow_ambiguous=True)
        assert str(record.seq) == "NRYWSMKHBVD"

    def test_rna_ambiguous_allowed(self):
        """Test that ambiguous RNA characters are accepted when allow_ambiguous=True."""
        record = coerce_record("ACGUN", "RNA", allow_ambiguous=True)
        assert str(record.seq) == "ACGUN"

    def test_protein_ambiguous_allowed(self):
        """Test that ambiguous protein characters are accepted when allow_ambiguous=True."""
        record = coerce_record("ACDEFGHIKLMNPQRSTVWYX", "protein", allow_ambiguous=True)
        assert str(record.seq) == "ACDEFGHIKLMNPQRSTVWYX"

    def test_dna_ambiguous_rejected(self):
        """Test that ambiguous DNA characters are rejected when allow_ambiguous=False."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACGTN", "DNA", allow_ambiguous=False)

        assert exc_info.value.char == "N"
        assert exc_info.value.index == 4
        assert exc_info.value.molecule_type == "DNA"
        assert "position 4" in str(exc_info.value)

    def test_rna_ambiguous_rejected(self):
        """Test that ambiguous RNA characters are rejected when allow_ambiguous=False."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACGUN", "RNA", allow_ambiguous=False)

        assert exc_info.value.char == "N"
        assert exc_info.value.index == 4
        assert exc_info.value.molecule_type == "RNA"

    def test_protein_ambiguous_rejected(self):
        """Test that ambiguous protein characters are rejected when allow_ambiguous=False."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACDEFGHIKLMNPQRSTVWYX", "protein", allow_ambiguous=False)

        assert exc_info.value.char == "X"
        assert exc_info.value.index == 20
        assert exc_info.value.molecule_type == "protein"

    def test_invalid_dna_character(self):
        """Test that invalid DNA characters raise InvalidSequenceError."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACGTU", "DNA")

        assert exc_info.value.char == "U"
        assert exc_info.value.index == 4
        assert exc_info.value.molecule_type == "DNA"
        assert "Invalid character 'U' at position 4 for DNA" in str(exc_info.value)

    def test_invalid_rna_character(self):
        """Test that invalid RNA characters raise InvalidSequenceError."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACGUT", "RNA")

        assert exc_info.value.char == "T"
        assert exc_info.value.index == 4
        assert exc_info.value.molecule_type == "RNA"

    def test_invalid_protein_character(self):
        """Test that invalid protein characters raise InvalidSequenceError."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACDEFGHIKLMNPQRSTVWY1", "protein")

        assert exc_info.value.char == "1"
        assert exc_info.value.index == 20
        assert exc_info.value.molecule_type == "protein"

    def test_invalid_character_correct_index(self):
        """Test that the reported index for invalid characters is correct."""
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACGTX", "DNA")

        assert exc_info.value.index == 4
        assert exc_info.value.char == "X"

    def test_invalid_character_after_normalization(self):
        """Test that invalid characters are detected after normalization."""
        # The 'x' will be uppercased to 'X' which is invalid for DNA
        with pytest.raises(InvalidSequenceError) as exc_info:
            coerce_record("ACGTx", "DNA")

        assert exc_info.value.char == "X"

    def test_complex_example(self):
        """Test a complex example with all features."""
        record = coerce_record(
            " ac gt n ",
            "DNA",
            id="seq1",
            description="Test sequence with whitespace and ambiguous character",
            allow_ambiguous=True,
        )
        assert str(record.seq) == "ACGTN"
        assert record.id == "seq1"
        assert record.description == "Test sequence with whitespace and ambiguous character"
        assert record.annotations["molecule_type"] == "DNA"


class TestGetMoleculeType:
    """Tests for get_molecule_type function."""

    def test_returns_dna(self):
        """Test that DNA molecule type is returned correctly."""
        record = coerce_record("ACGT", "DNA")
        assert get_molecule_type(record) == "DNA"

    def test_returns_rna(self):
        """Test that RNA molecule type is returned correctly."""
        record = coerce_record("ACGU", "RNA")
        assert get_molecule_type(record) == "RNA"

    def test_returns_protein(self):
        """Test that protein molecule type is returned correctly."""
        record = coerce_record("ACDEFG", "protein")
        assert get_molecule_type(record) == "protein"

    def test_raises_error_when_missing(self):
        """Test that KeyError is raised when molecule_type annotation is missing."""
        record = SeqRecord(Seq("ACGT"), id="test")
        with pytest.raises(KeyError) as exc_info:
            get_molecule_type(record)

        assert "missing 'molecule_type' annotation" in str(exc_info.value)
        assert "test" in str(exc_info.value)
        assert "coerce_record()" in str(exc_info.value)


class TestInvalidSequenceError:
    """Tests for InvalidSequenceError class."""

    def test_error_attributes(self):
        """Test that InvalidSequenceError stores attributes correctly."""
        error = InvalidSequenceError("Test error", char="X", index=5, molecule_type="DNA")
        assert str(error) == "Test error"
        assert error.char == "X"
        assert error.index == 5
        assert error.molecule_type == "DNA"

    def test_error_is_value_error(self):
        """Test that InvalidSequenceError is a subclass of ValueError."""
        error = InvalidSequenceError("Test")
        assert isinstance(error, ValueError)
