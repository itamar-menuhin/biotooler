"""Tests for sequence utilities."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_bytes, get_seq_str


class TestGetSeqStr:
    """Tests for get_seq_str function."""

    def test_returns_uppercase(self):
        """Test that sequence is converted to uppercase."""
        record = SeqRecord(Seq("acgt"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGT"

    def test_removes_whitespace(self):
        """Test that whitespace is removed from sequence."""
        record = SeqRecord(Seq("AC GT"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGT"

    def test_removes_newlines(self):
        """Test that newlines are removed from sequence."""
        record = SeqRecord(Seq("AC\nGT"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGT"

    def test_removes_tabs(self):
        """Test that tabs are removed from sequence."""
        record = SeqRecord(Seq("AC\tGT"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGT"

    def test_removes_carriage_returns(self):
        """Test that carriage returns are removed from sequence."""
        record = SeqRecord(Seq("AC\rGT"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGT"

    def test_normalizes_mixed_case_and_whitespace(self):
        """Test normalization with mixed case and whitespace."""
        record = SeqRecord(Seq("ac gt\n\tAA"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGTAA"

    def test_caches_result(self):
        """Test that result is cached in annotations."""
        record = SeqRecord(Seq("acgt"), id="seq1")
        result1 = get_seq_str(record)

        # Check that result is cached
        assert "_biotooler_seq_str" in record.annotations
        assert record.annotations["_biotooler_seq_str"] == "ACGT"

        # Second call should return cached value
        result2 = get_seq_str(record)
        assert result1 == result2
        assert result2 == "ACGT"

    def test_returns_cached_value(self):
        """Test that cached value is returned on subsequent calls."""
        record = SeqRecord(Seq("acgt"), id="seq1")

        # First call
        result1 = get_seq_str(record)
        assert result1 == "ACGT"

        # Modify the cached value to verify it's being used
        record.annotations["_biotooler_seq_str"] = "MODIFIED"

        # Second call should return modified cached value
        result2 = get_seq_str(record)
        assert result2 == "MODIFIED"

    def test_works_with_empty_sequence(self):
        """Test that empty sequences are handled correctly."""
        record = SeqRecord(Seq(""), id="seq1")
        result = get_seq_str(record)
        assert result == ""

    def test_preserves_valid_characters(self):
        """Test that valid sequence characters are preserved."""
        record = SeqRecord(Seq("ACGTRYSWKMBDHVN"), id="seq1")
        result = get_seq_str(record)
        assert result == "ACGTRYSWKMBDHVN"

    def test_different_records_have_separate_cache(self):
        """Test that different records have independent caches."""
        record1 = SeqRecord(Seq("acgt"), id="seq1")
        record2 = SeqRecord(Seq("cgta"), id="seq2")

        result1 = get_seq_str(record1)
        result2 = get_seq_str(record2)

        assert result1 == "ACGT"
        assert result2 == "CGTA"
        assert record1.annotations["_biotooler_seq_str"] == "ACGT"
        assert record2.annotations["_biotooler_seq_str"] == "CGTA"


class TestGetSeqBytes:
    """Tests for get_seq_bytes function."""

    def test_returns_bytes(self):
        """Test that result is bytes type."""
        record = SeqRecord(Seq("acgt"), id="seq1")
        result = get_seq_bytes(record)
        assert isinstance(result, bytes)
        assert result == b"ACGT"

    def test_normalizes_to_uppercase(self):
        """Test that sequence is converted to uppercase bytes."""
        record = SeqRecord(Seq("acgt"), id="seq1")
        result = get_seq_bytes(record)
        assert result == b"ACGT"

    def test_removes_whitespace(self):
        """Test that whitespace is removed from sequence."""
        record = SeqRecord(Seq("AC GT"), id="seq1")
        result = get_seq_bytes(record)
        assert result == b"ACGT"

    def test_caches_result(self):
        """Test that result is cached in annotations."""
        record = SeqRecord(Seq("acgt"), id="seq1")
        result1 = get_seq_bytes(record)

        # Check that result is cached
        assert "_biotooler_seq_bytes" in record.annotations
        assert record.annotations["_biotooler_seq_bytes"] == b"ACGT"

        # Second call should return cached value
        result2 = get_seq_bytes(record)
        assert result1 == result2
        assert result2 == b"ACGT"

    def test_returns_cached_value(self):
        """Test that cached value is returned on subsequent calls."""
        record = SeqRecord(Seq("acgt"), id="seq1")

        # First call
        result1 = get_seq_bytes(record)
        assert result1 == b"ACGT"

        # Modify the cached value to verify it's being used
        record.annotations["_biotooler_seq_bytes"] = b"MODIFIED"

        # Second call should return modified cached value
        result2 = get_seq_bytes(record)
        assert result2 == b"MODIFIED"

    def test_uses_get_seq_str_internally(self):
        """Test that get_seq_bytes uses get_seq_str for normalization."""
        record = SeqRecord(Seq("acgt"), id="seq1")

        # Call get_seq_bytes
        result = get_seq_bytes(record)
        assert result == b"ACGT"

        # Verify that get_seq_str cache was also created
        assert "_biotooler_seq_str" in record.annotations
        assert record.annotations["_biotooler_seq_str"] == "ACGT"

    def test_works_with_empty_sequence(self):
        """Test that empty sequences are handled correctly."""
        record = SeqRecord(Seq(""), id="seq1")
        result = get_seq_bytes(record)
        assert result == b""

    def test_different_records_have_separate_cache(self):
        """Test that different records have independent caches."""
        record1 = SeqRecord(Seq("acgt"), id="seq1")
        record2 = SeqRecord(Seq("cgta"), id="seq2")

        result1 = get_seq_bytes(record1)
        result2 = get_seq_bytes(record2)

        assert result1 == b"ACGT"
        assert result2 == b"CGTA"
        assert record1.annotations["_biotooler_seq_bytes"] == b"ACGT"
        assert record2.annotations["_biotooler_seq_bytes"] == b"CGTA"
