"""Sequence utilities for efficient sequence string operations with caching."""

from typing import cast

from Bio.SeqRecord import SeqRecord


def get_seq_str(record: SeqRecord) -> str:
    """Get normalized sequence string with caching.

    Returns the sequence as an uppercase string with all whitespace and newlines removed.
    The result is cached in record.annotations["_biotooler_seq_str"] for efficient reuse.

    Args:
        record: A SeqRecord object

    Returns:
        Normalized sequence as uppercase string

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("acgt"), id="seq1")
        >>> get_seq_str(record)
        'ACGT'
        >>> # Second call returns cached value
        >>> get_seq_str(record)
        'ACGT'
    """
    # Check if cached value exists
    cache_key = "_biotooler_seq_str"
    if cache_key in record.annotations:
        return cast(str, record.annotations[cache_key])

    # Get sequence string from record
    seq_str = str(record.seq)

    # Normalize: remove whitespace and convert to uppercase
    seq_str = seq_str.replace(" ", "").replace("\t", "").replace("\n", "").replace("\r", "")
    seq_str = seq_str.upper()

    # Cache the result
    record.annotations[cache_key] = seq_str  # type: ignore[assignment]

    return seq_str


def get_seq_bytes(record: SeqRecord) -> bytes:
    """Get normalized sequence as bytes with caching.

    Returns the sequence as uppercase bytes with all whitespace and newlines removed.
    The result is cached in record.annotations["_biotooler_seq_bytes"] for efficient reuse.

    Args:
        record: A SeqRecord object

    Returns:
        Normalized sequence as uppercase bytes

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("acgt"), id="seq1")
        >>> get_seq_bytes(record)
        b'ACGT'
        >>> # Second call returns cached value
        >>> get_seq_bytes(record)
        b'ACGT'
    """
    # Check if cached value exists
    cache_key = "_biotooler_seq_bytes"
    if cache_key in record.annotations:
        return cast(bytes, record.annotations[cache_key])

    # Get sequence string and convert to bytes
    seq_str = get_seq_str(record)
    seq_bytes = seq_str.encode("ascii")

    # Cache the result
    record.annotations[cache_key] = seq_bytes  # type: ignore[assignment]

    return seq_bytes
