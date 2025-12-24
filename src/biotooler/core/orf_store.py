"""Utilities for selecting and attaching ORF spans to sequence records.

This module provides functions for working with ORF (Open Reading Frame) spans
in the context of batch workflows and single-sequence selection.
"""

from typing import Any

from Bio.SeqRecord import SeqRecord

from biotooler.core.record import get_molecule_type

# Type alias for ORF span: (start, end) where coordinates are 0-based, end-exclusive
OrfSpan = tuple[int, int]


def select_orf_by_index(candidates: list[OrfSpan], orf_index: int) -> OrfSpan:
    """Select a single ORF from a list of candidates by index.

    Args:
        candidates: List of ORF spans as (start, end) tuples
        orf_index: Zero-based index of the ORF to select

    Returns:
        The selected ORF span as a (start, end) tuple

    Raises:
        IndexError: If orf_index is out of range

    Examples:
        >>> candidates = [(0, 9), (9, 18), (18, 27)]
        >>> select_orf_by_index(candidates, 0)
        (0, 9)
        >>> select_orf_by_index(candidates, 1)
        (9, 18)
        >>> select_orf_by_index(candidates, 3)
        Traceback (most recent call last):
            ...
        IndexError: ORF index 3 is out of range for 3 candidates (valid indices: 0-2)
    """
    if not 0 <= orf_index < len(candidates):
        if len(candidates) == 0:
            raise IndexError(f"ORF index {orf_index} is out of range: no candidates available")
        else:
            raise IndexError(
                f"ORF index {orf_index} is out of range for {len(candidates)} candidates "
                f"(valid indices: 0-{len(candidates) - 1})"
            )
    return candidates[orf_index]


def attach_orf(record: SeqRecord, orf: OrfSpan, *, key: str = "biotooler.orf") -> SeqRecord:
    """Attach an ORF span to a SeqRecord's annotations.

    This function stores the ORF span in the record's annotations dictionary
    for later retrieval. The record is modified in place and also returned
    for convenience.

    Args:
        record: SeqRecord to attach the ORF to
        orf: ORF span as (start, end) tuple (0-based, end-exclusive)
        key: Annotation key to store the ORF under (default: "biotooler.orf")

    Returns:
        The same SeqRecord with the ORF attached (modified in place)

    Raises:
        ValueError: If the record is a protein sequence

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("ATGAAATAA"), id="test")
        >>> record = attach_orf(record, (0, 9))
        >>> record.annotations["biotooler.orf"]
        (0, 9)
    """
    # Check if this is a protein sequence
    if "molecule_type" in record.annotations:
        mol_type = get_molecule_type(record)
        if mol_type.upper() == "PROTEIN":
            raise ValueError(
                "ORF operations are only supported for DNA/RNA sequences, not protein sequences"
            )

    record.annotations[key] = orf  # type: ignore[assignment]
    return record


def get_orf(record: SeqRecord, *, key: str = "biotooler.orf") -> OrfSpan:
    """Retrieve an ORF span from a SeqRecord's annotations.

    Args:
        record: SeqRecord to retrieve the ORF from
        key: Annotation key where the ORF is stored (default: "biotooler.orf")

    Returns:
        The ORF span as a (start, end) tuple (0-based, end-exclusive)

    Raises:
        KeyError: If the ORF is not found in the record's annotations

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("ATGAAATAA"), id="test")
        >>> record = attach_orf(record, (0, 9))
        >>> get_orf(record)
        (0, 9)
        >>> get_orf(record, key="nonexistent")
        Traceback (most recent call last):
            ...
        KeyError: "ORF not found in record 'test': no annotation with key 'nonexistent'..."
    """
    try:
        orf: Any = record.annotations[key]
        return orf  # type: ignore[return-value]
    except KeyError as e:
        raise KeyError(
            f"ORF not found in record {record.id!r}: no annotation with key {key!r}. "
            "Use attach_orf() to store an ORF first."
        ) from e
