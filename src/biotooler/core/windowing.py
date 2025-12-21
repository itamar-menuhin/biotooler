"""Sliding window utilities for sequence analysis."""

from collections.abc import Iterator

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_str


def iter_windows(
    record: SeqRecord,
    window_size: int,
    step: int,
    *,
    drop_partial: bool = True,
) -> Iterator[SeqRecord]:
    """Iterate over sliding windows of a sequence.

    Efficiently generates sliding windows from a SeqRecord by extracting the sequence
    string once and slicing it for each window. Window metadata is stored in each
    window's annotations.

    Args:
        record: The SeqRecord to extract windows from
        window_size: Size of each window (must be positive)
        step: Step size between windows (must be positive)
        drop_partial: If True, drops the last window if it's smaller than window_size.
                     If False, includes partial windows at the end.

    Yields:
        SeqRecord objects for each window with the following annotations:
        - parent_id: ID of the parent sequence
        - start: Start position in parent sequence (0-indexed)
        - end: End position in parent sequence (exclusive)
        - window_index: Index of this window (0-indexed)
        - All annotations from the parent record are also copied

    Raises:
        ValueError: If window_size or step is not positive

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> record = SeqRecord(Seq("ACGTACGT"), id="seq1")
        >>> windows = list(iter_windows(record, window_size=4, step=2))
        >>> len(windows)
        3
        >>> str(windows[0].seq)
        'ACGT'
        >>> windows[0].annotations["start"]
        0
        >>> windows[0].annotations["end"]
        4

        >>> # With drop_partial=False
        >>> record = SeqRecord(Seq("ACGTACG"), id="seq1")
        >>> windows = list(iter_windows(record, window_size=4, step=2, drop_partial=False))
        >>> len(windows)
        4
        >>> str(windows[-1].seq)
        'CG'
    """
    # Validate arguments
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}")
    if step <= 0:
        raise ValueError(f"step must be positive, got {step}")

    # Get the sequence string once (uses caching internally)
    seq_str = get_seq_str(record)
    seq_len = len(seq_str)

    # Generate windows
    window_index = 0
    start = 0

    while start < seq_len:
        end = min(start + window_size, seq_len)
        window_seq_str = seq_str[start:end]

        # Skip partial windows if drop_partial is True
        if drop_partial and len(window_seq_str) < window_size:
            break

        # Create new SeqRecord for window
        window_record = SeqRecord(
            Seq(window_seq_str),
            id=f"{record.id}_window_{window_index}",
            description=f"Window {window_index} from {record.id} [{start}:{end}]",
        )

        # Copy parent annotations
        window_record.annotations.update(record.annotations)

        # Add window-specific metadata
        window_record.annotations["parent_id"] = record.id  # type: ignore[assignment]
        window_record.annotations["start"] = start
        window_record.annotations["end"] = end
        window_record.annotations["window_index"] = window_index

        yield window_record

        window_index += 1
        start += step
