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

    This is a generic windowing function that works with DNA, RNA, and protein sequences.
    It does NOT enforce any restrictions on step or window_size being multiples of 3.
    For codon-specific or ORF windowing with such restrictions, use iter_orf_codon_windows
    (see Ticket 5).

    Args:
        record: The SeqRecord to extract windows from (DNA, RNA, or protein)
        window_size: Size of each window in residues/bases (must be positive, no restrictions)
        step: Step size between windows in residues/bases (must be positive, no restrictions)
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
        >>> # DNA sequence with step=2
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

        >>> # Protein sequence with step=1 (sliding window)
        >>> protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein1")
        >>> protein_windows = list(iter_windows(protein_record, window_size=5, step=1))
        >>> len(protein_windows)
        5
        >>> str(protein_windows[0].seq)
        'MKALV'
        >>> str(protein_windows[1].seq)
        'KALVS'

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


def iter_orf_codon_windows(
    record: SeqRecord,
    *,
    orf: tuple[int, int],
    window_nt: int,
    step_nt: int,
    drop_partial: bool = True,
) -> Iterator[SeqRecord]:
    """Iterate over codon-aligned sliding windows within an ORF.

    This function extracts sliding windows from a specified ORF region of a DNA/RNA sequence,
    with enforcement that step_nt must be a multiple of 3 (codon-aligned). This is the
    function referenced in Ticket 5.

    Args:
        record: DNA or RNA SeqRecord
        orf: Tuple of (start, end) coordinates for the ORF (0-based, end-exclusive)
        window_nt: Size of each window in nucleotides (must be positive)
        step_nt: Step size between windows in nucleotides (must be positive and multiple of 3)
        drop_partial: If True, drops windows smaller than window_nt at the end

    Yields:
        SeqRecord objects for each window with annotations:
        - parent_id: ID of the parent sequence
        - orf_start: Start position of ORF in parent sequence
        - orf_end: End position of ORF in parent sequence
        - window_start: Start position of window relative to ORF start (0-indexed)
        - window_end: End position of window relative to ORF start (exclusive)
        - window_index: Index of this window (0-indexed)
        - start: Absolute start position in parent sequence
        - end: Absolute end position in parent sequence
        - All annotations from the parent record are also copied

    Raises:
        ValueError: If step_nt is not a multiple of 3, or if window_nt or step_nt are not positive

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> # Create record
        >>> record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
        >>> # Extract windows with step=3 (codon-aligned)
        >>> windows = list(iter_orf_codon_windows(
        ...     record, orf=(0, 15), window_nt=9, step_nt=3
        ... ))
        >>> len(windows)
        3
        >>> str(windows[0].seq)
        'ATGAAACCC'
        >>> windows[0].annotations["window_start"]
        0
        >>> windows[0].annotations["window_index"]
        0
        >>> str(windows[1].seq)
        'AAACCCGGG'
        >>> windows[1].annotations["window_start"]
        3
        >>> windows[1].annotations["window_index"]
        1
    """
    # Validate step is multiple of 3
    if step_nt % 3 != 0:
        raise ValueError(
            f"step_nt must be a multiple of 3 for codon-aligned windows, got {step_nt}"
        )

    if window_nt <= 0:
        raise ValueError(f"window_nt must be positive, got {window_nt}")

    if step_nt <= 0:
        raise ValueError(f"step_nt must be positive, got {step_nt}")

    # Extract ORF coordinates
    orf_start, orf_end = orf

    # Extract ORF sequence
    seq_str = get_seq_str(record)
    orf_seq_str = seq_str[orf_start:orf_end]
    orf_len = len(orf_seq_str)

    # Generate windows within the ORF
    window_start_in_orf = 0
    window_index = 0

    while window_start_in_orf < orf_len:
        window_end_in_orf = min(window_start_in_orf + window_nt, orf_len)
        window_seq_str = orf_seq_str[window_start_in_orf:window_end_in_orf]

        # Skip partial windows if drop_partial is True
        if drop_partial and len(window_seq_str) < window_nt:
            break

        # Calculate absolute positions in parent sequence
        abs_start = orf_start + window_start_in_orf
        abs_end = orf_start + window_end_in_orf

        # Create new SeqRecord for window
        window_record = SeqRecord(
            Seq(window_seq_str),
            id=f"{record.id}_orf_window_{window_index}",
            description=f"ORF window {window_index} at {window_start_in_orf} from {record.id}",
        )

        # Copy parent annotations
        window_record.annotations.update(record.annotations)

        # Add window-specific metadata
        window_record.annotations["parent_id"] = record.id  # type: ignore[assignment]
        window_record.annotations["orf_start"] = orf_start
        window_record.annotations["orf_end"] = orf_end
        window_record.annotations["window_start"] = window_start_in_orf
        window_record.annotations["window_end"] = window_end_in_orf
        window_record.annotations["window_index"] = window_index
        window_record.annotations["start"] = abs_start
        window_record.annotations["end"] = abs_end

        yield window_record

        window_index += 1
        window_start_in_orf += step_nt
