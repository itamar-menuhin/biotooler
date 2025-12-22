"""Sliding window utilities for sequence analysis."""

from collections.abc import Iterator

import numpy as np
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


def compute_window_indices(
    sequence_length: int,
    window_size: int,
    step: int,
    *,
    drop_partial: bool = True,
    position_space: str = "nucleotide",
    start_offset: int = 0,
) -> list[tuple[int, int]]:
    """Compute window start/end index pairs for sliding windows.

    This is a shared indexing logic helper that generates window boundaries
    with proper validation for different position spaces (nucleotide codon-space
    vs amino acid/residue space).

    For nucleotide codon-space windows (position_space="codon"):
    - step and start_offset must be multiples of 3 (enforced)
    - window_size must be positive
    - indices are in nucleotide space (not codon space)

    For amino acid/residue windows (position_space="residue"):
    - step=1 is allowed (and any positive step)
    - window_size must be positive
    - no codon alignment restrictions

    Args:
        sequence_length: Length of the sequence in nucleotides or residues
        window_size: Size of each window in nucleotides or residues (must be positive)
        step: Step size between windows in nucleotides or residues (must be positive)
        drop_partial: If True, drops windows smaller than window_size at the end
        position_space: Either "codon" (enforces multiples of 3) or "residue" (no restrictions)
        start_offset: Starting offset for the first window (default: 0)

    Returns:
        List of (start, end) tuples representing window boundaries in nucleotide/residue space

    Raises:
        ValueError: If parameters are invalid, or if position_space="codon" and
                   step or start_offset are not multiples of 3

    Examples:
        >>> # Codon-space windows (step must be multiple of 3)
        >>> compute_window_indices(15, window_size=9, step=3, position_space="codon")
        [(0, 9), (3, 12), (6, 15)]

        >>> # Residue-space windows (step=1 allowed)
        >>> compute_window_indices(10, window_size=5, step=1, position_space="residue")
        [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9), (5, 10)]

        >>> # Codon-space with start_offset
        >>> compute_window_indices(
        ...     15, window_size=9, step=3, position_space="codon", start_offset=3
        ... )
        [(3, 12), (6, 15)]

        >>> # Error: step not multiple of 3 for codon space
        >>> compute_window_indices(15, window_size=9, step=1, position_space="codon")
        Traceback (most recent call last):
        ...
        ValueError: For codon position space, step must be a multiple of 3, got 1
    """
    # Validate window_size
    if window_size <= 0:
        raise ValueError(f"window_size must be positive, got {window_size}")

    # Validate step
    if step <= 0:
        raise ValueError(f"step must be positive, got {step}")

    # Validate position_space
    if position_space not in ("codon", "residue"):
        raise ValueError(
            f"position_space must be 'codon' or 'residue', got {position_space!r}"
        )

    # Enforce codon alignment for codon space
    if position_space == "codon":
        if step % 3 != 0:
            raise ValueError(
                f"For codon position space, step must be a multiple of 3, got {step}"
            )
        if start_offset % 3 != 0:
            raise ValueError(
                f"For codon position space, start_offset must be a multiple of 3, "
                f"got {start_offset}"
            )

    # Generate window boundaries
    windows = []
    start = start_offset

    while start < sequence_length:
        end = min(start + window_size, sequence_length)

        # Skip partial windows if drop_partial is True
        if drop_partial and (end - start) < window_size:
            break

        windows.append((start, end))
        start += step

    return windows


def compute_window_index_arrays(
    sequence_length: int,
    window_size: int,
    step: int,
    *,
    drop_partial: bool = True,
    position_space: str = "nucleotide",
    start_offset: int = 0,
    output_space: str | None = None,
) -> list[np.ndarray]:
    """Compute index arrays for each window, optionally converting to codon/residue space.

    This helper produces numpy index arrays that can be used to slice position-space
    vectors (e.g., per-codon or per-residue feature vectors) for windowing operations.

    For nucleotide sequences with codon-space features:
    - position_space="codon" enforces step/start multiples of 3
    - output_space="codon" converts nucleotide indices to codon indices (divides by 3)
    - Returns arrays of codon indices for slicing codon-space vectors

    For amino acid sequences with residue-space features:
    - position_space="residue" allows any step (including step=1)
    - output_space="residue" (or None) keeps residue indices
    - Returns arrays of residue indices for slicing residue-space vectors

    Args:
        sequence_length: Length of the sequence in nucleotides or residues
        window_size: Size of each window in nucleotides or residues
        step: Step size between windows in nucleotides or residues
        drop_partial: If True, drops windows smaller than window_size
        position_space: Either "codon" or "residue" (for validation)
        start_offset: Starting offset for the first window (default: 0)
        output_space: Optional conversion space. If "codon", converts nucleotide
                     indices to codon indices by dividing by 3. If None or "residue",
                     keeps original indices.

    Returns:
        List of numpy arrays, one per window, containing indices in the output space

    Raises:
        ValueError: If parameters are invalid or incompatible

    Examples:
        >>> # Codon-space windows with codon output
        >>> arrays = compute_window_index_arrays(
        ...     15, window_size=9, step=3, position_space="codon", output_space="codon"
        ... )
        >>> [arr.tolist() for arr in arrays]
        [[0, 1, 2], [1, 2, 3], [2, 3, 4]]

        >>> # Residue-space windows with step=1
        >>> arrays = compute_window_index_arrays(
        ...     10, window_size=3, step=1, position_space="residue"
        ... )
        >>> [arr.tolist() for arr in arrays]
        [[0, 1, 2], [1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9]]
    """
    # Get window boundaries in nucleotide/residue space
    window_boundaries = compute_window_indices(
        sequence_length=sequence_length,
        window_size=window_size,
        step=step,
        drop_partial=drop_partial,
        position_space=position_space,
        start_offset=start_offset,
    )

    # Convert to index arrays
    index_arrays = []
    for start, end in window_boundaries:
        # Convert to output space if requested
        if output_space == "codon":
            # Convert nucleotide boundaries to codon indices
            # For codon space, start and end should be divisible by 3
            codon_start = start // 3
            codon_end = end // 3
            indices = np.arange(codon_start, codon_end)
        elif output_space is None or output_space == "residue":
            # Generate indices in the original space
            indices = np.arange(start, end)
        else:
            raise ValueError(
                f"output_space must be 'codon', 'residue', or None, got {output_space!r}"
            )

        index_arrays.append(indices)

    return index_arrays
