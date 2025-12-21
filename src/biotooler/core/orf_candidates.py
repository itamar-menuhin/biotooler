"""Helper for enumerating ORF candidate spans based on start/stop codons.

This module provides functionality to identify all in-frame open reading frame (ORF)
candidates in a DNA or RNA sequence by finding start and stop codon pairs.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def find_orf_candidates(
    seq_or_record: str | Seq | SeqRecord,
    *,
    start: str = "ATG",
    stops: tuple[str, ...] = ("TAA", "TAG", "TGA"),
) -> list[tuple[int, int]]:
    """Find all in-frame ORF candidate spans based on start/stop codons.

    This function enumerates ALL in-frame start/stop codon pairs in the sequence.
    The output can be large for sequences with many start and stop codons.

    The function treats input as DNA. If the sequence contains 'U' characters (RNA),
    they are automatically converted to 'T' for scanning purposes.

    Args:
        seq_or_record: Input sequence as string, Seq, or SeqRecord
        start: Start codon sequence (default: "ATG")
        stops: Tuple of stop codon sequences (default: ("TAA", "TAG", "TGA"))

    Returns:
        List of (start_index, end_index) tuples where:
        - start_index is the position of the first base of the start codon
        - end_index is the position after the last base of the stop codon (stop_idx + 3)
        - Only pairs where stop_idx > start_idx and (stop_idx - start_idx) % 3 == 0
        - Sorted by start position, then by end position

    Examples:
        >>> from Bio.Seq import Seq
        >>> # Simple example with one ORF
        >>> find_orf_candidates("ATGAAATAA")
        [(0, 9)]

        >>> # Multiple starts and stops
        >>> find_orf_candidates("ATGAAATAGATGCCCTAGTAA")
        [(0, 9), (0, 18), (9, 18)]

        >>> # RNA example (U converted to T)
        >>> find_orf_candidates("AUGUGUUAA")
        [(0, 9)]

        >>> # No candidates
        >>> find_orf_candidates("AAACCCGGG")
        []
    """
    # Extract sequence string
    if isinstance(seq_or_record, SeqRecord):
        seq_str = str(seq_or_record.seq)
    elif isinstance(seq_or_record, Seq):
        seq_str = str(seq_or_record)
    else:
        seq_str = seq_or_record

    # Normalize: uppercase and convert U to T (RNA to DNA)
    seq_str = seq_str.upper().replace("U", "T")

    # Normalize start and stop codons
    start = start.upper().replace("U", "T")
    stops = tuple(s.upper().replace("U", "T") for s in stops)

    # Find all start codon positions
    start_positions = []
    for i in range(len(seq_str) - len(start) + 1):
        if seq_str[i : i + len(start)] == start:
            start_positions.append(i)

    # Find all stop codon positions
    stop_positions = []
    for stop in stops:
        for i in range(len(seq_str) - len(stop) + 1):
            if seq_str[i : i + len(stop)] == stop:
                stop_positions.append(i)

    # Find all valid in-frame pairs
    candidates = []
    for start_idx in start_positions:
        for stop_idx in stop_positions:
            # Check if stop is after start and in-frame
            if stop_idx > start_idx and (stop_idx - start_idx) % 3 == 0:
                # Return (start_idx, stop_idx + 3) as specified
                candidates.append((start_idx, stop_idx + 3))

    # Sort by start, then end
    candidates.sort()

    return candidates
