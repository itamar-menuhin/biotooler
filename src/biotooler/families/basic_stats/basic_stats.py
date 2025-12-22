"""Basic statistics feature computation for DNA/RNA/protein sequences."""

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_str
from biotooler.core.types import Scalar

# Standard alphabets
DNA_ALPHABET = "ACGT"
RNA_ALPHABET = "ACGU"
PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


class BasicStatsFeature:
    """Feature that computes basic sequence statistics.

    This feature computes nucleotide/amino acid counts, fractions, and length.
    For DNA/RNA sequences, it also computes GC content.

    The feature supports three sequence types:
    - DNA: Counts A, C, G, T (and N if present)
    - RNA: Counts A, C, G, U (and N if present)
    - Protein: Counts 20 standard amino acids (and X, * if present)

    Supports both baseline and incremental computation modes:
    - Baseline: Counts entire sequence for each window
    - Incremental: Maintains running counts and updates by subtracting outgoing
      and adding incoming subsequences

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = BasicStatsFeature()
        >>> # DNA sequence
        >>> record = SeqRecord(Seq("ACGTACGT"), id="test")
        >>> result = feature(record)
        >>> result["length"]
        8
        >>> result["count_a"]
        2
        >>> result["gc_fraction"]
        0.5
    """

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute basic statistics for the entire sequence.

        Args:
            record: DNA, RNA, or protein SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to scalar values:
            - length: Total sequence length
            - count_<char>: Count of each character (lowercase)
            - fraction_<char>: Fraction of each character (lowercase)
            - gc_fraction: GC content (DNA/RNA only)
        """
        # Get sequence string - validate cache to avoid window annotation issue
        # Windows created by iter_windows may copy parent's cached sequence
        cache_key = "_biotooler_seq_str"
        if cache_key in record.annotations:
            cached_seq = record.annotations[cache_key]
            actual_seq = str(record.seq).upper()
            if isinstance(cached_seq, str) and cached_seq != actual_seq:
                # Cache is stale (from parent record) - use actual sequence
                seq_str = actual_seq
            elif isinstance(cached_seq, str):
                seq_str = cached_seq
            else:
                # Cache value is not a string - use actual sequence
                seq_str = actual_seq
        else:
            seq_str = get_seq_str(record)
        return self._compute_stats(seq_str)

    def _compute_stats(self, seq_str: str) -> dict[str, Scalar]:
        """Compute statistics from a sequence string.

        Args:
            seq_str: Uppercase sequence string

        Returns:
            Dictionary of statistics
        """
        length = len(seq_str)

        # Detect sequence type
        alphabet = self._detect_alphabet(seq_str)

        # Count each character in alphabet
        counts = {}
        for char in alphabet:
            counts[char] = seq_str.count(char)

        # Build result dictionary
        result: dict[str, Scalar] = {"length": length}

        # Add counts
        for char, count in counts.items():
            result[f"count_{char.lower()}"] = count

        # Add fractions
        if length > 0:
            for char, count in counts.items():
                result[f"fraction_{char.lower()}"] = count / length
        else:
            for char in alphabet:
                result[f"fraction_{char.lower()}"] = 0.0

        # Add GC fraction for DNA/RNA
        base_alphabet = alphabet.rstrip("N")  # Remove N if present
        if base_alphabet in (DNA_ALPHABET, RNA_ALPHABET):
            if length > 0:
                gc_count = counts.get("G", 0) + counts.get("C", 0)
                result["gc_fraction"] = gc_count / length
            else:
                result["gc_fraction"] = 0.0

        return result

    def _detect_alphabet(self, seq_str: str) -> str:
        """Detect sequence alphabet from content.

        Args:
            seq_str: Uppercase sequence string

        Returns:
            Alphabet string (DNA_ALPHABET, RNA_ALPHABET, or PROTEIN_ALPHABET)
        """
        # Check for U -> RNA
        if "U" in seq_str:
            # RNA may also have N
            if "N" in seq_str:
                return RNA_ALPHABET + "N"
            return RNA_ALPHABET

        # Check for T -> DNA
        if "T" in seq_str:
            # DNA may also have N
            if "N" in seq_str:
                return DNA_ALPHABET + "N"
            return DNA_ALPHABET

        # Check for protein-specific amino acids
        protein_specific = set(seq_str) - set(DNA_ALPHABET) - {"N"}
        if protein_specific:
            # Protein sequence
            alphabet_set = set(PROTEIN_ALPHABET)
            # Add X and * if present
            if "X" in seq_str:
                alphabet_set.add("X")
            if "*" in seq_str:
                alphabet_set.add("*")
            return "".join(sorted(alphabet_set))

        # Ambiguous - could be DNA without T or protein with only ACGN
        # Default to DNA if only ACGN present
        if set(seq_str) <= set(DNA_ALPHABET + "N"):
            if "N" in seq_str:
                return DNA_ALPHABET + "N"
            return DNA_ALPHABET

        # Default to protein
        alphabet_set = set(PROTEIN_ALPHABET)
        if "X" in seq_str:
            alphabet_set.add("X")
        if "*" in seq_str:
            alphabet_set.add("*")
        return "".join(sorted(alphabet_set))

    def init_state(
        self,
        record: SeqRecord,
        *,
        window_start: int,
        window_end: int,
        **kwargs,
    ) -> dict:
        """Initialize state for the first window.

        Args:
            record: The SeqRecord containing the sequence
            window_start: Start position of window
            window_end: End position of window
            **kwargs: Additional parameters (region, orf - ignored for basic_stats)

        Returns:
            State dictionary containing:
            - record: The SeqRecord
            - counts: numpy array of character counts
            - alphabet: String of characters being counted
            - char_to_idx: Mapping from character to index in counts array
        """
        # Extract window sequence
        seq_str = get_seq_str(record)
        window_seq = seq_str[window_start:window_end]

        # Detect alphabet from the full sequence (not just window)
        # This ensures consistent alphabet across all windows
        alphabet = self._detect_alphabet(seq_str)

        # Create character to index mapping
        char_to_idx = {char: idx for idx, char in enumerate(alphabet)}

        # Initialize counts array
        counts = np.zeros(len(alphabet), dtype=np.int64)

        # Count characters in the window
        for char in window_seq:
            if char in char_to_idx:
                counts[char_to_idx[char]] += 1

        state = {
            "record": record,
            "counts": counts,
            "alphabet": alphabet,
            "char_to_idx": char_to_idx,
        }

        return state

    def step_state(
        self,
        state: dict,
        *,
        out_start: int,
        out_end: int,
        in_start: int,
        in_end: int,
        **kwargs,
    ) -> None:
        """Update state for the next window by updating counts.

        Args:
            state: The state object to update
            out_start: Start position of bases leaving the window
            out_end: End position of bases leaving the window
            in_start: Start position of bases entering the window
            in_end: End position of bases entering the window
            **kwargs: Additional parameters
        """
        record = state["record"]
        counts = state["counts"]
        char_to_idx = state["char_to_idx"]
        seq_str = get_seq_str(record)

        # Remove outgoing characters
        if out_end > out_start:
            out_seq = seq_str[out_start:out_end]
            for char in out_seq:
                if char in char_to_idx:
                    counts[char_to_idx[char]] -= 1

        # Add incoming characters
        if in_end > in_start:
            in_seq = seq_str[in_start:in_end]
            for char in in_seq:
                if char in char_to_idx:
                    counts[char_to_idx[char]] += 1

    def emit(self, state: dict) -> dict[str, Scalar]:
        """Emit feature values from current state.

        Args:
            state: The state object

        Returns:
            Dictionary mapping feature names to scalar values
        """
        counts = state["counts"]
        alphabet = state["alphabet"]

        # Calculate length
        length = int(np.sum(counts))

        # Build result dictionary
        result: dict[str, Scalar] = {"length": length}

        # Add counts
        for idx, char in enumerate(alphabet):
            result[f"count_{char.lower()}"] = int(counts[idx])

        # Add fractions
        if length > 0:
            for idx, char in enumerate(alphabet):
                result[f"fraction_{char.lower()}"] = float(counts[idx]) / length
        else:
            for char in alphabet:
                result[f"fraction_{char.lower()}"] = 0.0

        # Add GC fraction for DNA/RNA
        base_alphabet = alphabet.rstrip("N")  # Remove N if present
        if base_alphabet in (DNA_ALPHABET, RNA_ALPHABET):
            if length > 0:
                g_idx = alphabet.index("G")
                c_idx = alphabet.index("C")
                gc_count = int(counts[g_idx]) + int(counts[c_idx])
                result["gc_fraction"] = float(gc_count) / length
            else:
                result["gc_fraction"] = 0.0

        return result
