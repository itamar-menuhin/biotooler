"""Chimera index caching for efficient reuse of suffix arrays.

This module provides deterministic caching of chimera suffix arrays built from
ReferenceSequenceSet instances. Building suffix arrays is expensive, so this
caching layer enables reuse across multiple feature computations.
"""

import hashlib
from collections.abc import MutableMapping
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from biotooler.core.reference_sequences import ReferenceSequenceSet


@dataclass(frozen=True)
class ChimeraIndexKey:
    """Cache key for chimera reference index.

    Attributes:
        alphabet: Alphabet type for normalization ("nt", "aa", "codon")
        genetic_code_table: Genetic code table ID for translation (relevant for "aa" and "codon")
        ref_hash: Hash of normalized reference sequences
    """

    alphabet: str
    genetic_code_table: int
    ref_hash: str

    def __str__(self) -> str:
        """Return string representation of key."""
        return f"{self.alphabet}:gc{self.genetic_code_table}:{self.ref_hash[:16]}"


class ChimeraReferenceIndex:
    """Container for chimera suffix array and associated metadata.

    Attributes:
        suffix_array: The built suffix array from chimera.build_suffix_array
        key: The ChimeraIndexKey used to generate this index
    """

    def __init__(self, suffix_array: Any, key: ChimeraIndexKey):
        """Initialize ChimeraReferenceIndex.

        Args:
            suffix_array: The chimera suffix array
            key: The cache key for this index
        """
        self.suffix_array = suffix_array
        self.key = key


def get_reference_index(
    reference_set: "ReferenceSequenceSet",
    alphabet: str = "codon",
    genetic_code_table: int = 1,
    cache: MutableMapping[ChimeraIndexKey, ChimeraReferenceIndex] | None = None,
) -> ChimeraReferenceIndex:
    """Get or build a chimera reference index with caching.

    This function normalizes reference sequences according to the specified alphabet,
    computes a deterministic hash, and builds a suffix array using chimera's
    build_suffix_array function. Results are cached to avoid expensive rebuilds.

    Args:
        reference_set: ReferenceSequenceSet containing reference sequences
        alphabet: Alphabet type for normalization - "nt", "aa", or "codon" (default: "codon")
        genetic_code_table: Genetic code table ID for translation (default: 1)
        cache: Optional cache dict for storing built indices. If None, no caching occurs.

    Returns:
        ChimeraReferenceIndex containing the suffix array and cache key

    Raises:
        ValueError: If alphabet is not one of "nt", "aa", "codon"
        ValueError: If reference_set cannot provide required sequences
        ImportError: If pychimera is not installed

    Examples:
        >>> from biotooler.core.reference_sequences import ReferenceSequenceSet
        >>> ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATG"})
        >>> # Without cache (builds fresh each time)
        >>> index1 = get_reference_index(ref_set, alphabet="codon")
        >>> # With cache (reuses on second call)
        >>> cache = {}
        >>> index2 = get_reference_index(ref_set, alphabet="codon", cache=cache)
        >>> index3 = get_reference_index(ref_set, alphabet="codon", cache=cache)
        >>> # index3 reuses the suffix array from index2
    """
    # Lazy import chimera inside function
    from biotooler.families.chimera.integration import require_pychimera

    chimera = require_pychimera()

    # Validate alphabet
    if alphabet not in ("nt", "aa", "codon"):
        raise ValueError(f"alphabet must be 'nt', 'aa', or 'codon', got: {alphabet}")

    # Normalize reference sequences based on alphabet
    if alphabet == "nt":
        # Use CDS sequences directly (nucleotide)
        ref_sequences = reference_set.cds_strings(require_multiple_of_three=False)
    elif alphabet == "aa":
        # Use protein sequences (amino acid)
        ref_sequences = reference_set.protein_strings(
            strip_terminal_stop=True, error_on_internal_stop=False
        )
    else:  # codon
        # Get CDS strings for codon conversion
        cds_sequences = reference_set.cds_strings(require_multiple_of_three=True)
        # Convert to codon format using chimera's nt2codon
        ref_sequences = chimera.nt2codon(cds_sequences)

    # Compute deterministic hash from normalized sequences
    # Sort sequences for deterministic ordering, then concatenate
    sorted_sequences = sorted(ref_sequences)
    concatenated = "".join(sorted_sequences)
    ref_hash = hashlib.sha256(concatenated.encode()).hexdigest()

    # Create cache key
    # genetic_code_table is included for API consistency and future extensibility,
    # even though ref_hash already captures the sequence content
    key = ChimeraIndexKey(
        alphabet=alphabet, genetic_code_table=genetic_code_table, ref_hash=ref_hash
    )

    # Check cache if provided
    if cache is not None and key in cache:
        return cache[key]

    # Build suffix array using chimera
    if alphabet == "codon":
        # ref_sequences is already in codon format from nt2codon
        suffix_array = chimera.build_suffix_array(ref_sequences)
    elif alphabet == "aa":
        # Convert amino acids to codon format for chimera
        ref_sequences_codon = chimera.aa2codon(ref_sequences)
        suffix_array = chimera.build_suffix_array(ref_sequences_codon)
    else:  # nt
        # Convert nucleotides to codon format
        ref_sequences_codon = chimera.nt2codon(ref_sequences)
        suffix_array = chimera.build_suffix_array(ref_sequences_codon)

    # Create index
    index = ChimeraReferenceIndex(suffix_array=suffix_array, key=key)

    # Store in cache if provided
    if cache is not None:
        cache[key] = index

    return index
