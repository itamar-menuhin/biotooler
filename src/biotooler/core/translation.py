"""Translation utilities for converting DNA/RNA records to protein.

This module provides utilities for ensuring sequences are in protein form,
with automatic translation of DNA/RNA sequences when needed.
"""

from typing import Any, Literal

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.core.orf_store import OrfSpan, get_orf
from biotooler.core.record import get_molecule_type
from biotooler.core.seq_utils import get_seq_str


def ensure_protein_record(
    record: SeqRecord,
    *,
    table: int = 1,
    use_orf_if_present: bool = True,
    orf: OrfSpan | None = None,
    strip_terminal_stop: bool = True,
    on_internal_stop: str = "error",
    orf_policy: Literal["default", "longest_orf"] = "default",
) -> SeqRecord:
    """Ensure a SeqRecord is in protein form, translating if necessary.

    This function checks if the input record is already a protein sequence.
    If it is, it returns the record as-is (or a shallow copy). If it's a
    DNA/RNA sequence, it translates it to protein using the specified options.

    Args:
        record: Input SeqRecord (DNA, RNA, or protein)
        table: NCBI genetic code table number for translation (default: 1, standard code)
        use_orf_if_present: If True and no explicit orf given, use attached orf
            from record.annotations["biotooler.orf"] if present (default: True)
        orf: Explicit ORF region to translate as (start, end) tuple. If provided,
            overrides any attached ORF. If None and use_orf_if_present=True,
            uses attached ORF if present. Otherwise uses full sequence in frame 0.
        strip_terminal_stop: If True, removes terminal stop codon (*) from the
            translated sequence (default: True)
        on_internal_stop: Action to take if internal stop codons are found after
            translation. Must be "error" (raises ValueError) or "ignore" (keeps stops).
            Default: "error"
        orf_policy: ORF selection policy when no explicit or attached ORF is provided.
            Must be "default" or "longest_orf". Default: "default"
            - "default": Use full sequence in frame 0 (preserves existing behavior)
            - "longest_orf": Find ORF candidates and select the longest one. If multiple
              ORFs have the same length, selects the one with the earliest start position.
              If still tied, selects the one with the earliest stop position. Falls back
              to frame 0 if no ORF candidates are found.

    Returns:
        SeqRecord containing protein sequence with preserved id and description.
        If translation was performed, adds annotations about translation details.

    Raises:
        ValueError: If on_internal_stop="error" and internal stop codons are found,
            or if orf_policy has an invalid value
        KeyError: If molecule_type annotation is missing from record

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> # Protein record passes through unchanged
        >>> protein = SeqRecord(Seq("MKALV"), id="test")
        >>> protein.annotations["molecule_type"] = "protein"
        >>> result = ensure_protein_record(protein)
        >>> str(result.seq)
        'MKALV'
        >>> # DNA record is translated (default policy uses frame 0)
        >>> dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        >>> dna.annotations["molecule_type"] = "DNA"
        >>> result = ensure_protein_record(dna)
        >>> str(result.seq)
        'MKALV'
        >>> result.annotations["translation_performed"]
        True
        >>> result.annotations["translation_table"]
        1
        >>> # DNA record with longest_orf policy
        >>> dna = SeqRecord(Seq("NNNNATGAAAGCCCTGGTGTAA"), id="test")
        >>> dna.annotations["molecule_type"] = "DNA"
        >>> result = ensure_protein_record(dna, orf_policy="longest_orf")
        >>> str(result.seq)
        'MKALV'
        >>> result.annotations["translation_region_source"]
        'longest_orf'
    """
    # Validate on_internal_stop parameter
    if on_internal_stop not in ("error", "ignore"):
        raise ValueError(f"on_internal_stop must be 'error' or 'ignore', got {on_internal_stop!r}")

    # Validate orf_policy parameter
    if orf_policy not in ("default", "longest_orf"):
        raise ValueError(f"orf_policy must be 'default' or 'longest_orf', got {orf_policy!r}")

    # Get molecule type
    mol_type = get_molecule_type(record)

    # If already protein, return as-is (or shallow copy)
    if mol_type.upper() == "PROTEIN":
        # Create a shallow copy to avoid modifying the input
        return SeqRecord(
            record.seq,
            id=record.id,
            name=record.name,
            description=record.description,
            annotations=dict(record.annotations),
        )

    # DNA or RNA - needs translation
    if mol_type.upper() not in ("DNA", "RNA"):
        raise ValueError(
            f"Record molecule_type must be 'DNA', 'RNA', or 'protein', got {mol_type!r}"
        )

    # Determine region to translate
    region_start = 0
    region_end = len(record.seq) if record.seq else 0
    region_source = "full_sequence_frame0"

    if orf is not None:
        # Explicit ORF takes priority
        region_start, region_end = orf
        region_source = "explicit_orf"
    elif use_orf_if_present:
        # Check for attached ORF
        try:
            attached_orf = get_orf(record)
            region_start, region_end = attached_orf
            region_source = "attached_orf"
        except KeyError:
            # No attached ORF, check orf_policy
            if orf_policy == "longest_orf":
                # Find ORF candidates and select the longest one
                candidates = find_orf_candidates(record)
                if candidates:
                    # Select longest ORF with deterministic tie-breaking
                    # Candidates are already sorted by (start, end)
                    # Find the longest by length, then earliest start, then earliest stop
                    longest_orf = max(
                        candidates, key=lambda span: (span[1] - span[0], -span[0], -span[1])
                    )
                    region_start, region_end = longest_orf
                    region_source = "longest_orf"
                # else: no candidates, keep frame0 default
    elif orf_policy == "longest_orf":
        # use_orf_if_present=False but orf_policy="longest_orf" specified
        # Find ORF candidates and select the longest one
        candidates = find_orf_candidates(record)
        if candidates:
            longest_orf = max(candidates, key=lambda span: (span[1] - span[0], -span[0], -span[1]))
            region_start, region_end = longest_orf
            region_source = "longest_orf"
        # else: no candidates, keep frame0 default

    # Get sequence string and extract region
    seq_str = get_seq_str(record)
    region_seq = seq_str[region_start:region_end]

    # Normalize RNA U->T for translation (Biopython expects DNA for translate)
    if mol_type.upper() == "RNA":
        region_seq = region_seq.replace("U", "T")

    # Translate using Biopython
    # Use to_stop=False to translate through stops
    # Use cds=False to not require start/stop codons
    try:
        translated = Seq(region_seq).translate(
            table=table,
            to_stop=False,
            cds=False,  # type: ignore[arg-type]
        )
    except Exception as e:
        raise ValueError(
            f"Translation failed for record {record.id!r} using table {table}: {e}"
        ) from e

    # Convert to string for processing
    protein_str = str(translated)

    # Check for internal stops (before stripping terminal stop)
    # Internal stops are any stops except the terminal one
    has_terminal_stop = protein_str.endswith("*")
    if has_terminal_stop:
        # Check if there are stops before the terminal one
        internal_part = protein_str[:-1]
        if "*" in internal_part:
            if on_internal_stop == "error":
                stop_positions = [i for i, char in enumerate(internal_part) if char == "*"]
                raise ValueError(
                    f"Internal stop codon(s) found at position(s) {stop_positions} "
                    f"in translated protein from record {record.id!r}. "
                    f"Use on_internal_stop='ignore' to keep internal stops."
                )
            # else: on_internal_stop == "ignore", keep the stops
    else:
        # No terminal stop, so all stops are internal
        if "*" in protein_str:
            if on_internal_stop == "error":
                stop_positions = [i for i, char in enumerate(protein_str) if char == "*"]
                raise ValueError(
                    f"Internal stop codon(s) found at position(s) {stop_positions} "
                    f"in translated protein from record {record.id!r}. "
                    f"Use on_internal_stop='ignore' to keep internal stops."
                )
            # else: on_internal_stop == "ignore", keep the stops

    # Strip terminal stop if requested
    if strip_terminal_stop and has_terminal_stop:
        protein_str = protein_str[:-1]

    # Create new protein record
    annotations: dict[str, Any] = {
        "molecule_type": "protein",
        "translation_performed": True,
        "translation_table": table,
        "translation_source": mol_type,
        "translation_region": (region_start, region_end),
        "translation_region_source": region_source,
    }

    protein_record = SeqRecord(
        Seq(protein_str),
        id=record.id,
        name=record.name,
        description=record.description,
        annotations=annotations,  # type: ignore[arg-type]
    )

    return protein_record
