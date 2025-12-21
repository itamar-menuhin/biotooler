"""Record utilities for working with Biopython SeqRecord objects."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.errors import InvalidSequenceError


def coerce_record(
    seq: str | Seq | SeqRecord,
    molecule_type: str,
    id: str | None = None,
    description: str | None = None,
    allow_ambiguous: bool = True,
) -> SeqRecord:
    """Coerce a sequence to a SeqRecord with validation and normalization.

    This function accepts various sequence representations (string, Seq, or SeqRecord)
    and returns a normalized SeqRecord object with proper metadata and validation.

    Normalization steps:
    - Strip whitespace from sequence
    - Convert to uppercase
    - Set molecule_type annotation
    - Validate characters based on molecule type

    Args:
        seq: Input sequence as string, Bio.Seq.Seq, or Bio.SeqRecord.SeqRecord
        molecule_type: Must be one of "DNA", "RNA", or "protein"
        id: Optional sequence identifier. If None and seq is a SeqRecord, preserves
            existing id. Otherwise defaults to "<unknown id>".
        description: Optional sequence description. If None and seq is a SeqRecord,
            preserves existing description. Otherwise defaults to "<unknown description>".
        allow_ambiguous: If True, allows ambiguous characters (N for DNA/RNA, X for protein).
            If False, only standard characters are allowed.

    Returns:
        A validated SeqRecord with normalized sequence and proper annotations

    Raises:
        InvalidSequenceError: If the sequence contains invalid characters for the
            specified molecule type
        ValueError: If molecule_type is not one of "DNA", "RNA", or "protein"

    Examples:
        >>> record = coerce_record("ACGT", "DNA")
        >>> record.seq
        Seq('ACGT')
        >>> record.annotations["molecule_type"]
        'DNA'

        >>> record = coerce_record("acgt", "DNA", id="seq1", description="Test sequence")
        >>> str(record.seq)
        'ACGT'
        >>> record.id
        'seq1'

        >>> coerce_record("ACGTN", "DNA", allow_ambiguous=False)
        Traceback (most recent call last):
            ...
        InvalidSequenceError: Invalid character 'N' at position 4 for DNA sequence
    """
    # Validate molecule_type
    if molecule_type not in ("DNA", "RNA", "protein"):
        raise ValueError(f"molecule_type must be 'DNA', 'RNA', or 'protein', got {molecule_type!r}")

    # Extract sequence string and existing metadata
    if isinstance(seq, SeqRecord):
        existing_id = seq.id
        existing_description = seq.description
        seq_str = str(seq.seq)
    elif isinstance(seq, Seq):
        existing_id = None
        existing_description = None
        seq_str = str(seq)
    else:
        existing_id = None
        existing_description = None
        seq_str = str(seq)

    # Normalize sequence: strip whitespace and uppercase
    seq_str = seq_str.replace(" ", "").replace("\t", "").replace("\n", "").replace("\r", "")
    seq_str = seq_str.upper()

    # Validate sequence characters
    _validate_sequence(seq_str, molecule_type, allow_ambiguous)

    # Create new SeqRecord
    record_id = id if id is not None else (existing_id if existing_id else "<unknown id>")
    record_description = (
        description
        if description is not None
        else (existing_description if existing_description else "<unknown description>")
    )
    record = SeqRecord(Seq(seq_str), id=record_id, description=record_description)

    # Set molecule_type annotation
    record.annotations["molecule_type"] = molecule_type

    return record


def get_molecule_type(record: SeqRecord) -> str:
    """Get the molecule type from a SeqRecord's annotations.

    Args:
        record: A SeqRecord object with molecule_type annotation

    Returns:
        The molecule type as a string ("DNA", "RNA", or "protein")

    Raises:
        KeyError: If the record does not have a molecule_type annotation

    Examples:
        >>> record = coerce_record("ACGT", "DNA")
        >>> get_molecule_type(record)
        'DNA'
    """
    try:
        mol_type = record.annotations["molecule_type"]
        if not isinstance(mol_type, str):
            raise TypeError(f"molecule_type annotation must be a string, got {type(mol_type)}")
        return mol_type
    except KeyError as e:
        raise KeyError(
            f"SeqRecord {record.id!r} is missing 'molecule_type' annotation. "
            "Use coerce_record() to ensure proper metadata."
        ) from e


def _validate_sequence(seq: str, molecule_type: str, allow_ambiguous: bool) -> None:
    """Validate sequence characters based on molecule type.

    Args:
        seq: Sequence string to validate (should be uppercase)
        molecule_type: "DNA", "RNA", or "protein"
        allow_ambiguous: Whether to allow ambiguous characters

    Raises:
        InvalidSequenceError: If invalid characters are found
    """
    if molecule_type == "DNA":
        valid_chars = set("ACGT")
        if allow_ambiguous:
            valid_chars.update("NRYWSMKHBVD")  # IUPAC ambiguous nucleotide codes
    elif molecule_type == "RNA":
        valid_chars = set("ACGU")
        if allow_ambiguous:
            valid_chars.update("NRYWSMKHBVD")  # IUPAC ambiguous nucleotide codes
    else:  # protein
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY")
        if allow_ambiguous:
            # X=unknown, B=Asx, Z=Glx, J=Leu/Ile, U=Sec, O=Pyl, *=stop
            valid_chars.update("XBZJUO*")

    for i, char in enumerate(seq):
        if char not in valid_chars:
            raise InvalidSequenceError(
                f"Invalid character {char!r} at position {i} for {molecule_type} sequence",
                char=char,
                index=i,
                molecule_type=molecule_type,
            )
