"""Base class for features that operate on protein sequences.

This module provides a base class for features that need to work with protein
sequences, automatically handling DNA/RNA to protein translation when needed.
"""

from abc import ABC, abstractmethod

from Bio.SeqRecord import SeqRecord

from biotooler.core.translation import ensure_protein_record
from biotooler.core.types import Scalar


class ProteinFamilyFeature(ABC):
    """Base class for features that operate on protein sequences.

    This class provides a standardized adapter pattern for features that need
    to work with protein sequences. It automatically handles translation of
    DNA/RNA sequences to protein using ensure_protein_record(), then delegates
    to the subclass's compute_protein() method.

    Subclasses should:
    1. Store translation parameters (table, use_orf_if_present, etc.) in __init__
    2. Implement compute_protein(protein_record, ...) to compute features on protein

    The base class handles:
    - Calling ensure_protein_record() with the stored translation parameters
    - Passing the translated protein record to compute_protein()

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> from biotooler.core.protein_family import ProteinFamilyFeature
        >>> from biotooler.core.types import Scalar
        >>>
        >>> class MyProteinFeature(ProteinFamilyFeature):
        ...     def __init__(self, table=1):
        ...         self.table = table
        ...         self.use_orf_if_present = True
        ...         self.strip_terminal_stop = True
        ...         self.on_internal_stop = "error"
        ...
        ...     def compute_protein(self, protein_record: SeqRecord) -> dict[str, Scalar]:
        ...         return {"length": len(protein_record.seq)}
        >>>
        >>> feature = MyProteinFeature()
        >>> # Works with protein input
        >>> protein = SeqRecord(Seq("MKALV"), id="test")
        >>> protein.annotations["molecule_type"] = "protein"
        >>> feature(protein)
        {'length': 5}
        >>> # Works with DNA input (automatically translated)
        >>> dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        >>> dna.annotations["molecule_type"] = "DNA"
        >>> feature(dna)
        {'length': 5}
    """

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute features for the sequence, translating to protein if needed.

        This method handles the translation of DNA/RNA sequences to protein
        using ensure_protein_record(), then delegates to compute_protein().

        Args:
            record: DNA, RNA, or protein SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to scalar values

        Raises:
            ValueError: If translation fails or other errors occur
            AttributeError: If required translation parameters are not set
        """
        # Get translation parameters from subclass instance
        # These should be set in the subclass's __init__
        table = getattr(self, "table", 1)
        use_orf_if_present = getattr(self, "use_orf_if_present", True)
        strip_terminal_stop = getattr(self, "strip_terminal_stop", True)
        on_internal_stop = getattr(self, "on_internal_stop", "error")

        # Ensure we have a protein sequence
        protein_record = ensure_protein_record(
            record,
            table=table,
            use_orf_if_present=use_orf_if_present,
            strip_terminal_stop=strip_terminal_stop,
            on_internal_stop=on_internal_stop,
        )

        # Delegate to subclass implementation
        return self.compute_protein(protein_record)

    @abstractmethod
    def compute_protein(self, protein_record: SeqRecord) -> dict[str, Scalar]:
        """Compute features on a protein sequence.

        Subclasses must implement this method to compute features on the
        provided protein sequence. The input is guaranteed to be a protein
        sequence (translation has already been performed if needed).

        Args:
            protein_record: Protein SeqRecord to analyze (molecule_type="protein")

        Returns:
            Dictionary mapping feature names to scalar values

        Raises:
            Subclass-specific exceptions as appropriate
        """
        pass
