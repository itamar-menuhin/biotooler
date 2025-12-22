"""ProtParam feature computation using Bio.SeqUtils.ProtParam."""

from Bio.SeqRecord import SeqRecord

from biotooler.core.translation import ensure_protein_record
from biotooler.core.types import Scalar


class ProtParamFeature:
    """Feature that computes protein physicochemical parameters.

    This feature uses Bio.SeqUtils.ProtParam.ProteinAnalysis to compute:
    - molecular_weight: Molecular weight in Daltons
    - isoelectric_point: Theoretical isoelectric point (pI)
    - gravy: Grand Average of Hydropathy (GRAVY)
    - instability_index: Instability index (>40 indicates unstable protein)
    - aromaticity: Aromaticity (fraction of aromatic amino acids)

    The feature automatically translates DNA/RNA sequences to protein using
    ensure_protein_record() before computing parameters.

    Args:
        table: NCBI genetic code table number for translation (default: 1)
        use_orf_if_present: If True, use attached ORF for translation (default: True)
        strip_terminal_stop: If True, strip terminal stop codon (default: True)
        on_internal_stop: Action for internal stops: "error" or "ignore" (default: "error")

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = ProtParamFeature()
        >>> # Protein sequence
        >>> record = SeqRecord(Seq("MKALVSWGR"), id="test")
        >>> record.annotations["molecule_type"] = "protein"
        >>> result = feature(record)
        >>> "molecular_weight" in result
        True
        >>> "isoelectric_point" in result
        True
        >>> # DNA sequence (will be translated)
        >>> dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        >>> dna.annotations["molecule_type"] = "DNA"
        >>> result = feature(dna)
        >>> "molecular_weight" in result
        True
    """

    def __init__(
        self,
        *,
        table: int = 1,
        use_orf_if_present: bool = True,
        strip_terminal_stop: bool = True,
        on_internal_stop: str = "error",
    ):
        """Initialize ProtParamFeature.

        Args:
            table: NCBI genetic code table for translation (default: 1)
            use_orf_if_present: Use attached ORF if present (default: True)
            strip_terminal_stop: Strip terminal stop codon (default: True)
            on_internal_stop: "error" or "ignore" for internal stops (default: "error")
        """
        self.table = table
        self.use_orf_if_present = use_orf_if_present
        self.strip_terminal_stop = strip_terminal_stop
        self.on_internal_stop = on_internal_stop

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute ProtParam features for the sequence.

        Args:
            record: DNA, RNA, or protein SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to scalar values:
            - molecular_weight: Molecular weight in Daltons
            - isoelectric_point: Theoretical isoelectric point
            - gravy: Grand Average of Hydropathy
            - instability_index: Instability index
            - aromaticity: Fraction of aromatic amino acids

        Raises:
            ValueError: If translation fails or internal stops found (when configured)
        """
        # Ensure we have a protein sequence
        protein_record = ensure_protein_record(
            record,
            table=self.table,
            use_orf_if_present=self.use_orf_if_present,
            strip_terminal_stop=self.strip_terminal_stop,
            on_internal_stop=self.on_internal_stop,
        )

        # Get protein sequence string
        protein_seq = str(protein_record.seq)

        # Lazy import ProteinAnalysis to avoid loading it at module import time
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        # Create ProteinAnalysis object
        try:
            pa = ProteinAnalysis(protein_seq)
        except Exception as e:
            raise ValueError(
                f"Failed to create ProteinAnalysis for record {record.id!r}: {e}"
            ) from e

        # Compute all features
        result: dict[str, Scalar] = {}

        try:
            result["molecular_weight"] = pa.molecular_weight()
        except Exception as e:
            raise ValueError(
                f"Failed to compute molecular_weight for record {record.id!r}: {e}"
            ) from e

        try:
            result["isoelectric_point"] = pa.isoelectric_point()
        except Exception as e:
            raise ValueError(
                f"Failed to compute isoelectric_point for record {record.id!r}: {e}"
            ) from e

        try:
            result["gravy"] = pa.gravy()
        except Exception as e:
            raise ValueError(f"Failed to compute gravy for record {record.id!r}: {e}") from e

        try:
            result["instability_index"] = pa.instability_index()
        except Exception as e:
            raise ValueError(
                f"Failed to compute instability_index for record {record.id!r}: {e}"
            ) from e

        try:
            result["aromaticity"] = pa.aromaticity()
        except Exception as e:
            raise ValueError(f"Failed to compute aromaticity for record {record.id!r}: {e}") from e

        return result
