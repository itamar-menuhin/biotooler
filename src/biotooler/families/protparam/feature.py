"""ProtParam feature computation using Bio.SeqUtils.ProtParam."""

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.protein_family import ProteinFamilyFeature
from biotooler.core.types import Scalar
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class ProtParamFeature(ProteinFamilyFeature):
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

    def compute_protein(self, protein_record: SeqRecord) -> dict[str, Scalar]:
        """Compute ProtParam features for a protein sequence.

        Args:
            protein_record: Protein SeqRecord to analyze (molecule_type="protein")

        Returns:
            Dictionary mapping feature names to scalar values:
            - molecular_weight: Molecular weight in Daltons
            - isoelectric_point: Theoretical isoelectric point
            - gravy: Grand Average of Hydropathy
            - instability_index: Instability index
            - aromaticity: Fraction of aromatic amino acids

        Raises:
            ValueError: If computation fails
        """
        # Get protein sequence string
        protein_seq = str(protein_record.seq)

        # Lazy import ProteinAnalysis to avoid loading it at module import time
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        # Create ProteinAnalysis object
        try:
            pa = ProteinAnalysis(protein_seq)
        except Exception as e:
            raise ValueError(
                f"Failed to create ProteinAnalysis for record {protein_record.id!r}: {e}"
            ) from e

        # Compute all features
        result: dict[str, Scalar] = {}

        try:
            result["molecular_weight"] = pa.molecular_weight()
        except Exception as e:
            raise ValueError(
                f"Failed to compute molecular_weight for record {protein_record.id!r}: {e}"
            ) from e

        try:
            result["isoelectric_point"] = pa.isoelectric_point()
        except Exception as e:
            raise ValueError(
                f"Failed to compute isoelectric_point for record {protein_record.id!r}: {e}"
            ) from e

        try:
            result["gravy"] = pa.gravy()
        except Exception as e:
            raise ValueError(
                f"Failed to compute gravy for record {protein_record.id!r}: {e}"
            ) from e

        try:
            result["instability_index"] = pa.instability_index()
        except Exception as e:
            raise ValueError(
                f"Failed to compute instability_index for record {protein_record.id!r}: {e}"
            ) from e

        try:
            result["aromaticity"] = pa.aromaticity()
        except Exception as e:
            raise ValueError(
                f"Failed to compute aromaticity for record {protein_record.id!r}: {e}"
            ) from e

        return result

    @property
    def position_space(self) -> PositionSpace:
        """The position space for this feature (RESIDUE).

        Returns:
            PositionSpace.RESIDUE indicating features are computed per amino acid residue
        """
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Mapping of feature keys to aggregation specifications.

        For ProtParam, we compute per-residue values for metrics that can be decomposed:
        - aromaticity: Mean of binary indicators (1 for aromatic AAs: F, W, Y)
        - gravy: Mean of Kyte-Doolittle hydropathy values per residue
        - molecular_weight: Sum of per-residue weights minus water losses

        Note: isoelectric_point and instability_index cannot be meaningfully
        decomposed per-position, so they are not included in vector_keys.

        Returns:
            Dictionary mapping feature names to AggregationSpec objects that define
            how per-residue values should be aggregated into window values.
        """
        # Custom aggregation function for molecular weight that accounts for peptide bonds
        def aggregate_molecular_weight(values: np.ndarray) -> float:
            """Aggregate molecular weights accounting for water loss in peptide bonds.

            For a window of n residues, the total molecular weight is:
            sum(residue_weights) - (n-1) * 18.01528 (water lost in peptide bonds)
            """
            if len(values) == 0:
                return 0.0
            # Sum residue weights and subtract water molecules lost in peptide bonds
            water_weight = 18.01528
            total = np.sum(values) - (len(values) - 1) * water_weight
            return float(total)

        return {
            "aromaticity": AggregationSpec(aggregation_fn=np.mean),
            "gravy": AggregationSpec(aggregation_fn=np.mean),
            "molecular_weight": AggregationSpec(aggregation_fn=aggregate_molecular_weight),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-residue feature values for the entire protein sequence.

        This method computes per-residue values for decomposable ProtParam metrics,
        allowing for efficient sliding window analysis without sequence slicing.

        Args:
            record: Protein SeqRecord to analyze (molecule_type="protein")
            **kwargs: Additional parameters (unused, for protocol compatibility)

        Returns:
            Dictionary mapping feature names to numpy arrays of per-residue values:
            - aromaticity: Binary array (1 for F/W/Y, 0 otherwise)
            - gravy: Kyte-Doolittle hydropathy value per residue
            - molecular_weight: Residue molecular weight (will be adjusted for
              peptide bonds during aggregation)
        """
        from Bio.SeqUtils import ProtParamData
        from Bio.SeqUtils.ProtParam import ProteinAnalysis

        protein_seq = str(record.seq)

        vectors = {}

        # Aromaticity: binary indicator for aromatic amino acids (F, W, Y)
        aromatic_aas = set("FWY")
        aromaticity_vec = np.array(
            [1.0 if aa in aromatic_aas else 0.0 for aa in protein_seq]
        )
        vectors["aromaticity"] = aromaticity_vec

        # GRAVY: Kyte-Doolittle hydropathy values
        # Use 0.0 for unknown amino acids (though ProteinAnalysis would error on them)
        gravy_vec = np.array(
            [ProtParamData.kd.get(aa, 0.0) for aa in protein_seq]
        )
        vectors["gravy"] = gravy_vec

        # Molecular weight: per-residue weights
        # We compute individual residue weights; the aggregation function will
        # handle subtracting water molecules for peptide bonds
        mw_vec = np.array(
            [ProteinAnalysis(aa).molecular_weight() for aa in protein_seq]
        )
        vectors["molecular_weight"] = mw_vec

        return vectors
