"""Metapredict backend for disorder prediction.

This module provides protein intrinsic disorder prediction using metapredict.
"""

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.record import get_molecule_type
from biotooler.core.translation import ensure_protein_record
from biotooler.families.disorder.integration import require_metapredict
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class DisorderProfileMetapredict:
    """Compute per-residue intrinsic disorder probabilities using metapredict.

    This feature implements the PositionalFeature protocol to compute protein
    intrinsic disorder predictions across the full sequence. It uses metapredict's
    neural network model to predict the probability that each amino acid is in
    a disordered region.

    The feature operates in RESIDUE position space (per amino acid) and returns a
    vector of disorder probabilities (DISORDER_P values) for each position in the
    protein sequence.

    Input handling:
    - Protein sequences: Used directly for prediction
    - DNA/RNA sequences: Automatically translated to protein using ensure_protein_record()

    Settings:
    - Uses default metapredict model and parameters
    - Returns disorder probabilities for all amino acids in the sequence
    - Probabilities are in range [0, 1] where higher values indicate higher disorder

    Args:
        table: NCBI genetic code table number for translation (default: 1)
        use_orf_if_present: If True, use attached ORF for translation (default: True)
        strip_terminal_stop: If True, strip terminal stop codon (default: True)
        on_internal_stop: Action for internal stops: "error" or "ignore" (default: "error")

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = DisorderProfileMetapredict()
        >>> # Protein sequence
        >>> record = SeqRecord(Seq("MKALVSWGR"), id="test")
        >>> record.annotations["molecule_type"] = "protein"
        >>> result = feature.compute_vector(record)
        >>> print(result["DISORDER_P"])  # Array of disorder probabilities
        >>> # DNA sequence (will be translated)
        >>> dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        >>> dna.annotations["molecule_type"] = "DNA"
        >>> result = feature.compute_vector(dna)
        >>> len(result["DISORDER_P"])  # Length matches translated protein

    Notes:
        - metapredict is lazily imported only when the feature is used
        - Results are deterministic for the same input sequence
        - Empty sequences return empty arrays
        - Translation follows the same rules as other protein family features
    """

    def __init__(
        self,
        *,
        table: int = 1,
        use_orf_if_present: bool = True,
        strip_terminal_stop: bool = True,
        on_internal_stop: str = "error",
    ):
        """Initialize DisorderProfileMetapredict.

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

    @property
    def position_space(self) -> PositionSpace:
        """Return RESIDUE position space.

        Disorder prediction is computed per amino acid residue.

        Returns:
            PositionSpace.RESIDUE
        """
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specifications for each feature key.

        The DISORDER_P (disorder probability) values are aggregated using mean,
        which gives the average disorder propensity over a window.

        Returns:
            Dictionary mapping "DISORDER_P" to AggregationSpec with np.mean aggregation
        """
        return {
            "DISORDER_P": AggregationSpec(aggregation_fn=np.mean),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-residue disorder probabilities for the entire sequence.

        This method uses metapredict to predict intrinsic disorder probabilities
        for each amino acid position. If the input is DNA/RNA, it is automatically
        translated to protein using ensure_protein_record().

        The computation uses metapredict's predict_disorder() function, which
        returns disorder probabilities for the full protein sequence.

        Args:
            record: SeqRecord containing the DNA, RNA, or protein sequence to analyze
            **kwargs: Additional parameters that may include:
                - positions: Optional array of 0-based position indices. If provided,
                  the feature still computes the full vector (for correct context),
                  and the returned array will have the full length with valid values
                  at all positions. The window engine will extract needed positions
                  from the full vector.

        Returns:
            Dictionary with key "DISORDER_P" mapping to numpy array of disorder
            probabilities. Array length equals the protein sequence length, with
            values in range [0, 1].

        Raises:
            ImportError: If metapredict is not installed
            ValueError: If translation fails or other errors occur

        Notes:
            - Uses default metapredict model and parameters
            - Automatically handles DNA/RNA to protein translation
            - Returns deterministic results for the same input sequence
            - Computation is performed on the full sequence for correct context
            - positions parameter is accepted but full computation is always performed
        """
        # Lazy import metapredict
        metapredict = require_metapredict()

        # Get molecule type to check if translation is needed
        mol_type = get_molecule_type(record)

        # If DNA or RNA, translate to protein
        if mol_type.upper() in ("DNA", "RNA"):
            protein_record = ensure_protein_record(
                record,
                table=self.table,
                use_orf_if_present=self.use_orf_if_present,
                strip_terminal_stop=self.strip_terminal_stop,
                on_internal_stop=self.on_internal_stop,
            )
        else:
            # Already protein, use as-is
            protein_record = record

        # Get protein sequence string
        protein_seq = str(protein_record.seq)

        # Handle empty sequence
        if len(protein_seq) == 0:
            return {"DISORDER_P": np.array([], dtype=np.float64)}

        # Compute disorder probabilities using metapredict
        # predict_disorder returns a numpy array of disorder probabilities
        try:
            disorder_probs = metapredict.predict_disorder(protein_seq)
        except Exception as e:
            raise ValueError(f"Failed to predict disorder for record {record.id!r}: {e}") from e

        # Ensure the result is a numpy array and has the expected length
        disorder_probs = np.array(disorder_probs, dtype=np.float64)
        if len(disorder_probs) != len(protein_seq):
            raise ValueError(
                f"metapredict returned array of length {len(disorder_probs)}, "
                f"but expected {len(protein_seq)} for record {record.id!r}"
            )

        return {"DISORDER_P": disorder_probs}
