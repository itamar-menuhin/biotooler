"""ViennaRNA accessibility feature for RNA unpaired probabilities.

This module computes per-nucleotide unpaired probabilities (accessibility) using
ViennaRNA's partition function calculations. The accessibility feature measures
how likely each nucleotide is to be unpaired in the ensemble of RNA structures.
"""

from typing import TYPE_CHECKING

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.features.aggregation import AggregationSpec, PositionSpace

if TYPE_CHECKING:
    pass


class ViennaRNAAccessibility:
    """Compute per-nucleotide unpaired probabilities using ViennaRNA.

    This feature implements the PositionalFeature protocol to compute RNA accessibility
    (unpaired probabilities) across the full sequence. It uses ViennaRNA's partition
    function to calculate the probability that each nucleotide is unpaired in the
    ensemble of possible RNA secondary structures.

    The feature operates in RESIDUE position space (per-nucleotide) and returns a
    vector of unpaired probabilities (PU values) for each position in the sequence.

    Settings:
        - Uses default ViennaRNA temperature (37°C)
        - Computes partition function with default energy parameters
        - Returns unpaired probabilities for all nucleotides in the sequence

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = ViennaRNAAccessibility()
        >>> record = SeqRecord(Seq("ACGUACGU"), id="test")
        >>> result = feature.compute_vector(record)
        >>> print(result["PU"])  # Array of unpaired probabilities
        [0.95, 0.82, 0.45, 0.38, 0.45, 0.82, 0.95, 0.95]

    Notes:
        - Input sequences are automatically normalized from DNA to RNA (T->U)
        - ViennaRNA is lazily imported only when the feature is used
        - The partition function is computed once per sequence
        - Results are deterministic for the same input sequence
    """

    @property
    def position_space(self) -> PositionSpace:
        """Return RESIDUE position space.

        ViennaRNA accessibility is computed per nucleotide (residue).

        Returns:
            PositionSpace.RESIDUE
        """
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specifications for each feature key.

        The PU (unpaired probability) values are aggregated using mean,
        which gives the average accessibility over a window.

        Returns:
            Dictionary mapping "PU" to AggregationSpec with np.mean aggregation
        """
        return {
            "PU": AggregationSpec(aggregation_fn=np.mean),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-nucleotide unpaired probabilities for the entire sequence.

        This method uses ViennaRNA's partition function to compute unpaired
        probabilities (accessibility) for each nucleotide position. The sequence
        is automatically normalized from DNA to RNA (T->U) before computation.

        The computation uses ViennaRNA's pf_fold function to calculate the partition
        function and base pairing probabilities, then derives unpaired probabilities
        as 1 - (sum of pairing probabilities with all other positions).

        Args:
            record: SeqRecord containing the DNA or RNA sequence to analyze
            **kwargs: Additional parameters (currently unused, reserved for future)

        Returns:
            Dictionary with key "PU" mapping to numpy array of unpaired probabilities.
            Array length equals the sequence length, with values in range [0, 1].

        Raises:
            ImportError: If ViennaRNA is not installed

        Notes:
            - Uses default ViennaRNA settings (37°C, default energy parameters)
            - Automatically converts DNA (with T) to RNA (with U)
            - Returns deterministic results for the same input sequence
            - Computation is performed on the full sequence for correct context
        """
        # Lazy import ViennaRNA
        from biotooler.families.viennarna.integration import require_viennarna

        RNA = require_viennarna()

        # Get sequence and normalize DNA to RNA (T -> U)
        seq_str = str(record.seq).upper().replace("T", "U")

        # Handle empty sequence
        if len(seq_str) == 0:
            return {"PU": np.array([], dtype=np.float64)}

        # Compute partition function and base pairing probabilities
        # This uses ViennaRNA's RNA.pf_fold to get the full partition function
        fc = RNA.fold_compound(seq_str)
        fc.pf()  # Compute partition function

        # Get base pairing probability matrix
        # bpp is a tuple of tuples where bpp[i][j] is the probability of pairing between i and j
        # ViennaRNA uses 1-based indexing, so bpp[0] is padding
        bpp = fc.bpp()

        # Compute unpaired probabilities
        # For each position i, PU[i] = 1 - sum_j(P(i,j))
        # where P(i,j) is the probability that positions i and j are paired
        seq_len = len(seq_str)
        pu_values = np.ones(seq_len, dtype=np.float64)

        # Iterate through positions
        # bpp uses 1-based indexing: bpp[1..n] for n nucleotides
        for i in range(1, seq_len + 1):
            # Sum all pairing probabilities for position i
            # bpp[i][j] gives the probability that i pairs with j
            paired_prob = 0.0
            if i < len(bpp):
                for j in range(len(bpp[i])):
                    paired_prob += bpp[i][j]
            # Unpaired probability is 1 - sum of pairing probabilities
            pu_values[i - 1] = 1.0 - paired_prob

        return {"PU": pu_values}
