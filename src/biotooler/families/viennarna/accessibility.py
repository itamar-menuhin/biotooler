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
        >>> print(result["PU"])  # Array of unpaired probabilities (example values)
        # Output will be an array of floats between 0.0 and 1.0

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
            "PU": AggregationSpec(name="MEAN", aggregation_fn=np.mean),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-nucleotide unpaired probabilities for the entire sequence.

        This method uses ViennaRNA's partition function to compute unpaired
        probabilities (accessibility) for each nucleotide position. The sequence
        is automatically normalized from DNA to RNA (T->U) before computation.

        The computation uses ViennaRNA's partition function and plist_from_probs
        to obtain base pairing probabilities, then derives unpaired probabilities
        as 1 - (sum of pairing probabilities with all other positions).

        Args:
            record: SeqRecord containing the DNA or RNA sequence to analyze
            **kwargs: Additional parameters that may include:
                - positions: Optional array of 0-based position indices. If provided,
                  the feature still computes the full vector (for correct structural
                  context), and the returned array will have the full length with
                  valid values at all positions. The window engine will extract
                  needed positions from the full vector.

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
            - positions parameter is accepted but full computation is always performed
        """
        # Lazy import ViennaRNA
        from biotooler.families.viennarna.integration import require_viennarna

        RNA = require_viennarna()

        # Get sequence and normalize DNA to RNA (T -> U)
        seq_str = str(record.seq).upper().replace("T", "U")

        # Handle empty sequence
        if len(seq_str) == 0:
            return {"PU": np.array([], dtype=np.float64)}

        # Compute partition function
        fc = RNA.fold_compound(seq_str)
        fc.pf()  # Compute partition function

        # Get base pairing probabilities as a list of (i, j, p) tuples
        # ViennaRNA uses 1-based indexing for i and j
        # plist_from_probs(cutoff) returns pairs with probability >= cutoff
        # Use cutoff=0.0 to get all pairs
        pairs = fc.plist_from_probs(0.0)

        # Compute unpaired probabilities
        # For each position i, PU[i] = 1 - sum_j(P(i,j))
        # where P(i,j) is the probability that positions i and j are paired
        seq_len = len(seq_str)
        paired_prob = np.zeros(seq_len, dtype=np.float64)

        # Accumulate paired probabilities from the pair list
        # Each pair (i, j, p) contributes probability p to both positions i and j
        for pair in pairs:
            if pair.p > 0:
                # Convert from 1-based to 0-based indexing
                i_idx = pair.i - 1
                j_idx = pair.j - 1
                # Ensure indices are valid
                if 0 <= i_idx < seq_len and 0 <= j_idx < seq_len:
                    paired_prob[i_idx] += pair.p
                    paired_prob[j_idx] += pair.p

        # Compute unpaired probabilities
        pu_values = 1.0 - paired_prob

        return {"PU": pu_values}
