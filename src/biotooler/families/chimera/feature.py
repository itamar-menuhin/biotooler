"""ChimeraFeature implementation for gene expression prediction.

Gene expression prediction and sequence optimization using Chimera algorithms.
"""

from typing import TYPE_CHECKING

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.lazy_import import lazy_import
from biotooler.core.types import Scalar
from biotooler.features.aggregation import AggregationSpec, PositionSpace

if TYPE_CHECKING:
    from biotooler.core.reference_sequences import ReferenceSequenceSet

# Lazy import pyChimera at module level - will raise ImportError if not installed
pychimera = lazy_import(  # type: ignore[misc]
    "pychimera",
    extra="chimera",
    purpose="computing gene expression features using Chimera algorithms",
)


class ChimeraFeature:
    """Feature that computes gene expression predictions using Chimera algorithms.

    This feature wraps the pyChimera package to compute adaptation scores (cARS/PScARS)
    or optimize sequences (cMap/PScMap) using the Chimera algorithms. These algorithms
    predict gene expression by measuring similarity to reference gene sets.

    The pyChimera dependency is lazily loaded only when the feature module is imported
    to keep family-level imports lightweight.

    This feature implements the PositionalFeature protocol, which means it computes
    per-position values across the full sequence and then aggregates them into windows.
    This is the correct way to compute cARS/PScARS for windowed sequences, as slicing
    before computation would incorrectly change the maximal common subsequences.

    Args:
        reference_seqs: List of reference DNA sequences (host genes). Deprecated in favor
            of reference_set parameter.
        reference_set: Optional ReferenceSequenceSet containing reference sequences.
            If provided, will use protein_strings() if available, otherwise cds_strings().
        algorithm: Which algorithm to use ("cARS", "PScARS", "cMap", "PScMap")
        max_len: Maximum substring length for homolog filtering (default: 40)
        max_pos: Maximum position difference for position-specific algorithms (default: 0.5)
        win_params: Window parameters for position-specific algorithms (optional)

    Examples:
        >>> from biotooler.families.chimera import ChimeraFeature
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> # Create feature with reference genes
        >>> ref_genes = ["ATGAAACGC...", "ATGGCAAGC..."]
        >>> feature = ChimeraFeature(ref_genes, algorithm="cARS")
        >>> # Compute on a target gene
        >>> record = SeqRecord(Seq("ATGACCGTT..."), id="test")
        >>> result = feature(record)
        >>> print(result)  # {'cARS_score': 35.2}
    """

    def __init__(
        self,
        reference_seqs: list[str] | None = None,
        *,
        reference_set: "ReferenceSequenceSet | None" = None,
        algorithm: str = "cARS",
        max_len: int = 40,
        max_pos: float = 0.5,
        win_params: dict | None = None,
    ):
        """Initialize ChimeraFeature.

        Args:
            reference_seqs: List of reference DNA sequences (deprecated, use reference_set)
            reference_set: Optional ReferenceSequenceSet for reference sequences
            algorithm: Algorithm to use ("cARS", "PScARS", "cMap", "PScMap")
            max_len: Maximum substring length for homolog filtering
            max_pos: Maximum position difference for position-specific algorithms
            win_params: Window parameters dict (for position-specific algorithms)

        Raises:
            ValueError: If neither reference_seqs nor reference_set is provided,
                or if reference_set cannot provide required sequence type
        """
        # Store reference to pychimera module for future use
        self._chimera = pychimera

        # Handle reference sequences
        if reference_set is not None:
            # Validate that reference_set can provide the needed sequences
            # Prefer protein_strings, but accept cds_strings if proteins can be derived
            try:
                # Try to get protein strings (will derive from CDS if needed)
                self.reference_seqs = reference_set.protein_strings(
                    strip_terminal_stop=True, error_on_internal_stop=False
                )
                self._sequence_type = "protein"
            except Exception:
                # Fall back to CDS strings
                try:
                    self.reference_seqs = reference_set.cds_strings(
                        require_multiple_of_three=False
                    )
                    self._sequence_type = "cds"
                except Exception as e:
                    raise ValueError(
                        f"reference_set cannot provide required sequences: {e}"
                    ) from e
        elif reference_seqs is not None:
            self.reference_seqs = reference_seqs
            self._sequence_type = "cds"
        else:
            raise ValueError(
                "Either reference_seqs or reference_set must be provided"
            )

        self.algorithm = algorithm
        self.max_len = max_len
        self.max_pos = max_pos
        self.win_params = win_params

        # Initialize suffix array storage
        self._suffix_array = None

    @property
    def position_space(self) -> PositionSpace:
        """Return the position space for this feature.

        Chimera features compute per-codon values, as cARS operates on
        codon-level sequences.

        Returns:
            PositionSpace.CODON
        """
        return PositionSpace.CODON

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specifications for each feature key.

        For cARS/PScARS, we return per-position maximal common substring lengths
        that should be aggregated using mean (to match the scalar cARS definition).

        Returns:
            Dictionary mapping feature names to AggregationSpec objects
        """
        feature_name = f"{self.algorithm}_score"
        return {
            feature_name: AggregationSpec(aggregation_fn=np.mean),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-position cARS values for the entire sequence.

        This method computes the cARS vector (maximal common substring length at each
        position) across the full sequence. This is the correct way to handle windowing
        for cARS, as slicing before computation would incorrectly change the maximal
        common subsequences.

        Args:
            record: SeqRecord containing the DNA sequence to analyze
            **kwargs: Additional parameters (unused)

        Returns:
            Dictionary mapping feature names to numpy arrays of per-position values.
            For cARS, this is the maximal common substring length at each codon position.

        Raises:
            NotImplementedError: Feature computation not yet implemented
        """
        # TODO: Implement vector computation using self._chimera.calc_cARS with return_vec=True
        # Example implementation:
        # from chimera import build_suffix_array, calc_cARS, nt2codon
        # if self._suffix_array is None:
        #     ref_cod = nt2codon(self.reference_seqs)
        #     self._suffix_array = build_suffix_array(ref_cod)
        # target_cod = nt2codon([str(record.seq)])
        # cars_vec = calc_cARS(target_cod[0], self._suffix_array,
        #                      win_params=self.win_params if "PS" in self.algorithm else None,
        #                      max_len=self.max_len, max_pos=self.max_pos,
        #                      return_vec=True)
        # feature_name = f"{self.algorithm}_score"
        # return {feature_name: cars_vec}

        raise NotImplementedError(
            "ChimeraFeature compute_vector is not yet implemented. "
            "This is a stub for future implementation using pyChimera's return_vec=True API."
        )

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute Chimera features for the target sequence.

        Args:
            record: DNA SeqRecord to analyze (for cARS) or protein (for cMap)

        Returns:
            Dictionary mapping feature names to scalar values
            Currently returns empty dict (stub implementation)

        Raises:
            NotImplementedError: Feature computation not yet implemented
        """
        # TODO: Implement feature computation using self._chimera
        # Example implementation:
        # if self.algorithm in ("cARS", "PScARS"):
        #     from chimera import build_suffix_array, calc_cARS, nt2codon
        #     ref_cod = nt2codon(self.reference_seqs)
        #     SA = build_suffix_array(ref_cod)
        #     target_cod = nt2codon([str(record.seq)])
        #     score = calc_cARS(target_cod[0], SA,
        #                      win_params=self.win_params if "PS" in self.algorithm else None,
        #                      max_len=self.max_len, max_pos=self.max_pos)
        #     return {f"{self.algorithm}_score": score}

        raise NotImplementedError(
            "ChimeraFeature computation is not yet implemented. "
            "This is a stub for future implementation."
        )
