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

# Lazy import chimera at module level - will raise ImportError if not installed
# Note: The package is called "chimera-ugem" but the module is "chimera"
# To avoid conflicts with test packages named "chimera", we check if we got
# the correct module and attempt recovery if needed
try:
    pychimera = lazy_import(  # type: ignore[misc]
        "chimera",
        extra="chimera",
        purpose="computing gene expression features using Chimera algorithms",
    )

    # Verify we got the correct chimera module (not a test package)
    # The real chimera module should have calc_cARS
    if not hasattr(pychimera, "calc_cARS"):
        # We got the wrong module (probably tests/families/chimera)
        # Raise ImportError so tests can skip gracefully
        raise ImportError(
            "Module 'chimera' (pychimera) found but it doesn't have calc_cARS. "
            "This might be a naming conflict with test packages. "
            'Install chimera-ugem (pychimera) with: pip install "biotooler[chimera]"'
        )
except ImportError:
    # Re-raise to preserve original error message
    raise


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
            # cARS/PScARS algorithms work on DNA sequences, not proteins
            # Get CDS strings from the reference set
            try:
                self.reference_seqs = reference_set.cds_strings(
                    require_multiple_of_three=False
                )
                self._sequence_type = "cds"
            except Exception as e:
                raise ValueError(
                    f"reference_set cannot provide CDS sequences: {e}"
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
            record: SeqRecord containing the DNA/RNA sequence to analyze
            **kwargs: Additional parameters (unused)

        Returns:
            Dictionary mapping feature names to numpy arrays of per-position values.
            For cARS, this is the maximal common substring length at each codon position.

        Raises:
            ValueError: If algorithm is not supported
        """
        # Import from the pychimera lazy import to get the correct module
        # (avoids conflicts with local test package names)
        build_suffix_array = pychimera.build_suffix_array
        calc_cARS = pychimera.calc_cARS
        nt2codon = pychimera.nt2codon

        # Only cARS/PScARS are supported for now
        if self.algorithm not in ("cARS", "PScARS"):
            raise NotImplementedError(
                f"Algorithm '{self.algorithm}' is not yet implemented for vector computation. "
                f"Only 'cARS' and 'PScARS' are currently supported."
            )

        # Normalize RNA to DNA (replace U with T)
        target_seq = str(record.seq).upper().replace("U", "T")

        # Normalize reference sequences (RNA to DNA)
        normalized_refs = [ref.upper().replace("U", "T") for ref in self.reference_seqs]

        # Build suffix array from reference sequences (cache if not already built)
        if self._suffix_array is None:
            ref_cod = nt2codon(normalized_refs)
            # Include position-specific data if using PScARS
            self._suffix_array = build_suffix_array(
                ref_cod, pos_spec=("PS" in self.algorithm)
            )

        # Convert target sequence to codon representation
        target_cod = nt2codon([target_seq])

        # Handle empty sequence after codon conversion
        if not target_cod or not target_cod[0]:
            feature_name = f"{self.algorithm}_score"
            return {feature_name: np.array([], dtype=np.float64)}

        # Compute cARS vector (per-position values)
        cars_vec = calc_cARS(
            target_cod[0],
            self._suffix_array,
            win_params=self.win_params if "PS" in self.algorithm else None,
            max_len=self.max_len,
            max_pos=self.max_pos,
            return_vec=True,
        )

        feature_name = f"{self.algorithm}_score"
        return {feature_name: np.array(cars_vec, dtype=np.float64)}

    def __call__(
        self, record: SeqRecord, orf_nt_span: tuple[int, int] | None = None
    ) -> dict[str, Scalar]:
        """Compute Chimera features for the target sequence.

        Args:
            record: DNA/RNA SeqRecord to analyze (for cARS) or protein (for cMap)
            orf_nt_span: Optional tuple (start, end) for extracting subsequence
                in nucleotide coordinates. If provided, extracts record.seq[start:end]
                before computing features.

        Returns:
            Dictionary mapping feature names to scalar values.
            For cARS/PScARS: returns {f"{algorithm}_score": float}

        Raises:
            ValueError: If reference_set was not provided or algorithm not supported
        """
        # Import from the pychimera lazy import to get the correct module
        # (avoids conflicts with local test package names)
        build_suffix_array = pychimera.build_suffix_array
        calc_cARS = pychimera.calc_cARS
        nt2codon = pychimera.nt2codon

        # Only cARS/PScARS are supported for now
        if self.algorithm not in ("cARS", "PScARS"):
            raise NotImplementedError(
                f"Algorithm '{self.algorithm}' is not yet implemented. "
                f"Only 'cARS' and 'PScARS' are currently supported."
            )

        # Extract subsequence if orf_nt_span is provided
        target_seq = str(record.seq)
        if orf_nt_span is not None:
            start, end = orf_nt_span
            target_seq = target_seq[start:end]

        # Normalize RNA to DNA (replace U with T)
        target_seq = target_seq.upper().replace("U", "T")

        # Normalize reference sequences (RNA to DNA)
        normalized_refs = [ref.upper().replace("U", "T") for ref in self.reference_seqs]

        # Build suffix array from reference sequences (cache if not already built)
        if self._suffix_array is None:
            ref_cod = nt2codon(normalized_refs)
            # Include position-specific data if using PScARS
            self._suffix_array = build_suffix_array(
                ref_cod, pos_spec=("PS" in self.algorithm)
            )

        # Convert target sequence to codon representation
        target_cod = nt2codon([target_seq])

        # Handle empty sequence after codon conversion
        if not target_cod or not target_cod[0]:
            feature_name = f"{self.algorithm}_score"
            return {feature_name: 0.0}

        # Compute cARS score
        score = calc_cARS(
            target_cod[0],
            self._suffix_array,
            win_params=self.win_params if "PS" in self.algorithm else None,
            max_len=self.max_len,
            max_pos=self.max_pos,
            return_vec=False,
        )

        feature_name = f"{self.algorithm}_score"
        return {feature_name: float(score)}
