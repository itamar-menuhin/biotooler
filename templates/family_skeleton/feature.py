"""{{FAMILY_NAME_TITLE}} feature implementation.

{{FAMILY_DESCRIPTION}}
"""

from Bio.SeqRecord import SeqRecord
import numpy as np

from biotooler.core.types import Scalar
from biotooler.features.aggregation import AggregationSpec, PositionSpace

# TODO: If you have optional dependencies, use lazy_import:
# from biotooler.core.lazy_import import lazy_import
#
# optional_package = lazy_import(
#     "package_name",
#     extra="{{FAMILY_NAME}}",
#     purpose="computing {{FAMILY_NAME}} features"
# )


class {{FAMILY_NAME_TITLE}}Feature:
    """Feature that computes {{FAMILY_NAME}} metrics.

    TODO: Add detailed class documentation explaining:
    - What this feature computes
    - What algorithms or methods are used
    - Expected input formats
    - Output format and interpretation

    Args:
        TODO: Document initialization parameters

    Examples:
        TODO: Add usage examples
        >>> # Example usage here
        >>> pass
    """

    def __init__(self):
        """Initialize {{FAMILY_NAME_TITLE}}Feature.

        TODO: Add initialization parameters and logic
        """
        pass

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute {{FAMILY_NAME}} features for the entire sequence.

        Args:
            record: DNA, RNA, or protein SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to scalar values

        Raises:
            TODO: Document exceptions that may be raised
        """
        # TODO: Implement feature computation
        return {}

    def init_state(
        self,
        record: SeqRecord,
        *,
        orf: tuple[int, int],
        window_start: int,
        window_end: int,
        **kwargs,
    ) -> dict:
        """Initialize state for the first window (optional, for incremental computation).

        This method is optional. Implement it if your feature supports incremental
        computation over sliding windows for better performance.

        Args:
            record: The SeqRecord containing the sequence
            orf: Tuple of (start, end) coordinates for the ORF
            window_start: Start position of window relative to ORF start
            window_end: End position of window relative to ORF start
            **kwargs: Additional parameters

        Returns:
            State dictionary for tracking incremental computation
        """
        # TODO: Implement state initialization if needed
        return {}

    def step_state(
        self,
        state: dict,
        *,
        out_start: int,
        out_end: int,
        in_start: int,
        in_end: int,
        **kwargs,
    ) -> None:
        """Update state for the next window (optional, for incremental computation).

        This method is optional. Implement it if your feature supports incremental
        computation over sliding windows.

        Args:
            state: The state object to update
            out_start: Start position (relative to ORF) of bases leaving the window
            out_end: End position (relative to ORF) of bases leaving the window
            in_start: Start position (relative to ORF) of bases entering the window
            in_end: End position (relative to ORF) of bases entering the window
            **kwargs: Additional parameters
        """
        # TODO: Implement state update if needed
        pass

    def emit(self, state: dict) -> dict[str, Scalar]:
        """Emit feature values from current state (optional, for incremental computation).

        This method is optional. Implement it if your feature supports incremental
        computation over sliding windows.

        Args:
            state: The state object

        Returns:
            Dictionary mapping feature names to scalar values
        """
        # TODO: Implement feature emission from state if needed
        return {}

    # ========================================================================
    # Positional Feature Protocol (NEW WINDOWING SEMANTICS)
    # ========================================================================
    # The following methods implement the PositionalFeature protocol, which
    # enforces full-context vector computation followed by window aggregation.
    # This is the RECOMMENDED approach for new families.
    #
    # To implement this protocol, you must provide:
    # 1. position_space property: defines granularity (RESIDUE or CODON)
    # 2. vector_keys property: maps feature keys to aggregation specs
    # 3. compute_vector method: computes per-position values across full sequence
    #
    # See docs/dev/adding_a_family.md for detailed guidance.
    # ========================================================================

    @property
    def position_space(self) -> PositionSpace:
        """Define the position space for this feature (REQUIRED for positional features).

        Returns:
            PositionSpace.RESIDUE for per-nucleotide/amino-acid features
            PositionSpace.CODON for per-codon features

        Example:
            return PositionSpace.RESIDUE
        """
        # TODO: Uncomment and return the appropriate position space
        # return PositionSpace.RESIDUE  # or PositionSpace.CODON
        raise NotImplementedError(
            "position_space must be implemented for positional features. "
            "See docs/dev/adding_a_family.md for guidance."
        )

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Define aggregation specs for each feature key (REQUIRED for positional features).

        Returns:
            Dictionary mapping feature names to AggregationSpec objects.
            Each spec defines how per-position values are aggregated into windows.

        Example:
            return {
                "gc_mean": AggregationSpec(aggregation_fn=np.mean),
                "gc_sum": AggregationSpec(aggregation_fn=np.sum),
            }
        """
        # TODO: Uncomment and define your feature keys with aggregation functions
        # return {
        #     "feature_key": AggregationSpec(aggregation_fn=np.mean),
        # }
        raise NotImplementedError(
            "vector_keys must be implemented for positional features. "
            "See docs/dev/adding_a_family.md for guidance."
        )

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-position feature values (REQUIRED for positional features).

        This method computes feature values at each position across the FULL sequence.
        The windowing framework will then slice and aggregate these values into windows.

        IMPORTANT: When wrapping upstream libraries, ensure you:
        1. Call the library with the FULL sequence (not individual windows)
        2. Transform library output into per-position numpy arrays
        3. Return arrays that match the position_space granularity

        Args:
            record: The SeqRecord containing the full sequence
            **kwargs: Additional parameters (e.g., orf coordinates)

        Returns:
            Dictionary mapping feature keys (from vector_keys) to numpy arrays.
            Array length must match sequence length in the position_space.

        Example:
            seq = str(record.seq).upper()
            gc_vector = np.array([1.0 if base in "GC" else 0.0 for base in seq])
            return {"gc_mean": gc_vector, "gc_sum": gc_vector}

        See Also:
            - docs/dev/adding_a_family.md: How to use upstream positional APIs
            - README.md: "Windowing correctness" section for testing guidance
        """
        # TODO: Implement per-position computation across full sequence
        # When wrapping upstream libraries, call them with the full sequence here
        raise NotImplementedError(
            "compute_vector must be implemented for positional features. "
            "See docs/dev/adding_a_family.md for guidance on wrapping upstream libraries."
        )
