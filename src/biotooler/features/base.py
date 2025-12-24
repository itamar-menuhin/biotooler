"""Base interfaces for feature computation."""

from __future__ import annotations

from typing import TYPE_CHECKING, Protocol

from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar
from biotooler.features.aggregation import AggregationSpec, PositionSpace

if TYPE_CHECKING:
    import numpy as np


class IncrementalFeature(Protocol):
    """Protocol for features that support incremental window computation.

    Features implementing this protocol can compute windows incrementally by:
    1. Initializing state for the first window (init_state)
    2. Updating state for subsequent windows (step_state)
    3. Emitting feature values from state (emit)

    This allows for significant performance improvements when computing many
    overlapping windows, as shared computations can be reused.

    Example:
        >>> class SumFeature:
        ...     def init_state(self, record, *, orf, window_start, window_end, **kwargs):
        ...         # Initialize state with sum of first window
        ...         seq_str = str(record.seq)[orf[0] + window_start:orf[0] + window_end]
        ...         return {"sum": sum(ord(c) for c in seq_str), "window_start": window_start}
        ...
        ...     def step_state(self, state, *, out_start, out_end, in_start, in_end, **kwargs):
        ...         # Update state by removing outgoing bases and adding incoming bases
        ...         # (implementation would go here)
        ...         pass
        ...
        ...     def emit(self, state):
        ...         return {"sum": state["sum"]}
    """

    def init_state(
        self,
        record: SeqRecord,
        *,
        orf: tuple[int, int],
        window_start: int,
        window_end: int,
        **kwargs,
    ) -> object:
        """Initialize state for the first window.

        Args:
            record: The SeqRecord containing the sequence
            orf: Tuple of (start, end) coordinates for the ORF
            window_start: Start position of window relative to ORF start
            window_end: End position of window relative to ORF start
            **kwargs: Additional parameters that may be needed

        Returns:
            State object (implementation-defined) for this window
        """
        ...

    def step_state(
        self,
        state: object,
        *,
        out_start: int,
        out_end: int,
        in_start: int,
        in_end: int,
        **kwargs,
    ) -> None:
        """Update state for the next window (in-place).

        Args:
            state: The state object to update (from init_state or previous step_state)
            out_start: Start position (relative to ORF) of bases leaving the window
            out_end: End position (relative to ORF) of bases leaving the window
            in_start: Start position (relative to ORF) of bases entering the window
            in_end: End position (relative to ORF) of bases entering the window
            **kwargs: Additional parameters that may be needed

        Note:
            This method modifies state in-place and returns None.
        """
        ...

    def emit(self, state: object) -> dict[str, Scalar]:
        """Emit feature values from current state.

        Args:
            state: The state object (from init_state or after step_state)

        Returns:
            Dictionary mapping feature names to scalar values
        """
        ...


class PositionalFeature(Protocol):
    """Protocol for features that compute per-position values before window aggregation.

    The new windowing semantics:
    1. Compute per-position values across the full sequence
    2. Aggregate per-position values into windows using specified aggregation functions

    This protocol requires:
    - position_space: The space in which positions are computed (RESIDUE or CODON)
    - vector_keys: Mapping of feature keys to their aggregation specifications
    - compute_vector: Method to compute per-position feature values

    Example:
        >>> class GCContentFeature:
        ...     @property
        ...     def position_space(self):
        ...         return PositionSpace.RESIDUE
        ...
        ...     @property
        ...     def vector_keys(self):
        ...         return {"gc": AggregationSpec(name="MEAN", aggregation_fn=np.mean)}
        ...
        ...     def compute_vector(self, record, **kwargs):
        ...         seq = str(record.seq).upper()
        ...         gc_vector = np.array([1.0 if b in 'GC' else 0.0 for b in seq])
        ...         return {"gc": gc_vector}
    """

    @property
    def position_space(self) -> PositionSpace:
        """The position space for this feature (RESIDUE or CODON).

        Returns:
            PositionSpace indicating whether features are computed per residue or codon
        """
        ...

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Mapping of feature keys to aggregation specifications.

        Returns:
            Dictionary mapping feature names to AggregationSpec objects that define
            how per-position values should be aggregated into window values
        """
        ...

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-position feature values for the entire sequence.

        Args:
            record: The SeqRecord containing the sequence
            **kwargs: Additional parameters that may be needed, including:
                - positions: Optional numpy array of position indices to compute.
                  If provided, only compute features for these positions in the
                  position_space (e.g., codon indices or residue indices).
                  If not provided, compute for all positions in the sequence.

        Returns:
            Dictionary mapping feature names to numpy arrays of per-position values.
            Array length should match the sequence length in the position_space (if
            positions not specified) or the length of the positions array (if specified).
            When positions is specified, the returned arrays should be indexed such that
            result[i] corresponds to position positions[i].
        """
        ...
