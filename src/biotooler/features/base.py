"""Base interfaces for feature computation."""

from typing import Protocol

from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar


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
