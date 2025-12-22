"""{{FAMILY_NAME_TITLE}} feature implementation.

{{FAMILY_DESCRIPTION}}
"""

from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar

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
