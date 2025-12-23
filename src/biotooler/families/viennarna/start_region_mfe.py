"""Start region MFE summary wrapper for ViennaRNA family.

Thin wrapper that generates window starts for a region and computes MFE values.
"""

from typing import TYPE_CHECKING

from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar

if TYPE_CHECKING:
    pass  # RNA module imported lazily


class WindowMFEStartRegion:
    """Thin wrapper that computes MFE for windows in a start region.

    This wrapper generates window starts for a specified region (start, end)
    with a given step size, then delegates to WindowMFEFeature for computation.
    It provides a convenient interface for computing MFE values across a region
    without hardcoding window parameters.

    The wrapper:
    - Generates window starts: [region_start, region_start+step, ..., < region_end]
    - Delegates to WindowMFEFeature for MFE computation
    - Returns dict with stable keys "MFE_<window_start_nt>"
    - Avoids recomputation by calling underlying feature once

    Args:
        region_start: Start position of the region (0-based, inclusive)
        region_end: End position of the region (0-based, exclusive)
        window_size: Size of each window in nucleotides
        step: Step size between window starts in nucleotides

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = WindowMFEStartRegion(
        ...     region_start=0, region_end=30, window_size=15, step=10
        ... )
        >>> record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGT"), id="test")
        >>> result = feature(record)
        >>> "MFE_0" in result
        True
        >>> "MFE_10" in result
        True
        >>> "MFE_20" in result
        True
    """

    def __init__(
        self,
        region_start: int,
        region_end: int,
        window_size: int,
        step: int,
    ):
        """Initialize WindowMFEStartRegion.

        Args:
            region_start: Start position of the region (0-based, inclusive)
            region_end: End position of the region (0-based, exclusive)
            window_size: Size of each window in nucleotides
            step: Step size between window starts in nucleotides

        Raises:
            ValueError: If parameters are invalid
        """
        if region_start < 0:
            raise ValueError("region_start must be non-negative")
        if region_end <= region_start:
            raise ValueError("region_end must be greater than region_start")
        if window_size <= 0:
            raise ValueError("window_size must be positive")
        if step <= 0:
            raise ValueError("step must be positive")

        self.region_start = region_start
        self.region_end = region_end
        self.window_size = window_size
        self.step = step

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute MFE for windows in the start region.

        Args:
            record: DNA or RNA SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to MFE values:
            - MFE_<window_start>: MFE for window starting at position <window_start>

        Raises:
            ImportError: If ViennaRNA is not installed
        """
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        # Generate window starts for the region
        window_starts = list(range(self.region_start, self.region_end, self.step))

        # Handle empty window_starts
        if not window_starts:
            return {}

        # Delegate to WindowMFEFeature
        feature = WindowMFEFeature(window_starts=window_starts, window_size=self.window_size)
        return feature(record)
