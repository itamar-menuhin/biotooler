"""Start region accessibility summary wrapper for ViennaRNA family.

Thin wrapper that generates window starts for a region and computes aggregated
unpaired probability (PU) values.
"""

from typing import TYPE_CHECKING

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar

if TYPE_CHECKING:
    pass  # RNA module imported lazily


class AccessibilityStartRegion:
    """Thin wrapper that computes aggregated PU for windows in a start region.

    This wrapper generates window starts for a specified region (start, end)
    with a given step size, then delegates to ViennaRNAAccessibility for
    per-nucleotide PU computation and aggregates into windows. It provides
    a convenient interface for computing accessibility summaries across a region
    without hardcoding window parameters.

    The wrapper:
    - Generates window starts: [region_start, region_start+step, ..., < region_end]
    - Computes full PU vector once via ViennaRNAAccessibility
    - Aggregates PU values into windows using mean
    - Returns dict with stable keys "PU_<window_start_nt>"
    - Avoids recomputation by computing PU vector once

    Args:
        region_start: Start position of the region (0-based, inclusive)
        region_end: End position of the region (0-based, exclusive)
        window_size: Size of each window in nucleotides
        step: Step size between window starts in nucleotides

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = AccessibilityStartRegion(
        ...     region_start=0, region_end=30, window_size=15, step=10
        ... )
        >>> record = SeqRecord(Seq("ACGUACGUACGUACGUACGUACGUACGUACGU"), id="test")
        >>> result = feature(record)
        >>> "PU_0" in result
        True
        >>> "PU_10" in result
        True
        >>> "PU_20" in result
        True
    """

    def __init__(
        self,
        region_start: int,
        region_end: int,
        window_size: int,
        step: int,
    ):
        """Initialize AccessibilityStartRegion.

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
        """Compute aggregated PU for windows in the start region.

        Args:
            record: DNA or RNA SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to aggregated PU values:
            - PU_<window_start>: Mean PU for window starting at position <window_start>

        Raises:
            ImportError: If ViennaRNA is not installed
        """
        from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility

        # Generate window starts for the region
        window_starts = list(range(self.region_start, self.region_end, self.step))

        # Handle empty window_starts
        if not window_starts:
            return {}

        # Compute full PU vector once (no recomputation)
        feature = ViennaRNAAccessibility()
        vectors = feature.compute_vector(record)
        pu_vector = vectors["PU"]

        # Aggregate PU values into windows
        result: dict[str, Scalar] = {}
        seq_len = len(pu_vector)

        for window_start in window_starts:
            window_end = window_start + self.window_size

            # Skip if window extends beyond sequence
            if window_start >= seq_len:
                continue

            # Clip window_end to sequence length
            window_end = min(window_end, seq_len)

            # Aggregate PU values using mean
            window_pu = pu_vector[window_start:window_end]
            if len(window_pu) > 0:
                result[f"PU_{window_start}"] = float(np.mean(window_pu))

        return result
