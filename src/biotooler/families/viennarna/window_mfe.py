"""Window MFE feature for ViennaRNA family.

Computes minimum free energy (MFE) for requested window substrings.
"""

from typing import TYPE_CHECKING

from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_str
from biotooler.core.types import Scalar

if TYPE_CHECKING:
    from biotooler.families.viennarna.cache import ViennaFoldCache


class WindowMFEFeature:
    """Feature that computes MFE for specific window substrings.

    This feature computes the minimum free energy (MFE) for RNA secondary
    structure prediction using ViennaRNA, but only for specifically requested
    window positions. Unlike per-residue features, MFE is a scalar value
    per window substring.

    The feature:
    - Normalizes DNA input to RNA (T->U) before folding
    - Computes MFE via RNA.fold_compound(seq).mfe()
    - Returns dict with wide keys "MFE_<window_start_nt>"
    - Only computes the requested window starts (no full scan)

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = WindowMFEFeature(window_starts=[0, 10], window_size=20)
        >>> record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGT"), id="test")
        >>> result = feature(record)
        >>> "MFE_0" in result
        True
        >>> "MFE_10" in result
        True
        >>> "MFE_5" in result  # Not requested
        False
    """

    def __init__(
        self,
        window_starts: list[int],
        window_size: int,
        cache: "ViennaFoldCache | None" = None,
    ):
        """Initialize WindowMFEFeature.

        Args:
            window_starts: List of start positions for windows to compute
            window_size: Size of each window in nucleotides
            cache: Optional ViennaFoldCache for reusing fold results (default: None)

        Raises:
            ValueError: If window_starts is empty or window_size is not positive
        """
        if not window_starts:
            raise ValueError("window_starts must not be empty")
        if window_size <= 0:
            raise ValueError("window_size must be positive")

        self.window_starts = list(window_starts)
        self.window_size = window_size
        self.cache = cache

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute MFE for requested window substrings.

        Args:
            record: DNA or RNA SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to MFE values:
            - MFE_<window_start>: MFE for window starting at position <window_start>

        Raises:
            ImportError: If ViennaRNA is not installed
        """
        from biotooler.families.viennarna.cache import get_context_mfe_cached

        # Get sequence string
        seq_str = get_seq_str(record)

        result: dict[str, Scalar] = {}

        # Compute MFE for each requested window
        for window_start in self.window_starts:
            # Skip if window extends beyond sequence
            if window_start + self.window_size > len(seq_str):
                continue

            # Use cached fold computation (no flanks for this feature)
            _, mfe = get_context_mfe_cached(
                seq_str=seq_str,
                window_start=window_start,
                window_size=self.window_size,
                flank_left=0,
                flank_right=0,
                cache=self.cache,
            )

            # Store with wide key
            result[f"MFE_{window_start}"] = mfe

        return result
