"""Window MFE feature for ViennaRNA family.

Computes minimum free energy (MFE) for requested window substrings.
"""

from typing import TYPE_CHECKING

from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_str
from biotooler.core.types import Scalar

if TYPE_CHECKING:
    pass  # RNA module imported lazily


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

    def __init__(self, window_starts: list[int], window_size: int):
        """Initialize WindowMFEFeature.

        Args:
            window_starts: List of start positions for windows to compute
            window_size: Size of each window in nucleotides

        Raises:
            ValueError: If window_starts is empty or window_size is not positive
        """
        if not window_starts:
            raise ValueError("window_starts must not be empty")
        if window_size <= 0:
            raise ValueError("window_size must be positive")

        self.window_starts = list(window_starts)
        self.window_size = window_size

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
        # Import RNA module lazily
        from biotooler.families.viennarna.integration import require_viennarna

        RNA = require_viennarna()

        # Get sequence string
        seq_str = get_seq_str(record)

        result: dict[str, Scalar] = {}

        # Compute MFE for each requested window
        for window_start in self.window_starts:
            window_end = window_start + self.window_size

            # Skip if window extends beyond sequence
            if window_end > len(seq_str):
                continue

            # Extract window substring
            window_seq = seq_str[window_start:window_end]

            # Normalize DNA to RNA (T->U)
            rna_seq = window_seq.replace("T", "U")

            # Compute MFE using ViennaRNA
            fc = RNA.fold_compound(rna_seq)
            _, mfe = fc.mfe()

            # Store with wide key
            result[f"MFE_{window_start}"] = mfe

        return result
