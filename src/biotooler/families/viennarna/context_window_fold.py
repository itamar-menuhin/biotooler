"""Context Window Fold feature for ViennaRNA family.

Computes MFE for a window in the context of flanking regions and reports
window-local pairing behavior in that context.
"""

from typing import TYPE_CHECKING

from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_str
from biotooler.core.types import Scalar

if TYPE_CHECKING:
    from biotooler.families.viennarna.cache import ViennaFoldCache


class ContextWindowFoldFeature:
    """Feature that computes MFE for windows with flanking context.

    This feature computes the minimum free energy (MFE) for RNA secondary
    structure prediction using ViennaRNA, but folds a larger context region
    (window + flanks) and reports both the context MFE and window-local metrics.

    Unlike WindowMFEFeature which folds only the window substring, this feature:
    - Folds window + left flank + right flank as a single context slice
    - Reports the context MFE (CTX_MFE_<start>)
    - Computes fraction of window nucleotides paired to positions outside the window
      (PAIR_OUT_FRAC_<start>)

    The feature:
    - Normalizes DNA input to RNA (T->U) before folding
    - Computes MFE via RNA.fold_compound(ctx_seq).mfe()
    - Returns dict with wide keys for each requested window start
    - Only computes the requested window starts (no full scan)

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = ContextWindowFoldFeature(
        ...     starts_nt=[10],
        ...     window_size_nt=20,
        ...     flank_left_nt=5,
        ...     flank_right_nt=5
        ... )
        >>> record = SeqRecord(Seq("ACGT" * 20), id="test")
        >>> result = feature(record)
        >>> "CTX_MFE_10" in result
        True
        >>> "PAIR_OUT_FRAC_10" in result
        True
    """

    def __init__(
        self,
        starts_nt: list[int],
        window_size_nt: int,
        flank_left_nt: int = 0,
        flank_right_nt: int = 0,
        mode: str = "mfe",
        cache: "ViennaFoldCache | None" = None,
    ):
        """Initialize ContextWindowFoldFeature.

        Args:
            starts_nt: List of start positions for windows to compute
            window_size_nt: Size of each window in nucleotides
            flank_left_nt: Number of nucleotides to include as left flank (default: 0)
            flank_right_nt: Number of nucleotides to include as right flank (default: 0)
            mode: Folding mode - currently only "mfe" is supported (default: "mfe")
            cache: Optional ViennaFoldCache for reusing fold results (default: None)

        Raises:
            ValueError: If starts_nt is empty, window_size_nt is not positive,
                        flank sizes are negative, or mode is not "mfe"
        """
        if not starts_nt:
            raise ValueError("starts_nt must not be empty")
        if window_size_nt <= 0:
            raise ValueError("window_size_nt must be positive")
        if flank_left_nt < 0:
            raise ValueError("flank_left_nt must be non-negative")
        if flank_right_nt < 0:
            raise ValueError("flank_right_nt must be non-negative")
        if mode != "mfe":
            raise ValueError(f"mode must be 'mfe', got {mode!r}")

        self.starts_nt = list(starts_nt)
        self.window_size_nt = window_size_nt
        self.flank_left_nt = flank_left_nt
        self.flank_right_nt = flank_right_nt
        self.mode = mode
        self.cache = cache

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute context fold metrics for requested window positions.

        Args:
            record: DNA or RNA SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to values:
            - CTX_MFE_<start>: MFE for context slice (window + flanks)
            - PAIR_OUT_FRAC_<start>: Fraction of window nts paired to outside window

        Raises:
            ImportError: If ViennaRNA is not installed
        """
        from biotooler.families.viennarna.cache import get_context_mfe_cached

        # Get sequence string
        seq_str = get_seq_str(record)
        seq_len = len(seq_str)

        result: dict[str, Scalar] = {}

        # Compute metrics for each requested window
        for start in self.starts_nt:
            window_end = start + self.window_size_nt

            # Skip if window extends beyond sequence
            if window_end > seq_len:
                continue

            # Use cached fold computation with flanks
            structure, mfe = get_context_mfe_cached(
                seq_str=seq_str,
                window_start=start,
                window_size=self.window_size_nt,
                flank_left=self.flank_left_nt,
                flank_right=self.flank_right_nt,
                cache=self.cache,
            )

            # Store context MFE
            result[f"CTX_MFE_{start}"] = mfe

            # Compute PAIR_OUT_FRAC: fraction of window nts paired to outside window
            # Window positions in context slice coordinates
            ctx_start = max(0, start - self.flank_left_nt)
            window_ctx_start = start - ctx_start
            window_ctx_end = window_ctx_start + self.window_size_nt

            pair_out_frac = self._compute_pair_out_frac(
                structure, window_ctx_start, window_ctx_end
            )
            result[f"PAIR_OUT_FRAC_{start}"] = pair_out_frac

        return result

    def _compute_pair_out_frac(
        self, structure: str, window_start: int, window_end: int
    ) -> float:
        """Compute fraction of window nucleotides paired to positions outside window.

        Args:
            structure: Dot-bracket structure string
            window_start: Start position of window in structure (0-based)
            window_end: End position of window in structure (exclusive)

        Returns:
            Fraction in [0, 1] of window nts paired to outside window
        """
        # Parse dot-bracket structure to find base pairs
        pairs = self._parse_dotbracket(structure)

        # Count window nucleotides paired to outside window
        count_paired_out = 0
        window_nts = 0

        for i in range(window_start, window_end):
            window_nts += 1
            if i in pairs:
                partner = pairs[i]
                # Check if partner is outside window
                if partner < window_start or partner >= window_end:
                    count_paired_out += 1

        if window_nts == 0:
            return 0.0

        return count_paired_out / window_nts

    def _parse_dotbracket(self, structure: str) -> dict[int, int]:
        """Parse dot-bracket structure to find base pairs.

        Args:
            structure: Dot-bracket structure string

        Returns:
            Dictionary mapping position to its paired position (0-based indexing)
        """
        pairs: dict[int, int] = {}
        stack: list[int] = []

        for i, char in enumerate(structure):
            if char == "(":
                stack.append(i)
            elif char == ")":
                if stack:
                    j = stack.pop()
                    pairs[j] = i
                    pairs[i] = j
            # '.' means unpaired, no action needed

        return pairs
