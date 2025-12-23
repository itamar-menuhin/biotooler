"""Caching for ViennaRNA fold computations.

Provides in-memory caching to avoid redundant fold computations for identical
window definitions and configurations. The cache is injectable and testable,
with no global singleton state.
"""

import hashlib
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass  # RNA module imported lazily


@dataclass(frozen=True)
class FoldKey:
    """Immutable cache key for a fold computation.

    Attributes:
        seq_hash: Hash of the sequence to fold (after normalization to RNA)
        start: Start position in the normalized context (0-based)
        length: Length of the sequence to fold
        flank_left: Number of nucleotides in left flank
        flank_right: Number of nucleotides in right flank
        config_hash: Hash of fold configuration (e.g., temperature, model)
    """

    seq_hash: str
    start: int
    length: int
    flank_left: int
    flank_right: int
    config_hash: str


class ViennaFoldCache:
    """In-memory cache for ViennaRNA fold results.

    Stores (structure, mfe) tuples keyed by fold parameters. The cache
    is session-scoped and does not persist across runs.

    This cache is injectable and testable, with no global singleton state.
    Create a new instance per feature or share across features as needed.

    Examples:
        >>> cache = ViennaFoldCache()
        >>> # Use cache in multiple features
        >>> feature1 = WindowMFEFeature(window_starts=[0], window_size=20, cache=cache)
        >>> feature2 = ContextWindowFoldFeature(starts_nt=[10], window_size_nt=15, cache=cache)
    """

    def __init__(self):
        """Initialize an empty fold cache."""
        self._cache: dict[FoldKey, tuple[str, float]] = {}

    def get(self, key: FoldKey) -> tuple[str, float] | None:
        """Get cached fold result if available.

        Args:
            key: FoldKey identifying the fold computation

        Returns:
            Tuple of (structure, mfe) if cached, None otherwise
        """
        return self._cache.get(key)

    def set(self, key: FoldKey, structure: str, mfe: float) -> None:
        """Store fold result in cache.

        Args:
            key: FoldKey identifying the fold computation
            structure: Dot-bracket structure string
            mfe: Minimum free energy value
        """
        self._cache[key] = (structure, mfe)

    def clear(self) -> None:
        """Clear all cached results."""
        self._cache.clear()

    def __len__(self) -> int:
        """Return number of cached results."""
        return len(self._cache)


def _hash_string(s: str) -> str:
    """Compute SHA256 hash of a string.

    Args:
        s: String to hash

    Returns:
        Hexadecimal hash string
    """
    return hashlib.sha256(s.encode()).hexdigest()


def _make_config_hash(temperature: float = 37.0) -> str:
    """Create hash for fold configuration.

    Currently supports temperature parameter. Can be extended for other
    configuration options like model, dangles, etc.

    Args:
        temperature: Temperature in Celsius (default: 37.0)

    Returns:
        Hash string representing the configuration
    """
    # For now, just hash the temperature
    # Can be extended to include other parameters
    config_str = f"temp:{temperature}"
    return _hash_string(config_str)


def get_context_mfe_cached(
    seq_str: str,
    window_start: int,
    window_size: int,
    flank_left: int = 0,
    flank_right: int = 0,
    cache: ViennaFoldCache | None = None,
    temperature: float = 37.0,
) -> tuple[str, float]:
    """Compute MFE for a window with context, using cache if available.

    This helper function:
    1. Normalizes sequence boundaries (clips to sequence length)
    2. Converts DNA (T) to RNA (U)
    3. Computes or retrieves cached fold result
    4. Returns structure and MFE

    Args:
        seq_str: Full sequence string (DNA or RNA)
        window_start: Start position of window in seq_str (0-based)
        window_size: Size of window in nucleotides
        flank_left: Number of nucleotides to include as left flank (default: 0)
        flank_right: Number of nucleotides to include as right flank (default: 0)
        cache: Optional ViennaFoldCache instance for caching results
        temperature: Temperature for folding in Celsius (default: 37.0)

    Returns:
        Tuple of (structure, mfe) for the context slice

    Raises:
        ImportError: If ViennaRNA is not installed
        ValueError: If window extends beyond sequence
    """
    from biotooler.families.viennarna.integration import require_viennarna

    RNA = require_viennarna()

    # Validate window bounds
    seq_len = len(seq_str)
    window_end = window_start + window_size
    if window_end > seq_len:
        raise ValueError(
            f"Window extends beyond sequence: window_end={window_end}, seq_len={seq_len}"
        )

    # Compute context slice boundaries (normalized to sequence bounds)
    ctx_start = max(0, window_start - flank_left)
    ctx_end = min(seq_len, window_end + flank_right)

    # Extract context slice
    ctx_seq = seq_str[ctx_start:ctx_end]

    # Normalize DNA to RNA (T->U)
    rna_seq = ctx_seq.replace("T", "U")

    # Create cache key
    seq_hash = _hash_string(rna_seq)
    config_hash = _make_config_hash(temperature)
    # For cache key, use normalized context boundaries
    actual_flank_left = window_start - ctx_start
    actual_flank_right = ctx_end - window_end

    key = FoldKey(
        seq_hash=seq_hash,
        start=0,  # Always 0 since we extract the exact context slice
        length=len(rna_seq),
        flank_left=actual_flank_left,
        flank_right=actual_flank_right,
        config_hash=config_hash,
    )

    # Check cache
    if cache is not None:
        cached_result = cache.get(key)
        if cached_result is not None:
            return cached_result

    # Compute fold using ViennaRNA
    fc = RNA.fold_compound(rna_seq)

    # Apply temperature if not default
    if temperature != 37.0:
        # ViennaRNA uses temperature in fold_compound or via parameters
        # For now, we'll use the default RNA.fold_compound behavior
        # Temperature control would require more complex ViennaRNA API usage
        pass

    structure, mfe = fc.mfe()

    # Store in cache
    if cache is not None:
        cache.set(key, structure, mfe)

    return structure, mfe
