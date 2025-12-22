"""Feature set management and computation."""

from typing import TYPE_CHECKING

__all__ = ["CodonBiasFeature", "FeatureSet", "IncrementalFeature"]

from biotooler.features.base import IncrementalFeature
from biotooler.features.sets import FeatureSet

# Type checking import - doesn't execute at runtime
if TYPE_CHECKING:
    from biotooler.families.codon_bias import CodonBiasFeature


def __dir__():
    """Support for dir() to include lazily imported names."""
    return __all__


def __getattr__(name):
    """Lazy import for CodonBiasFeature to avoid loading codonbias unnecessarily.

    This implements PEP 562 to defer importing CodonBiasFeature (and its codonbias
    dependency) until it's actually used.
    """
    if name == "CodonBiasFeature":
        from biotooler.families.codon_bias import CodonBiasFeature

        return CodonBiasFeature
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
