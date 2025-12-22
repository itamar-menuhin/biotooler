"""Feature set management and computation."""

__all__ = ["CodonBiasFeature", "FeatureSet", "IncrementalFeature"]

from biotooler.features.base import IncrementalFeature
from biotooler.features.sets import FeatureSet


def __getattr__(name):
    """Lazy import for CodonBiasFeature to avoid loading codonbias unnecessarily.

    This implements PEP 562 to defer importing CodonBiasFeature (and its codonbias
    dependency) until it's actually used.
    """
    if name == "CodonBiasFeature":
        from biotooler.features.codon_bias import CodonBiasFeature

        return CodonBiasFeature
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
