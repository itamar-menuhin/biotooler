"""Codon usage bias feature computation using codonbias package."""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotooler.families.codon_bias.codon_bias import CodonBiasFeature

__all__ = ["FAMILY_META", "get_features", "CodonBiasFeature"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "codon_bias",
    "extra": "codon_bias",  # pip install biotooler[codon_bias]
    "owner": "@itamar-menuhin",  # GitHub username or team for CODEOWNERS
    "summary": "Codon usage bias features using external codonbias package",
}


def get_features() -> list[type["CodonBiasFeature"]]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight. The codonbias dependency is only loaded when
    this function is called.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed.

    Raises:
        ImportError: If codonbias is not installed.
            The error message includes installation instructions.
    """
    from biotooler.families.codon_bias.codon_bias import CodonBiasFeature

    return [CodonBiasFeature]


# Maintain backward compatibility: lazy import of CodonBiasFeature
def __getattr__(name: str):
    """Lazy import for CodonBiasFeature to maintain backward compatibility.

    This allows existing code like `from biotooler.families.codon_bias import CodonBiasFeature`
    to continue working while keeping imports lightweight.
    """
    if name == "CodonBiasFeature":
        from biotooler.families.codon_bias.codon_bias import CodonBiasFeature

        return CodonBiasFeature
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
