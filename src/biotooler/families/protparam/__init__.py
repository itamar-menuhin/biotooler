"""Protein parameter feature family using Bio.SeqUtils.ProtParam.

This family provides protein physicochemical property calculations including
molecular weight, isoelectric point, GRAVY, instability index, and aromaticity.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotooler.families.protparam.feature import ProtParamFeature

__all__ = ["FAMILY_META", "get_features", "ProtParamFeature"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "protparam",
    "extra": None,  # No extra dependencies needed (Biopython is already in core deps)
    "owner": "@itamar-menuhin",
    "summary": "Protein physicochemical parameters using Bio.SeqUtils.ProtParam",
}


def get_features() -> list[type["ProtParamFeature"]]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed.
    """
    from biotooler.families.protparam.feature import ProtParamFeature

    return [ProtParamFeature]


# Maintain backward compatibility: lazy import of ProtParamFeature
def __getattr__(name: str):
    """Lazy import for ProtParamFeature to maintain backward compatibility.

    This allows code like `from biotooler.families.protparam import ProtParamFeature`
    to work while keeping imports lightweight.
    """
    if name == "ProtParamFeature":
        from biotooler.families.protparam.feature import ProtParamFeature

        return ProtParamFeature
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
