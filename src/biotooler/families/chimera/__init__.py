"""Chimera feature family for biotooler.

Feature family for chimeric protein structure analysis using pyChimera.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotooler.families.chimera.feature import ChimeraFeature

__all__ = ["FAMILY_META", "get_features"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "chimera",
    "extra": "chimera",  # pip install biotooler[chimera]
    "owner": "@itamar-menuhin",  # GitHub username or team for CODEOWNERS
    "summary": "Chimeric protein structure analysis using pyChimera",
}


def get_features() -> list[type["ChimeraFeature"]]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight. The pyChimera dependency is only loaded when
    this function is called.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed.

    Raises:
        ImportError: If pyChimera is not installed.
            The error message includes installation instructions.
    """
    from biotooler.families.chimera.feature import ChimeraFeature

    return [ChimeraFeature]
