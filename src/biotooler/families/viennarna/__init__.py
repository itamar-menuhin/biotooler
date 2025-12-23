"""ViennaRNA feature family for biotooler.

RNA secondary structure prediction and analysis using ViennaRNA.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    pass  # No features implemented yet

__all__ = ["FAMILY_META", "get_features"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "viennarna",
    "extra": "viennarna",  # pip install biotooler[viennarna]
    "owner": "@itamar-menuhin",  # GitHub username or team for CODEOWNERS
    "summary": "RNA secondary structure prediction and analysis using ViennaRNA",
}


def get_features() -> list[type]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight. The ViennaRNA dependency is only loaded when
    this function is called.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed.
        Currently returns an empty list as no features are implemented yet.

    Raises:
        ImportError: If ViennaRNA is not installed.
            The error message includes installation instructions.
    """
    # Import integration to verify ViennaRNA is available
    from biotooler.families.viennarna.integration import require_viennarna

    # This will raise ImportError with helpful message if RNA is not installed
    require_viennarna()

    # No features implemented yet, but the above ensures proper error handling
    return []
