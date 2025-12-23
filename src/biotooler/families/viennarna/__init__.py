"""ViennaRNA feature family for biotooler.

RNA secondary structure prediction and analysis using ViennaRNA.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotooler.families.viennarna.window_mfe import WindowMFEFeature

__all__ = ["FAMILY_META", "get_features", "WindowMFEFeature"]

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

    Raises:
        ImportError:  If ViennaRNA is not installed.
            The error message includes installation instructions.
    """
    # Import integration to verify ViennaRNA is available
    from biotooler.families.viennarna.integration import require_viennarna

    # This will raise ImportError with helpful message if RNA is not installed
    require_viennarna()

    # Import feature classes
    from biotooler.families.viennarna.accessibility import ViennaRNAAccessibility
    from biotooler.families.viennarna.window_mfe import WindowMFEFeature

    return [ViennaRNAAccessibility, WindowMFEFeature]


def __getattr__(name: str):
    """Lazy import for WindowMFEFeature to maintain backward compatibility.

    This allows code like `from biotooler.families.viennarna import WindowMFEFeature`
    to work while keeping imports lightweight and not loading ViennaRNA until needed.
    """
    if name == "WindowMFEFeature":
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        return WindowMFEFeature
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
