"""Disorder feature family for biotooler.

Protein disorder prediction using metapredict backend with optional future backends.
"""

__all__ = ["FAMILY_META", "get_features"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "disorder",
    "extra": "disorder",  # pip install biotooler[disorder]
    "owner": "@itamar-menuhin",  # GitHub username or team for CODEOWNERS
    "summary": "Protein disorder prediction using metapredict",
}


def get_features() -> list[type]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight. The metapredict dependency is only loaded when
    this function is called.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed.

    Raises:
        ImportError: If metapredict is not installed.
            The error message includes installation instructions.
    """
    # Import integration to verify metapredict is available
    from biotooler.families.disorder.integration import require_metapredict

    # This will raise ImportError with helpful message if metapredict is not installed
    require_metapredict()

    # Future feature classes will be imported here
    # For now, return an empty list as no features are implemented yet
    return []
