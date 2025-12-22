"""{{FAMILY_NAME_TITLE}} feature family for biotooler.

{{FAMILY_DESCRIPTION}}
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biotooler.families.{{FAMILY_NAME}}.feature import {{FAMILY_NAME_TITLE}}Feature

__all__ = ["FAMILY_META", "get_features"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "{{FAMILY_NAME}}",
    "extra": "{{FAMILY_NAME}}",  # pip install extra name if this family has optional deps
    "owner": "{{OWNER}}",  # GitHub username or team for CODEOWNERS
    "summary": "{{FAMILY_DESCRIPTION}}",  # Short description of what this family provides
}


def get_features() -> list[type["{{FAMILY_NAME_TITLE}}Feature"]]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight. Heavy optional dependencies are only loaded when
    this function is called.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed with their own parameters.

    Raises:
        ImportError: If optional dependencies are not installed.
            The error message includes installation instructions.
    """
    from biotooler.families.{{FAMILY_NAME}}.feature import {{FAMILY_NAME_TITLE}}Feature

    # Return list of feature classes (users will instantiate with their own params)
    # Customize this based on your family's needs
    return [{{FAMILY_NAME_TITLE}}Feature]
