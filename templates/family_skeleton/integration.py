"""Lazy import integration for {{FAMILY_NAME}} feature family.

This module provides lazy import support for the {{FAMILY_NAME_TITLE}}Feature
to be accessible through biotooler.features while avoiding loading heavy
dependencies until explicitly used.

To add this feature to the main features module, add the following to
src/biotooler/features/__init__.py:

1. Add "{{FAMILY_NAME_TITLE}}Feature" to __all__
2. Add TYPE_CHECKING import:
   if TYPE_CHECKING:
       from biotooler.families.{{FAMILY_NAME}} import {{FAMILY_NAME_TITLE}}Feature

3. Add case to __getattr__:
   if name == "{{FAMILY_NAME_TITLE}}Feature":
       from biotooler.families.{{FAMILY_NAME}} import {{FAMILY_NAME_TITLE}}Feature
       return {{FAMILY_NAME_TITLE}}Feature
"""

# This file serves as documentation for the integration.
# No code is required here unless there are specific integration utilities needed.
