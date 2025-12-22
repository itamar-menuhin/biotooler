"""Lazy import integration for {{FAMILY_NAME}} feature family.

This module provides documentation for integrating the {{FAMILY_NAME_TITLE}}Feature
into biotooler.features with lazy loading of optional dependencies.

## Using lazy_import for Optional Dependencies

If your feature has optional dependencies (e.g., external packages not in core biotooler),
use the lazy_import helper in your feature implementation:

```python
from biotooler.core.lazy_import import lazy_import

# Lazy import optional dependency
optional_package = lazy_import(
    "package_name",
    extra="{{FAMILY_NAME}}",
    purpose="computing {{FAMILY_NAME}} features"
)
```

Then update pyproject.toml to add the optional dependency:

```toml
[project.optional-dependencies]
{{FAMILY_NAME}} = [
    "package-name",
]
```

## Adding to biotooler.features Module (Optional)

To make this feature accessible via `from biotooler.features import {{FAMILY_NAME_TITLE}}Feature`,
add the following to src/biotooler/features/__init__.py:

1. Add "{{FAMILY_NAME_TITLE}}Feature" to __all__

2. Add TYPE_CHECKING import:
   ```python
   if TYPE_CHECKING:
       from biotooler.families.{{FAMILY_NAME}} import {{FAMILY_NAME_TITLE}}Feature
   ```

3. Add case to __getattr__:
   ```python
   if name == "{{FAMILY_NAME_TITLE}}Feature":
       from biotooler.families.{{FAMILY_NAME}} import {{FAMILY_NAME_TITLE}}Feature
       return {{FAMILY_NAME_TITLE}}Feature
   ```

This enables lazy loading where the feature family is only imported when explicitly accessed.
"""

# This file serves as documentation for the integration.
# No code is required here unless there are specific integration utilities needed.
