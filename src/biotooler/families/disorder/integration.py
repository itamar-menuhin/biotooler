"""Lazy import integration for disorder feature family.

This module provides helper functions for lazy loading of metapredict dependency.
"""

from biotooler.core.lazy_import import lazy_import


def require_metapredict():
    """Require metapredict dependency with lazy loading.

    This function lazily imports metapredict and raises a helpful error message
    if it's not installed.

    Returns:
        The metapredict module

    Raises:
        ImportError: If metapredict is not installed, with installation instructions
    """
    return lazy_import(
        "metapredict",
        extra="disorder",
        purpose="computing protein disorder features using metapredict",
    )
