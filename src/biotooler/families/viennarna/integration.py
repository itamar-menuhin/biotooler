"""Lazy import integration for viennarna feature family.

This module provides helper functions for lazy loading of ViennaRNA dependency.
"""

from biotooler.core.lazy_import import lazy_import


def require_viennarna():
    """Require ViennaRNA dependency with lazy loading.

    This function lazily imports ViennaRNA (imported as RNA) and raises a helpful
    error message if it's not installed.

    Returns:
        The RNA module from ViennaRNA package

    Raises:
        ImportError: If ViennaRNA is not installed, with installation instructions
    """
    return lazy_import(
        "RNA",
        extra="viennarna",
        purpose="computing RNA secondary structure features using ViennaRNA",
    )
