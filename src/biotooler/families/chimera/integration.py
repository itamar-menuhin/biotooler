"""Lazy import integration for chimera feature family.

This module provides helper functions for lazy loading of pyChimera dependency.
"""

from biotooler.core.lazy_import import lazy_import


def require_chimera_dep():
    """Require pyChimera dependency with lazy loading.

    This function lazily imports pyChimera and raises a helpful error message
    if it's not installed.

    Returns:
        The pyChimera module

    Raises:
        ImportError: If pyChimera is not installed, with installation instructions
    """
    return lazy_import(
        "pyChimera",
        extra="chimera",
        purpose="computing protein structure features",
    )
