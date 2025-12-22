"""Lazy import helper for optional dependencies."""

import importlib
from types import ModuleType


def lazy_import(module: str, *, extra: str, purpose: str) -> ModuleType:
    """Lazily import a module with helpful error message if missing.

    Args:
        module: The module name to import (e.g., "codonbias")
        extra: The extra name for installation (e.g., "codon_bias")
        purpose: Description of what the module is used for

    Returns:
        The imported module

    Raises:
        ImportError: If the module cannot be imported, with installation instructions

    Examples:
        >>> # Import codonbias with helpful error if missing
        >>> codonbias = lazy_import(
        ...     "codonbias",
        ...     extra="codon_bias",
        ...     purpose="computing codon usage bias features"
        ... )
    """
    try:
        return importlib.import_module(module)
    except ImportError as e:
        msg = (
            f"Module '{module}' is required for {purpose}. "
            f'Install it with: pip install "biotooler[{extra}]"'
        )
        raise ImportError(msg) from e
