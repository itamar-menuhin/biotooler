"""Lazy import integration for disorder feature family.

This module provides helper functions for lazy loading of metapredict and idrpred dependencies.
"""

import shutil

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


def require_idrpred_cli():
    """Require idrpred CLI tool availability.

    This function checks if the idrpred command-line tool is available on PATH
    and raises a helpful error message if it's not found.

    Raises:
        ImportError: If idrpred CLI is not available, with installation instructions
    """
    if shutil.which("idrpred") is None:
        msg = (
            "The 'idrpred' command-line tool is required for computing protein disorder "
            "features using IDRPred. "
            'Install it with: pip install "biotooler[disorder-idrpred]" '
            "or install IDRPred manually from https://github.com/matthiasblum/idrpred "
            "and ensure 'idrpred' is on PATH."
        )
        raise ImportError(msg)


def require_idrpred_pkg():
    """Require idrpred package with lazy loading.

    This function lazily imports idrpred and raises a helpful error message
    if it's not installed.

    Returns:
        The idrpred module

    Raises:
        ImportError: If idrpred is not installed, with installation instructions
    """
    return lazy_import(
        "idrpred",
        extra="disorder-idrpred",
        purpose="computing protein disorder features using IDRPred",
    )
