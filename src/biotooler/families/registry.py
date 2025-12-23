"""Static registry of feature families without importing family modules.

This module provides lightweight discovery of available feature families and their
metadata. It does NOT import any family modules or their dependencies, keeping the
import fast and deterministic.

The registry is manually maintained and should be updated when new families are added.
"""

from typing import Any

# Static registry of families with their metadata
# This is manually maintained to avoid importing any family modules
FAMILIES: dict[str, dict[str, Any]] = {
    "chimera": {
        "summary": "Gene expression prediction and optimization using Chimera algorithms",
        "extra": "chimera",  # pip install extra name
        "owner": "@itamar-menuhin",
        "heavy": True,  # Has heavy optional dependencies
    },
    "codon_bias": {
        "summary": "Codon usage bias features using external codonbias package",
        "extra": "codon_bias",  # pip install extra name
        "owner": "@itamar-menuhin",
        "heavy": True,  # Has heavy optional dependencies
    },
    "disorder": {
        "summary": "Protein disorder prediction using metapredict",
        "extra": "disorder",  # pip install extra name
        "owner": "@itamar-menuhin",
        "heavy": True,  # Has heavy optional dependencies
    },
    "protparam": {
        "summary": "Protein physicochemical parameters using Bio.SeqUtils.ProtParam",
        "extra": None,  # No extra dependencies needed
        "owner": "@itamar-menuhin",
        "heavy": False,  # No heavy optional dependencies
    },
    "viennarna": {
        "summary": "RNA secondary structure prediction and analysis using ViennaRNA",
        "extra": "viennarna",  # pip install extra name
        "owner": "@itamar-menuhin",
        "heavy": True,  # Has heavy optional dependencies
    },
}


def list_families() -> list[dict[str, Any]]:
    """List all available feature families with their metadata.

    Returns a sorted list of family dictionaries, each containing:
    - name: Family name (str)
    - summary: Brief description (str)
    - extra: Optional dependency extra name or None (str | None)
    - owner: GitHub username of the family owner (str)
    - heavy: Whether family has heavy optional dependencies (bool)

    Returns:
        List of family metadata dictionaries, sorted by name

    Examples:
        >>> families = list_families()
        >>> len(families) >= 1
        True
        >>> families[0]['name']
        'codon_bias'
    """
    result = []
    for name, meta in sorted(FAMILIES.items()):
        family_dict = {"name": name, **meta}
        result.append(family_dict)
    return result


def get_family_meta(name: str) -> dict[str, Any]:
    """Get metadata for a specific family.

    Args:
        name: Family name (e.g., "codon_bias")

    Returns:
        Dictionary with family metadata including:
        - name: Family name
        - summary: Brief description
        - extra: Optional dependency extra name or None
        - owner: GitHub username of the owner
        - heavy: Whether family has heavy optional dependencies

    Raises:
        ValueError: If family name is not found in registry

    Examples:
        >>> meta = get_family_meta("codon_bias")
        >>> meta['name']
        'codon_bias'
        >>> get_family_meta("nonexistent")
        Traceback (most recent call last):
        ...
        ValueError: Unknown family 'nonexistent'. Available families: codon_bias
    """
    if name not in FAMILIES:
        available = ", ".join(sorted(FAMILIES.keys()))
        raise ValueError(f"Unknown family '{name}'. Available families: {available}")

    return {"name": name, **FAMILIES[name]}


def get_install_extras() -> dict[str, list[str]]:
    """Get mapping of install extras to families that require them.

    Returns a dictionary mapping pip install extra names to lists of family names
    that use those extras. Families without extras are not included.

    Returns:
        Dictionary mapping extra name to list of family names

    Examples:
        >>> extras = get_install_extras()
        >>> 'codon_bias' in extras
        True
        >>> 'codon_bias' in extras['codon_bias']
        True
    """
    result: dict[str, list[str]] = {}
    for name, meta in FAMILIES.items():
        extra = meta.get("extra")
        if extra is not None:
            if extra not in result:
                result[extra] = []
            result[extra].append(name)
    return result
