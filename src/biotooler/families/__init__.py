"""Feature family modules for biotooler.

This module provides lightweight discovery of feature families without importing
any family-specific code. Each family lives in its own submodule for isolation
and to enable CODEOWNERS.

Available families:
    - codon_bias: Codon usage bias features using external codonbias package

To use a family, import it explicitly:
    >>> from biotooler.families.codon_bias import CodonBiasFeature

Note: This module does NOT import any family modules to keep imports lightweight.
Heavy optional dependencies (e.g., codonbias) are only loaded when explicitly imported.
"""

__all__: list[str] = []

