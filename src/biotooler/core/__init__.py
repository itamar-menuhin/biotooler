"""Core utilities for biotooler."""

from biotooler.core.errors import InvalidSequenceError
from biotooler.core.orf_store import OrfSpan, attach_orf, get_orf, select_orf_by_index
from biotooler.core.record import coerce_record, get_molecule_type
from biotooler.core.reference_sequences import ReferenceSequenceSet
from biotooler.core.types import FeatureOutput, Scalar
from biotooler.core.windowing import iter_orf_codon_windows, iter_windows

__all__ = [
    "InvalidSequenceError",
    "Scalar",
    "FeatureOutput",
    "coerce_record",
    "get_molecule_type",
    "OrfSpan",
    "select_orf_by_index",
    "attach_orf",
    "get_orf",
    "iter_windows",
    "iter_orf_codon_windows",
    "ReferenceSequenceSet",
]
