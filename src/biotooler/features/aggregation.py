"""Aggregation specifications for positional features."""

from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

import numpy as np


class PositionSpace(Enum):
    """Position space for per-position feature computation.

    Attributes:
        RESIDUE: Features computed per residue position (nucleotide or amino acid)
        CODON: Features computed per codon position
    """

    RESIDUE = "residue"
    CODON = "codon"


@dataclass(frozen=True)
class AggregationSpec:
    """Specification for how to aggregate per-position values into windows.

    Attributes:
        name: Human-readable name for the aggregation (e.g., "MEAN", "GEOMEAN", "MAX")
        aggregation_fn: Function that aggregates an array of values into a scalar.
                       Common examples: np.mean, geometric_mean, np.sum, np.max
    """

    name: str
    aggregation_fn: Callable[[np.ndarray], float]


def geometric_mean(values: np.ndarray) -> float:
    """Compute geometric mean of values.

    Args:
        values: Array of values to aggregate

    Returns:
        Geometric mean of values, or 0.0 if any value is non-positive
    """
    if len(values) == 0:
        return 0.0
    if np.any(values <= 0):
        return 0.0
    return float(np.exp(np.mean(np.log(values))))
