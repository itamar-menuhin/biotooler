# Windowing Semantics

This document describes the new windowing semantics for positional features in biotooler.

## Overview

We are introducing a new approach to computing features over sequence windows:

1. **Compute per-position values** on the full sequence
2. **Aggregate per-position values** into windows using specified aggregation functions

This is an **additive change** - existing scalar features continue to work as before.

## Core Concepts

### PositionSpace

The `PositionSpace` enum defines the granularity at which per-position values are computed:

- `PositionSpace.RESIDUE`: Features are computed per residue position (nucleotide or amino acid)
- `PositionSpace.CODON`: Features are computed per codon position (groups of 3 nucleotides)

### AggregationSpec

The `AggregationSpec` dataclass specifies how per-position values should be aggregated into window-level values:

```python
from biotooler.features.aggregation import AggregationSpec
import numpy as np

# Use mean aggregation
spec = AggregationSpec(aggregation_fn=np.mean)

# Use geometric mean
from biotooler.features.aggregation import geometric_mean
spec = AggregationSpec(aggregation_fn=geometric_mean)

# Use sum
spec = AggregationSpec(aggregation_fn=np.sum)
```

### PositionalFeature Protocol

Features implementing the new `PositionalFeature` protocol must provide:

1. **`position_space`** property: Returns a `PositionSpace` enum value
2. **`vector_keys`** property: Returns a dict mapping feature keys to `AggregationSpec` objects
3. **`compute_vector(...)`** method: Returns a dict mapping feature keys to numpy arrays of per-position values

## Example

Here's a simple example of a positional feature that computes GC content:

```python
import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.features.aggregation import AggregationSpec, PositionSpace


class GCContentFeature:
    """Compute GC content as a positional feature."""

    @property
    def position_space(self) -> PositionSpace:
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        return {
            "gc_mean": AggregationSpec(aggregation_fn=np.mean),
            "gc_sum": AggregationSpec(aggregation_fn=np.sum),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-residue GC indicator."""
        seq = str(record.seq).upper()
        gc_vector = np.array([1.0 if base in "GC" else 0.0 for base in seq])
        return {
            "gc_mean": gc_vector,
            "gc_sum": gc_vector,
        }
```

## Usage in Windowing

When a window is processed:

1. The feature's `compute_vector()` method is called once for the full sequence
2. The resulting per-position arrays are sliced to the window boundaries
3. Each window slice is aggregated using the corresponding `AggregationSpec.aggregation_fn`

For example, given a sequence `"ACGTGCTA"` and a window from position 0 to 4:

1. `compute_vector()` returns `{"gc_mean": [0, 1, 1, 0, 1, 1, 0, 0]}`
2. Window slice: `[0, 1, 1, 0]`
3. Aggregation: `np.mean([0, 1, 1, 0])` = `0.5`

## Compatibility

**Important:** Existing scalar features that compute values directly on windows remain unchanged and fully supported. The new `PositionalFeature` protocol is an optional, additive enhancement for cases where per-position computation is beneficial.

## Benefits

This new approach provides several advantages:

1. **Separation of concerns**: Feature computation logic is separate from windowing logic
2. **Reusability**: Per-position values can be cached and reused across multiple windows
3. **Flexibility**: Different aggregation functions can be applied to the same per-position values
4. **Performance**: For dense window coverage, computing once and aggregating many times is more efficient

## Special Cases: Features Requiring Positional APIs

Some features have mathematical definitions that make windowing via slicing incorrect. These features **must** use the positional API to support windowing.

### Chimera Features (cARS/PScARS)

The Chimera Average Repetitive Substring (cARS) score is defined as an average over maximal common subsequences at every position in the target sequence. **Slicing the target sequence before computing cARS would be incorrect** because:

1. The maximal common subsequences are computed relative to the full target sequence
2. Slicing changes which subsequences are maximal at each position
3. This would produce mathematically incorrect results that don't reflect the true cARS definition

**Correct approach for Chimera windowing:**

The `ChimeraFeature` class implements the `PositionalFeature` protocol and uses pyChimera's `return_vec=True` API to compute per-position maximal common substring lengths across the full sequence. These per-position values are then aggregated into windows using mean aggregation (matching the scalar cARS definition).

```python
from biotooler.families.chimera import ChimeraFeature

# ChimeraFeature implements PositionalFeature protocol
feature = ChimeraFeature(reference_seqs=ref_genes, algorithm="cARS")

# position_space: CODON (operates on codon-level sequences)
assert feature.position_space == PositionSpace.CODON

# vector_keys: Returns aggregation spec using mean
assert "cARS_score" in feature.vector_keys
assert feature.vector_keys["cARS_score"].aggregation_fn == np.mean

# compute_vector: Computes per-codon maximal common substring lengths
# (uses pyChimera's calc_cARS with return_vec=True internally)
cars_vectors = feature.compute_vector(record)
# Returns: {"cARS_score": np.array([len1, len2, ..., lenN])}
```

**Why this matters:**

Without the positional API, windowing would require either:
1. Silently producing incorrect results by slicing before computation
2. Raising an error to prevent incorrect usage

By implementing the `PositionalFeature` protocol, Chimera features can support windowing correctly by construction, with no possibility of silent mathematical errors.
