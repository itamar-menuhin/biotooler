# Adding a New Feature Family

This guide explains how to add a new feature family to biotooler, with emphasis on implementing correct windowing semantics.

## Overview

A **feature family** is a collection of related features that compute bioinformatics metrics on sequences. Each family is isolated, has its own documentation, and must follow strict API contracts to ensure correctness, especially around windowing behavior.

## Quick Start

1. Use the family template in `templates/family_skeleton/`
2. Implement the required API surface (see below)
3. Document your family thoroughly following the template structure
4. Add tests validating windowing correctness
5. Register your family in `src/biotooler/families/__init__.py`

## Required API Surface

Every feature family must provide certain methods and properties to integrate with biotooler's windowing framework. The **recommended approach** for new families is to implement the **PositionalFeature protocol**.

### PositionalFeature Protocol (Recommended)

The PositionalFeature protocol enforces the new windowing semantics:
1. **Compute per-position values** on the full sequence
2. **Aggregate per-position values** into windows using specified aggregation functions

#### Required Components

**1. `position_space` property**

Defines the granularity at which per-position values are computed:

```python
from biotooler.features.aggregation import PositionSpace

@property
def position_space(self) -> PositionSpace:
    """The position space for this feature."""
    return PositionSpace.RESIDUE  # or PositionSpace.CODON
```

- `PositionSpace.RESIDUE`: Features computed per nucleotide or amino acid
- `PositionSpace.CODON`: Features computed per codon (groups of 3 nucleotides)

**2. `vector_keys` property**

Maps feature keys to their aggregation specifications:

```python
from biotooler.features.aggregation import AggregationSpec
import numpy as np

@property
def vector_keys(self) -> dict[str, AggregationSpec]:
    """Aggregation specs for each feature key."""
    return {
        "gc_mean": AggregationSpec(aggregation_fn=np.mean),
        "gc_sum": AggregationSpec(aggregation_fn=np.sum),
        "gc_max": AggregationSpec(aggregation_fn=np.max),
    }
```

Common aggregation functions:
- `np.mean`: Average value across window
- `np.sum`: Total sum across window
- `np.max` / `np.min`: Maximum/minimum value in window
- `geometric_mean`: Geometric mean (from `biotooler.features.aggregation`)

**3. `compute_vector` method**

Computes per-position feature values for the **entire sequence**:

```python
def compute_vector(
    self, record: SeqRecord, **kwargs
) -> dict[str, np.ndarray]:
    """Compute per-position feature values for the full sequence.
    
    Args:
        record: The SeqRecord containing the full sequence
        **kwargs: Additional parameters (e.g., orf coordinates)
    
    Returns:
        Dictionary mapping feature keys to numpy arrays of per-position values.
        Array length must match sequence length in the position_space.
    """
    seq = str(record.seq).upper()
    
    # Example: GC content as binary indicator per position
    gc_vector = np.array([1.0 if base in "GC" else 0.0 for base in seq])
    
    return {
        "gc_mean": gc_vector,
        "gc_sum": gc_vector,
        "gc_max": gc_vector,
    }
```

**IMPORTANT**: The `compute_vector` method must:
- Process the **full sequence** (not individual windows)
- Return arrays matching the length of the sequence in `position_space`
- Return the same keys as defined in `vector_keys`

### How Windowing Works

When a window is processed:

1. `compute_vector()` is called **once** for the full sequence
2. The windowing framework slices the per-position arrays to window boundaries
3. Each window slice is aggregated using the corresponding `aggregation_fn`

Example:
```python
# Full sequence: "ACGTGCTA"
# compute_vector returns: {"gc_mean": [0, 1, 1, 0, 1, 1, 0, 0]}

# For window [0:4]: "ACGT"
# Window slice: [0, 1, 1, 0]
# Aggregated value: np.mean([0, 1, 1, 0]) = 0.5
```

## Using Upstream Positional APIs

When wrapping external libraries, you must ensure the library is called with the **full sequence context**, not individual windows.

### Pattern: Library with Positional Output

If an upstream library naturally returns per-position values:

```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    """Wrap upstream library that returns positional data."""
    # ✅ CORRECT: Call library with full sequence
    seq_str = str(record.seq)
    upstream_result = upstream_library.analyze_full_sequence(seq_str)
    
    # Transform library output into per-position arrays
    feature_vector = np.array(upstream_result.per_position_scores)
    
    return {
        "feature_mean": feature_vector,
        "feature_sum": feature_vector,
    }
```

### Pattern: Library with Scalar Output Only

If the library only computes scalar values, you may need to compute windows yourself:

```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    """Wrap upstream library that only returns scalars."""
    seq_str = str(record.seq)
    
    # Option 1: If library has a per-position API, use it
    if hasattr(upstream_library, "per_position_analysis"):
        result = upstream_library.per_position_analysis(seq_str)
        return {"feature_mean": np.array(result)}
    
    # Option 2: Compute a sliding window over the full sequence yourself
    # (This should be a last resort; prefer libraries with positional APIs)
    window_size = 1  # Smallest unit for per-position
    values = []
    for i in range(len(seq_str)):
        window = seq_str[i:i+window_size]
        value = upstream_library.compute(window)
        values.append(value)
    
    return {"feature_mean": np.array(values)}
```

### Pattern: Using ORF Context

Some features need ORF coordinates for translation-aware analysis:

```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    """Feature that uses ORF context for codon-level analysis."""
    # Extract ORF coordinates if provided
    orf = kwargs.get("orf")
    if orf is None:
        # Default to full sequence as ORF
        orf = (0, len(record.seq))
    
    orf_start, orf_end = orf
    orf_seq = str(record.seq[orf_start:orf_end])
    
    # Call upstream library with full ORF sequence
    result = upstream_library.analyze_codons(orf_seq)
    
    # Return per-codon values (position_space = CODON)
    return {
        "codon_feature": np.array(result.per_codon_scores)
    }
```

### Anti-Patterns to Avoid

❌ **DON'T** compute features on individual windows:
```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    # ❌ WRONG: This defeats the purpose of full-context computation
    window_start = kwargs.get("window_start", 0)
    window_end = kwargs.get("window_end", len(record.seq))
    window_seq = str(record.seq[window_start:window_end])
    result = upstream_library.analyze(window_seq)
    return {"feature": np.array([result])}
```

❌ **DON'T** cache only window-specific results:
```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    # ❌ WRONG: This loses full-context information
    if not hasattr(self, "_cache"):
        window = kwargs.get("window_seq")
        self._cache = upstream_library.analyze(window)
    return {"feature": np.array([self._cache])}
```

✅ **DO** compute on full sequence and let framework handle windowing:
```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    # ✅ CORRECT: Full sequence analysis
    full_seq = str(record.seq)
    result = upstream_library.analyze_full_sequence(full_seq)
    return {"feature": np.array(result.per_position_values)}
```

## Documentation Requirements

Every family must have a comprehensive `README.md` following the template structure. Key sections:

### Required Sections

1. **What this family provides**: Concise summary of features
2. **Intuition**: Biological/computational context
3. **Mathematical formulation**: Precise definitions and algorithms
4. **Features and output schema**: Exact output format
5. **References**: Papers, algorithms, or resources (must include at least one http(s) link)
6. **Upstream library links**: External libraries used (or "N/A")
7. **Examples**: Comprehensive usage examples
8. **Edge cases and validation**: How edge cases are handled
9. **Windowing correctness**: How windowing semantics are enforced (NEW)
10. **Maintenance notes**: Dependencies, known issues, future work

### Windowing Correctness Section

This section must document:

- **Position space**: RESIDUE or CODON
- **Vector computation**: How full-context vectors are computed
- **Aggregation strategy**: Which aggregation functions are used and why
- **Testing approach**: How windowing correctness is validated

Example:
```markdown
## Windowing correctness

### Position space

This family computes features in **residue space** (per-nucleotide positions).

### Vector computation

The `compute_vector` method calls `upstream_library.analyze_sequence()` with 
the full sequence, ensuring complete context is available. Per-position scores 
are extracted from the library's `scores` attribute.

### Aggregation strategy

- `feature_mean`: Uses `np.mean` to compute average score across window
- `feature_max`: Uses `np.max` to identify peak score in window

### Testing approach

Tests validate that:
1. `compute_vector` returns arrays of length equal to sequence length
2. Window aggregation matches manual computation
3. Different window sizes produce consistent results
4. Boundary conditions are handled correctly

See `tests/families/my_family/test_windowing.py` for implementation.
```

If your family wraps an upstream library, you must also provide links to the upstream library's positional APIs in this section.

## Testing Requirements

Every family must include tests that validate windowing correctness:

```python
def test_compute_vector_full_context():
    """Verify that compute_vector uses full sequence context."""
    feature = MyFeature()
    record = SeqRecord(Seq("ACGTACGTACGT"), id="test")
    
    result = feature.compute_vector(record)
    
    # Vector length must match full sequence
    assert len(result["feature_key"]) == len(record.seq)

def test_window_aggregation():
    """Verify window aggregation matches expected values."""
    feature = MyFeature()
    record = SeqRecord(Seq("AAAACCCCGGGGTTTT"), id="test")
    
    # Compute full-context vector
    vectors = feature.compute_vector(record)
    
    # Manually aggregate a window
    window_start, window_end = 4, 8
    expected_mean = np.mean(vectors["feature_mean"][window_start:window_end])
    
    # Compare with framework aggregation (requires FeatureSet integration)
    # ... framework comparison code ...
    
    assert abs(actual_mean - expected_mean) < 1e-10

def test_window_boundaries():
    """Test windows at sequence boundaries."""
    feature = MyFeature()
    # Test empty window, single-base window, full-sequence window, etc.
    # ...
```

## Legacy Patterns

### IncrementalFeature Protocol (Legacy)

Older families may use the `IncrementalFeature` protocol with `init_state`, `step_state`, and `emit` methods. This pattern is still supported but **not recommended** for new families.

If you must use incremental computation:
- Document why the PositionalFeature approach is insufficient
- Ensure incremental results exactly match baseline computation
- Add tests comparing incremental vs. baseline modes

## Registration

After implementing your family, register it in `src/biotooler/families/__init__.py`:

```python
# In FAMILY_META or similar registry structure
FAMILIES = {
    "my_family": {
        "name": "my_family",
        "owner": "github_username",
        "summary": "Brief description",
    },
    # ... other families
}
```

## Checklist for New Families

- [ ] Implement `position_space` property
- [ ] Implement `vector_keys` property
- [ ] Implement `compute_vector` method using full sequence
- [ ] Document all required README sections
- [ ] Include "Windowing correctness" section in README
- [ ] Add tests validating windowing semantics
- [ ] Document upstream library links (or "N/A")
- [ ] Run `scripts/check_family_docs.py` and ensure it passes
- [ ] Register family in `src/biotooler/families/__init__.py`

## Resources

- [Windowing Semantics Documentation](./windowing.md): Technical details on windowing
- [Family Template](../../templates/family_skeleton/): Starter template
- [Example: basic_stats family](../../src/biotooler/families/basic_stats/): Reference implementation

## Getting Help

If you have questions about implementing a family:
1. Review existing families in `src/biotooler/families/`
2. Check the windowing documentation in `docs/dev/windowing.md`
3. Open a discussion issue on GitHub
4. Contact the biotooler maintainers
