# {{FAMILY_NAME_TITLE}} Feature Family

## What this family provides

TODO: Provide a concise description of what this feature family computes and why it's useful for biological sequence analysis.

## Intuition

TODO: Explain the biological or computational intuition behind these features. Help users understand:
- When to use these features
- What insights they provide
- How they relate to biological phenomena

## Mathematical formulation

TODO: Provide precise mathematical definitions of the features computed. Include:
- Formal notation
- Equations
- Algorithm descriptions
- Complexity analysis where relevant

## Features and output schema

TODO: Document the exact output format and feature keys. Include:
- Feature names and their data types
- Value ranges and units
- Example outputs
- Any sequence-type specific variations

## References

TODO: Add references to papers, algorithms, or external resources that describe the methods used.

Example:
- Author A et al. (Year) Title. Journal Volume(Issue):Pages. https://doi.org/...

## Upstream library links

TODO: Add links to any upstream libraries or APIs used by this family.

If wrapping an external library, include:
- Library name and version
- GitHub/documentation URL
- Specific API references used

If no upstream library is used, write "N/A" with a brief explanation.

## Examples

TODO: Provide comprehensive usage examples covering:
- Basic usage with different sequence types (DNA, RNA, protein)
- Window-based computation
- Integration with other biotooler features
- Common use cases and workflows

## Edge cases and validation

TODO: Document how the feature handles edge cases:
- Empty sequences
- Sequences with ambiguous characters
- Very short or very long sequences
- Invalid input handling
- Boundary conditions

## Windowing correctness

TODO: Document how this family ensures correct windowing semantics.

### Position space

Specify whether this family computes features in:
- **Residue space**: Per-nucleotide or per-amino-acid positions
- **Codon space**: Per-codon positions (groups of 3 nucleotides)

### Vector computation

Explain how this family implements the `compute_vector` method:
- Does it use the full sequence context?
- How are per-position values calculated?
- What upstream library calls are made (if any)?

### Aggregation strategy

Document the aggregation functions used for each feature:
- Feature key → aggregation function mapping
- Why each aggregation function is appropriate
- Examples of aggregated values

### Testing approach

Describe how windowing correctness is validated:
- Test cases covering edge cases (window boundaries, different window sizes)
- Validation that full-context vectors are used
- Comparison with expected baseline computations

Example test approach:
```python
# Verify that compute_vector uses full sequence context
full_vector = feature.compute_vector(full_record)
assert len(full_vector["feature_key"]) == len(full_record.seq)

# Verify window aggregation matches manual computation
window_value = aggregation_fn(full_vector["feature_key"][start:end])
assert abs(window_value - expected_value) < 1e-10
```

See [docs/dev/adding_a_family.md](../../docs/dev/adding_a_family.md) for detailed guidance on implementing and testing windowing semantics.

## Maintenance notes

TODO: Document maintenance-specific information:
- Key dependencies and version constraints
- Known issues or limitations
- Performance characteristics
- Future enhancement ideas
- Contact information for the family owner
