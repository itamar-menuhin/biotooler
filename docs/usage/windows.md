# Generic Window Analysis

This guide covers how to compute features across generic sliding windows for any sequence type (DNA/RNA/protein) using the `FeatureSet.compute_windows` method.

## Overview

Generic window analysis allows you to:
- Compute features across sliding windows for DNA, RNA, or protein sequences
- Use any step size (not restricted to codon boundaries like `compute_orf_windows`)
- Get results in **wide format**: one row per record with columns per window (e.g., `gc_content_0`, `gc_content_5`, `gc_content_10`)
- Optionally restrict windowing to a specific region of the sequence
- Work seamlessly with both incremental and non-incremental features

## Comparison with `compute_orf_windows`

| Feature | `compute_windows` | `compute_orf_windows` |
|---------|-------------------|----------------------|
| Sequence types | DNA/RNA/protein | DNA/RNA only |
| Step restrictions | Any positive integer | Must be multiple of 3 (codon-aligned) |
| Region specification | `region=(start, end)` parameter | `orf=(start, end)` or `orf_index` |
| Window numbering | Relative to region start (always starts at _0) | Relative to ORF start |
| Metadata columns | `record_id`, `region_start`, `region_end` | `record_id`, `orf_start`, `orf_end` |

## Basic Usage

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet

# Define a feature computation function
def compute_gc_content(record):
    """Compute GC content of a sequence."""
    seq = str(record.seq)
    gc_count = seq.count('G') + seq.count('C')
    return {"gc_content": gc_count / len(seq) if seq else 0.0}

# Create a FeatureSet
fs = FeatureSet(compute_gc_content, name="gc")

# Analyze a DNA sequence with sliding windows
record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

# Compute features over sliding windows (wide format)
result = fs.compute_windows(
    record,
    window_size=6,  # Window size: 6 bases
    step=3          # Step: 3 bases (50% overlap)
)

# Result is a single-row DataFrame
print(result.shape)  # (1, 6)

# Columns: record_id, region_start, region_end, gc.gc_content_0, gc.gc_content_3, gc.gc_content_6
print(result.columns.tolist())
# ['record_id', 'region_start', 'region_end', 'gc.gc_content_0', 'gc.gc_content_3', 'gc.gc_content_6']

# Access feature values
print(result['gc.gc_content_0'].iloc[0])  # GC content of first window (positions 0-5)
print(result['gc.gc_content_3'].iloc[0])  # GC content of second window (positions 3-8)
```

## Working with Protein Sequences

Unlike `compute_orf_windows`, `compute_windows` works seamlessly with protein sequences and supports step=1:

```python
# Analyze a protein sequence with sliding windows
protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein1")
protein_record.annotations["molecule_type"] = "protein"

# Use step=1 for true sliding window over amino acids
result = fs.compute_windows(
    protein_record,
    window_size=5,  # 5 amino acid window
    step=1          # Slide by 1 amino acid
)

# Creates windows at positions 0, 1, 2, 3, 4
print([c for c in result.columns if c.startswith('gc.')])
# ['gc.gc_content_0', 'gc.gc_content_1', 'gc.gc_content_2', 'gc.gc_content_3', 'gc.gc_content_4']
```

## Using the Region Parameter

The `region` parameter allows you to restrict windowing to a specific subsequence:

```python
# Analyze only a specific region of the sequence
full_record = SeqRecord(Seq("NNNNACGTACGTNNNN"), id="seq1")

# Only window over positions 4-12 (ACGTACGT)
result = fs.compute_windows(
    full_record,
    window_size=4,
    step=2,
    region=(4, 12)  # Only analyze this subsequence
)

# Windows are numbered relative to region start
# Window 0 is at position 4-8 in the full sequence
# Window 2 is at position 6-10 in the full sequence
# etc.

print(result['region_start'].iloc[0])  # 4
print(result['region_end'].iloc[0])    # 12
```

## Output Format

### Wide Format

`compute_windows` returns a single-row DataFrame in wide format:

```python
result = fs.compute_windows(record, window_size=6, step=3)
print(result.shape)  # (1, N) - single row, N columns
```

### Column Naming

Columns follow this pattern: `{set_name}.{feature_key}_{window_start}`

```python
fs = FeatureSet(compute_multi_feature, name="stats")
result = fs.compute_windows(record, window_size=6, step=3)

# Columns: stats.gc_content_0, stats.gc_content_3, stats.length_0, stats.length_3, etc.
```

### Column Ordering

Columns are deterministically ordered:
1. **Metadata columns first**: `record_id`, `region_start`, `region_end`
2. **Feature columns sorted by**: feature key (alphabetically), then window start (numerically)

```python
# Example ordering:
# ['record_id', 'region_start', 'region_end',
#  'stats.gc_content_0', 'stats.gc_content_3', 'stats.gc_content_6',  # gc_content sorted
#  'stats.length_0', 'stats.length_3', 'stats.length_6']               # then length
```

## Handling Partial Windows

Use the `drop_partial` parameter to control what happens with partial windows at the end:

```python
record = SeqRecord(Seq("ACGTACGTACG"), id="seq1")  # 11 bases

# drop_partial=True (default): only full windows
result_full = fs.compute_windows(
    record,
    window_size=4,
    step=2,
    drop_partial=True  # Default
)
# Windows: 0-4, 2-6, 4-8, 6-10 (4 full windows)

# drop_partial=False: include partial windows
result_all = fs.compute_windows(
    record,
    window_size=4,
    step=2,
    drop_partial=False
)
# Windows: 0-4, 2-6, 4-8, 6-10, 8-11, 10-11 (includes partial windows of size 3 and 1)
```

## Incremental Features

`compute_windows` supports incremental features, which can significantly improve performance for overlapping windows:

```python
class IncrementalGCContent:
    """Incremental GC content computation."""
    
    def init_state(self, record, *, window_start, window_end, **kwargs):
        """Initialize state for first window."""
        seq = str(record.seq)[window_start:window_end].upper()
        return {
            'record': record,
            'gc_count': seq.count('G') + seq.count('C'),
            'window_size': window_end - window_start
        }
    
    def step_state(self, state, *, out_start, out_end, in_start, in_end, **kwargs):
        """Update state incrementally."""
        record = state['record']
        # Remove outgoing bases
        if out_end > out_start:
            out_seq = str(record.seq)[out_start:out_end].upper()
            state['gc_count'] -= (out_seq.count('G') + out_seq.count('C'))
        # Add incoming bases
        if in_end > in_start:
            in_seq = str(record.seq)[in_start:in_end].upper()
            state['gc_count'] += (in_seq.count('G') + in_seq.count('C'))
    
    def emit(self, state):
        """Emit feature values."""
        return {
            'gc_content': state['gc_count'] / state['window_size']
        }

# Use incremental feature
fs = FeatureSet(IncrementalGCContent(), name="gc")
result = fs.compute_windows(record, window_size=100, step=1)
# Much faster for large overlapping windows!
```

## Multiple Features

You can compute multiple features simultaneously:

```python
def compute_stats(record):
    """Compute multiple sequence statistics."""
    seq = str(record.seq)
    return {
        "length": len(seq),
        "gc_count": seq.count('G') + seq.count('C'),
        "at_count": seq.count('A') + seq.count('T')
    }

fs = FeatureSet(compute_stats, name="stats")
result = fs.compute_windows(record, window_size=6, step=3)

# Features are grouped by key in column ordering:
# stats.at_count_0, stats.at_count_3, stats.at_count_6,
# stats.gc_count_0, stats.gc_count_3, stats.gc_count_6,
# stats.length_0, stats.length_3, stats.length_6
```

## Common Use Cases

### Sliding Window Analysis for Proteins

```python
# Analyze hydrophobicity along a protein sequence
protein = SeqRecord(Seq("MKALVSWGRGQPVQVFLD"), id="protein1")
protein.annotations["molecule_type"] = "protein"

fs = FeatureSet(compute_hydrophobicity, name="hydro")
result = fs.compute_windows(
    protein,
    window_size=7,  # 7 residue window
    step=1          # Slide by 1 residue
)
```

### DNA/RNA Analysis with Flexible Windows

```python
# Analyze GC content with 50% overlapping windows
dna = SeqRecord(Seq("ACGTACGTACGTACGT"), id="dna1")

fs = FeatureSet(compute_gc_content, name="gc")
result = fs.compute_windows(
    dna,
    window_size=10,
    step=5  # Not restricted to multiples of 3
)
```

### Region-Specific Analysis

```python
# Analyze only the promoter region
# Note: Using 'N' * 1000 as placeholder for a real gene sequence
full_gene = SeqRecord(Seq("N" * 1000), id="gene1")

# Focus on positions 0-500 (promoter region)
result = fs.compute_windows(
    full_gene,
    window_size=50,
    step=10,
    region=(0, 500)  # Only analyze promoter
)
```

## Notes

- Window starts in column names are always relative to the region start (or sequence start if no region specified)
- The first window always has suffix `_0` regardless of the region parameter
- For codon-aligned ORF analysis with step restricted to multiples of 3, use `compute_orf_windows` instead
- Incremental features can provide significant performance benefits for large windows with small steps
