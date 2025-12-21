# ORF Window Analysis

This guide covers how to compute features across codon-aligned windows within ORFs using the `FeatureSet.compute_orf_windows` method.

## Overview

ORF window analysis allows you to:
- Compute features (e.g., GC content, codon usage) across sliding windows within an ORF
- Get results in **wide format**: one row per record with columns per window (e.g., `CAI_0`, `CAI_3`, `CAI_6`)
- Enforce codon-aligned windowing (step_nt must be multiple of 3)
- Resolve ORF coordinates flexibly (explicit, by index, or attached)

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

# Analyze a sequence
record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

# Compute features over ORF windows (wide format)
result = fs.compute_orf_windows(
    record,
    orf=(0, 15),  # Explicit ORF coordinates
    window_nt=9,  # Window size: 9 nucleotides (3 codons)
    step_nt=3     # Step: 3 nucleotides (1 codon, must be multiple of 3)
)

# Result is a single-row DataFrame
print(result.shape)  # (1, 6)

# Columns: record_id, orf_start, orf_end, gc.gc_content_0, gc.gc_content_3, gc.gc_content_6
print(result.columns.tolist())
# ['record_id', 'orf_start', 'orf_end', 'gc.gc_content_0', 'gc.gc_content_3', 'gc.gc_content_6']

# Access feature values
print(result['gc.gc_content_0'].iloc[0])  # GC content of first window (positions 0-8)
print(result['gc.gc_content_3'].iloc[0])  # GC content of second window (positions 3-11)
print(result['gc.gc_content_6'].iloc[0])  # GC content of third window (positions 6-14)
```

## ORF Resolution

The `compute_orf_windows` method provides three ways to specify ORF coordinates:

### 1. Explicit ORF Coordinates

Provide the ORF coordinates directly using the `orf` parameter:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet

def count_a(record):
    return {"a_count": str(record.seq).count('A')}

fs = FeatureSet(count_a, name="counts")

# Sequence with flanking regions
record = SeqRecord(Seq("NNNNATGAAACCCGGGTTTNNNN"), id="seq1")

# Specify exact ORF region (positions 4-19)
result = fs.compute_orf_windows(
    record,
    orf=(4, 19),  # Explicit coordinates
    window_nt=9,
    step_nt=3
)

print(result['orf_start'].iloc[0])  # 4
print(result['orf_end'].iloc[0])    # 19
```

### 2. By ORF Index

Find all ORF candidates and select one by index:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet

def dummy_feature(record):
    return {"length": len(record.seq)}

fs = FeatureSet(dummy_feature, name="test")

# Sequence with multiple ORFs
record = SeqRecord(Seq("ATGAAATAGATGCCCTAGTAA"), id="seq1")
# This has ORF candidates: (0, 9), (0, 18), (0, 21), (9, 18), (9, 21)

# Select the first ORF candidate (index 0)
result = fs.compute_orf_windows(
    record,
    orf_index=0,  # Select by index
    window_nt=6,
    step_nt=3
)

# Will use ORF (0, 9)
print(result['orf_start'].iloc[0])  # 0
print(result['orf_end'].iloc[0])    # 9
```

### 3. Attached ORF

Use an ORF that was previously attached to the record:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.orf_store import attach_orf
from biotooler.features import FeatureSet

def dummy_feature(record):
    return {"length": len(record.seq)}

fs = FeatureSet(dummy_feature, name="test")

# Attach ORF to record
record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
attach_orf(record, (0, 15))

# No orf or orf_index needed - uses attached ORF
result = fs.compute_orf_windows(record, window_nt=9, step_nt=3)

print(result['orf_start'].iloc[0])  # 0
print(result['orf_end'].iloc[0])    # 15
```

## Multiple Features

Compute multiple features simultaneously:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet

def multi_features(record):
    """Compute multiple features for a sequence window."""
    seq = str(record.seq)
    return {
        "length": len(seq),
        "gc_count": seq.count('G') + seq.count('C'),
        "a_count": seq.count('A'),
        "t_count": seq.count('T')
    }

fs = FeatureSet(multi_features, name="features")
record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

# Result has multiple feature columns for each window
# Columns are sorted by feature name, then window position:
# features.a_count_0, features.a_count_3, features.a_count_6,
# features.gc_count_0, features.gc_count_3, features.gc_count_6,
# features.length_0, features.length_3, features.length_6,
# features.t_count_0, features.t_count_3, features.t_count_6
```

## Output Format

The output is always a **single-row DataFrame** with:

1. **Metadata columns** (first):
   - `record_id`: ID of the sequence
   - `orf_start`: Start position of ORF (0-indexed)
   - `orf_end`: End position of ORF (exclusive)

2. **Feature columns** (deterministically sorted):
   - Format: `{featureset_name}.{feature_key}_{window_start}`
   - Sorted by: feature key alphabetically, then window_start numerically
   - Example: `gc.gc_content_0`, `gc.gc_content_3`, `gc.gc_content_6`

## Partial Windows

Control whether partial windows at the end are included:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet

def count_bases(record):
    return {"length": len(record.seq)}

fs = FeatureSet(count_bases, name="counts")

# 14 nt sequence, windows of 9 nt with step 3 nt
record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

# With drop_partial=True (default): only full-size windows
result = fs.compute_orf_windows(
    record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=True
)
# Only windows at positions 0 and 3 (both size 9)
print([c for c in result.columns if c.startswith('counts.')])
# ['counts.length_0', 'counts.length_3']

# With drop_partial=False: include partial windows
result = fs.compute_orf_windows(
    record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=False
)
# Windows at positions 0, 3, 6, 9, 12 (last ones are partial)
print([c for c in result.columns if c.startswith('counts.')])
# ['counts.length_0', 'counts.length_3', 'counts.length_6', 'counts.length_9', 'counts.length_12']
```

## Error Handling

The method provides clear error messages:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet

def dummy_feature(record):
    return {"feat": 1}

fs = FeatureSet(dummy_feature, name="test")

# Error: Missing ORF information
record = SeqRecord(Seq("ATGAAACCC"), id="seq1")
try:
    result = fs.compute_orf_windows(record, window_nt=9, step_nt=3)
except ValueError as e:
    print(e)  # "No ORF information provided for record 'seq1'..."

# Error: Protein sequences not supported
protein_record = SeqRecord(Seq("MKLVLS"), id="protein1")
protein_record.annotations["molecule_type"] = "protein"
try:
    result = fs.compute_orf_windows(protein_record, orf=(0, 6), window_nt=3, step_nt=3)
except ValueError as e:
    print(e)  # "ORF window computation is only supported for DNA/RNA sequences..."

# Error: step_nt not multiple of 3
record = SeqRecord(Seq("ATGAAACCC"), id="seq1")
try:
    result = fs.compute_orf_windows(record, orf=(0, 9), window_nt=9, step_nt=4)
except ValueError as e:
    print(e)  # "step_nt must be a multiple of 3 for codon-aligned windows..."
```

## API Reference

### `FeatureSet.compute_orf_windows()`

Compute features for ORF windows and return in wide format.

**Parameters:**
- `record` (SeqRecord): DNA/RNA sequence to analyze
- `orf` (tuple[int, int] | None): Explicit ORF coordinates (start, end)
- `orf_index` (int | None): Index to select from ORF candidates
- `window_nt` (int): Window size in nucleotides
- `step_nt` (int): Step size in nucleotides (must be multiple of 3)
- `drop_partial` (bool): Whether to drop partial windows at end (default: True)

**Returns:**
- `pd.DataFrame`: Single-row DataFrame with wide-format features

**Raises:**
- `ValueError`: If step_nt not multiple of 3, protein sequence, or missing ORF info
