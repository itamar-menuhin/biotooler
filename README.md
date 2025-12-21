# biotooler

A Python package for bioinformatics tools.

## Installation

```bash
pip install biotooler
```

For development:

```bash
pip install -e ".[dev]"
```

## Quick Start

### Basic Sequence Handling

```python
from biotooler.core.record import coerce_record

# Create a validated DNA sequence
record = coerce_record("ACGTACGT", "DNA", id="seq1", description="My sequence")
print(record.seq)  # ACGTACGT
```

### Sliding Windows

The `iter_windows` function provides efficient, **generic** sliding-window extraction from sequences.
It works with DNA, RNA, and protein sequences without any restrictions on step or window size.
The function uses internal caching to avoid repeatedly converting sequence objects to strings.

**Important:** `iter_windows` is purely generic and does NOT enforce any codon-specific restrictions
(such as step or window_size being multiples of 3). For codon-specific ORF windowing with such
restrictions, use `iter_orf_codon_windows` from `biotooler.core.windowing`.

```python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from biotooler.core.windowing import iter_windows

# Works with DNA sequences
dna_record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

# Extract non-overlapping windows
for window in iter_windows(dna_record, window_size=4, step=4):
    print(f"Window {window.annotations['window_index']}: {window.seq}")
    # Window 0: ACGT
    # Window 1: ACGT
    # Window 2: ACGT

# Extract overlapping windows
for window in iter_windows(dna_record, window_size=4, step=2):
    print(f"Position {window.annotations['start']}-{window.annotations['end']}: {window.seq}")
    # Position 0-4: ACGT
    # Position 2-6: GTAC
    # Position 4-8: TACG
    # Position 6-10: CGTA
    # Position 8-12: GTAC

# Works with protein sequences using step=1 (common for motif detection)
protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein1")
for window in iter_windows(protein_record, window_size=5, step=1):
    print(f"{window.seq}")
    # MKALV
    # KALVS
    # ALVSW
    # LVSWG
    # VSWGR

# Step size doesn't need to be a multiple of 3 (generic windowing)
dna_record2 = SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="seq2")
windows = list(iter_windows(dna_record2, window_size=10, step=4))
# This works fine! 3 windows are generated at positions 0, 4, and 8

# Include partial windows at the end
record = SeqRecord(Seq("ACGTACG"), id="seq3")
windows = list(iter_windows(record, window_size=4, step=2, drop_partial=False))
# Last window will be "CG" even though it's smaller than window_size
```

### ORF Codon Windows

For analyzing features across codon-aligned windows within an ORF, use `iter_orf_codon_windows` 
and `FeatureSet.compute_orf_windows`. This enforces that step_nt is a multiple of 3 and produces
wide-format output with one row per record.

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.windowing import iter_orf_codon_windows
from biotooler.features import FeatureSet

# Define a feature computation function
def compute_gc_content(record):
    seq = str(record.seq)
    gc_count = seq.count('G') + seq.count('C')
    return {"gc_content": gc_count / len(seq) if seq else 0.0}

# Create a FeatureSet
fs = FeatureSet(compute_gc_content, name="gc")

# Compute features over ORF windows (wide format output)
record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
result = fs.compute_orf_windows(
    record,
    orf=(0, 15),  # Explicit ORF coordinates
    window_nt=9,  # 3 codons per window
    step_nt=3     # Step by 1 codon (must be multiple of 3)
)

# Result is a single-row DataFrame with columns:
# record_id, orf_start, orf_end, gc.gc_content_0, gc.gc_content_3, gc.gc_content_6
print(result.shape)  # (1, 6)
print(result.columns.tolist())
# ['record_id', 'orf_start', 'orf_end', 'gc.gc_content_0', 'gc.gc_content_3', 'gc.gc_content_6']
```

### Sequence Utilities

The package provides efficient sequence string/bytes extraction with caching:

```python
from biotooler.core.seq_utils import get_seq_str, get_seq_bytes

# Get normalized sequence string (uppercase, whitespace removed)
# Result is cached for efficiency
seq_str = get_seq_str(record)  # Returns "ACGTACGT"

# Get sequence as bytes (also cached)
seq_bytes = get_seq_bytes(record)  # Returns b"ACGTACGT"

# Cached values are reused on subsequent calls
# This is especially beneficial when processing the same sequence multiple times
```

## Development

### Running Tests

```bash
pytest
```

### Linting

```bash
ruff check .
```

### Type Checking

```bash
pyright
```

## Requirements

- Python >= 3.11
- biopython
- numpy
- pandas