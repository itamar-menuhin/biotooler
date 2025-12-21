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

The `iter_windows` function provides efficient sliding-window extraction from sequences.
It uses internal caching to avoid repeatedly converting sequence objects to strings.

```python
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from biotooler.core.windowing import iter_windows

# Create a sequence
record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

# Extract non-overlapping windows
for window in iter_windows(record, window_size=4, step=4):
    print(f"Window {window.annotations['window_index']}: {window.seq}")
    # Window 0: ACGT
    # Window 1: ACGT
    # Window 2: ACGT

# Extract overlapping windows
for window in iter_windows(record, window_size=4, step=2):
    print(f"Position {window.annotations['start']}-{window.annotations['end']}: {window.seq}")
    # Position 0-4: ACGT
    # Position 2-6: GTAC
    # Position 4-8: TACG
    # Position 6-10: CGTA
    # Position 8-12: GTAC

# Include partial windows at the end
record = SeqRecord(Seq("ACGTACG"), id="seq2")
windows = list(iter_windows(record, window_size=4, step=2, drop_partial=False))
# Last window will be "CG" even though it's smaller than window_size
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