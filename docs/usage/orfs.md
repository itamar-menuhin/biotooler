# Working with ORFs

This guide covers how to work with Open Reading Frames (ORFs) in biotooler, including finding candidates, selecting specific ORFs, and using them in batch workflows.

## Overview

biotooler provides two complementary approaches for working with ORFs:

1. **Single-sequence workflow**: Find ORF candidates and select one by index
2. **Batch workflow**: Pre-compute ORFs and attach them to records for downstream processing

## Finding ORF Candidates

Use `find_orf_candidates()` to enumerate all in-frame start/stop codon pairs in a DNA or RNA sequence:

```python
from biotooler.core.orf_candidates import find_orf_candidates

# Simple sequence with one ORF
seq = "ATGAAATAA"
candidates = find_orf_candidates(seq)
print(candidates)  # [(0, 9)]

# Sequence with multiple ORFs
seq = "ATGAAATAGATGCCCTAGTAA"
candidates = find_orf_candidates(seq)
print(candidates)  # [(0, 9), (0, 18), (0, 21), (9, 18), (9, 21)]
```

**Note:** ORF coordinates are 0-based and end-exclusive (consistent with Python slicing).

## Single-Sequence Workflow

For analyzing a single sequence, find candidates and select one by index:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.core.orf_store import select_orf_by_index

# Create a DNA sequence
record = SeqRecord(Seq("ATGAAATAGATGCCCTAGTAA"), id="seq1")

# Find all ORF candidates
candidates = find_orf_candidates(record)
# [(0, 9), (0, 18), (0, 21), (9, 18), (9, 21)]

# Select the longest ORF (index 2)
orf = select_orf_by_index(candidates, 2)
print(f"Selected ORF: {orf}")  # (0, 21)

# Extract the ORF sequence
start, end = orf
orf_seq = str(record.seq[start:end])
print(f"ORF sequence: {orf_seq}")  # ATGAAATAGATGCCCTAGTAA
```

### Error Handling

`select_orf_by_index()` provides clear error messages for invalid indices:

```python
from biotooler.core.orf_store import select_orf_by_index

candidates = [(0, 9), (9, 18)]

# Index out of range
try:
    select_orf_by_index(candidates, 5)
except IndexError as e:
    print(e)
    # "ORF index 5 is out of range for 2 candidates (valid indices: 0-1)"

# No candidates available
try:
    select_orf_by_index([], 0)
except IndexError as e:
    print(e)
    # "ORF index 0 is out of range: no candidates available"
```

## Batch Workflow

For processing multiple sequences with pre-determined ORFs, use `attach_orf()` to store ORFs in record annotations:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.orf_store import attach_orf, get_orf

# Prepare multiple records with pre-computed ORFs
records = [
    SeqRecord(Seq("ATGAAATAA"), id="seq1"),
    SeqRecord(Seq("ATGCCCTAG"), id="seq2"),
    SeqRecord(Seq("ATGGGGTGA"), id="seq3"),
]

# Pre-determined ORFs for each sequence
orfs = [(0, 9), (0, 9), (0, 9)]

# Attach ORFs to records
for record, orf in zip(records, orfs):
    attach_orf(record, orf)

# Later, retrieve and use the ORFs
for record in records:
    orf = get_orf(record)
    start, end = orf
    orf_seq = str(record.seq[start:end])
    print(f"{record.id}: {orf_seq}")
```

### Custom Annotation Keys

You can use custom keys to store multiple ORFs per record:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.orf_store import attach_orf, get_orf

record = SeqRecord(Seq("ATGAAATAAATGCCCTAG"), id="multi_orf")

# Attach multiple ORFs with different keys
attach_orf(record, (0, 9), key="orf.primary")
attach_orf(record, (9, 18), key="orf.secondary")

# Retrieve specific ORFs
primary = get_orf(record, key="orf.primary")
secondary = get_orf(record, key="orf.secondary")

print(f"Primary ORF: {primary}")    # (0, 9)
print(f"Secondary ORF: {secondary}") # (9, 18)
```

### Error Handling

`get_orf()` provides clear error messages when an ORF is not found:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.orf_store import get_orf

record = SeqRecord(Seq("ATGAAATAA"), id="no_orf_attached")

try:
    get_orf(record)
except KeyError as e:
    print(e)
    # "ORF not found in record 'no_orf_attached': no annotation with key 
    # 'biotooler.orf'. Use attach_orf() to store an ORF first."
```

## Constraints and Validation

### DNA/RNA Only

ORF operations only work with DNA and RNA sequences. Protein sequences will raise a clear error:

```python
from biotooler.core.record import coerce_record
from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.core.orf_store import attach_orf

# Protein sequence
protein_record = coerce_record("MKLVLS", "protein", id="protein1")

# This raises ValueError
try:
    find_orf_candidates(protein_record)
except ValueError as e:
    print(e)
    # "ORF candidate finding is only supported for DNA/RNA sequences, 
    # not protein sequences"

# This also raises ValueError
try:
    attach_orf(protein_record, (0, 6))
except ValueError as e:
    print(e)
    # "ORF operations are only supported for DNA/RNA sequences, 
    # not protein sequences"
```

### Coordinate System

All ORF coordinates use 0-based, end-exclusive indexing (consistent with Python slicing):

```python
from biotooler.core.orf_candidates import find_orf_candidates

seq = "ATGAAATAA"
#     012345678  # positions

candidates = find_orf_candidates(seq)
# [(0, 9)]  # Start at 0, end at 9 (exclusive)

# Direct slicing works
start, end = candidates[0]
assert seq[start:end] == "ATGAAATAA"
```

## API Reference

### `OrfSpan`

Type alias for ORF coordinates: `tuple[int, int]` (start, end)

### `select_orf_by_index(candidates, orf_index)`

Select a single ORF from a list of candidates.

**Parameters:**
- `candidates`: List of ORF spans as `(start, end)` tuples
- `orf_index`: Zero-based index of the ORF to select

**Returns:** The selected ORF span as `(start, end)` tuple

**Raises:** `IndexError` if `orf_index` is out of range

### `attach_orf(record, orf, *, key="biotooler.orf")`

Attach an ORF span to a SeqRecord's annotations.

**Parameters:**
- `record`: SeqRecord to attach the ORF to
- `orf`: ORF span as `(start, end)` tuple
- `key`: Annotation key to store the ORF under (default: `"biotooler.orf"`)

**Returns:** The same SeqRecord with the ORF attached (modified in place)

**Raises:** `ValueError` if the record is a protein sequence

### `get_orf(record, *, key="biotooler.orf")`

Retrieve an ORF span from a SeqRecord's annotations.

**Parameters:**
- `record`: SeqRecord to retrieve the ORF from
- `key`: Annotation key where the ORF is stored (default: `"biotooler.orf"`)

**Returns:** The ORF span as `(start, end)` tuple

**Raises:** `KeyError` if the ORF is not found

## ORF Selection Policies in Translation

When translating DNA/RNA sequences to protein using `ensure_protein_record()`, you can control how ORFs are selected using the `orf_policy` parameter. This is particularly useful for protein-family features where you want to automatically select the best ORF.

### Default Policy

The default policy (`orf_policy="default"`) preserves the existing predictable behavior:

1. Use explicit ORF if provided via `orf` parameter
2. Use attached ORF if present and `use_orf_if_present=True`
3. Translate full sequence from frame 0

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.translation import ensure_protein_record

# Sequence with multiple potential ORFs
dna = SeqRecord(Seq("NNNNATGAAAGCCCTGGTGTAA"), id="seq1")
dna.annotations["molecule_type"] = "DNA"

# Default policy translates from frame 0 (includes leading NNN)
result = ensure_protein_record(dna, orf_policy="default")
print(result.annotations["translation_region_source"])  # "full_sequence_frame0"
```

### Longest ORF Policy

The longest ORF policy (`orf_policy="longest_orf"`) automatically finds and selects the longest ORF candidate:

1. Use explicit ORF if provided via `orf` parameter (highest priority)
2. Use attached ORF if present and `use_orf_if_present=True`
3. **Find all ORF candidates and select the longest one**
4. Fall back to frame 0 if no ORF candidates found

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.translation import ensure_protein_record

# Sequence with 5' UTR followed by coding sequence
dna = SeqRecord(Seq("NNNNATGAAAGCCCTGGTGTAA"), id="seq1")
dna.annotations["molecule_type"] = "DNA"

# Longest ORF policy skips the NNN prefix
result = ensure_protein_record(dna, orf_policy="longest_orf")
print(result.annotations["translation_region"])  # (4, 22)
print(result.annotations["translation_region_source"])  # "longest_orf"
print(str(result.seq))  # "MKALV"
```

### Tie-Breaking Rules

When multiple ORF candidates have the same maximum length, the policy uses deterministic tie-breakers:

1. **Longest length** (primary criterion)
2. **Earliest start position** (first tie-breaker)
3. **Earliest stop position** (second tie-breaker)

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.translation import ensure_protein_record

# Multiple ORFs: (0,9), (12,21), (24,33) - all 9 bases
# The one starting at position 0 wins the tie
dna = SeqRecord(Seq("ATGAAATAAGGGGGGATGCCCTAGGGGGGGATGTTTTAA"), id="seq1")
dna.annotations["molecule_type"] = "DNA"

result = ensure_protein_record(dna, orf_policy="longest_orf")
print(result.annotations["translation_region"])  # (0, 9) - earliest start
```

### When to Use Longest ORF Policy

The `longest_orf` policy is useful when:

- Processing sequences with unknown or variable 5'/3' UTRs
- Analyzing genomic regions where the coding sequence needs to be identified
- Working with protein-family features that require translation from the best ORF
- You want automatic ORF selection without manual intervention

**Important:** This policy is opt-in to maintain predictable default behavior. Always test with your specific use case to ensure the selected ORFs are appropriate.

### Interaction with Explicit and Attached ORFs

The ORF selection priority remains unchanged regardless of policy:

1. **Explicit ORF** (via `orf` parameter) - always takes priority
2. **Attached ORF** (when `use_orf_if_present=True`) - takes priority over policy-based selection
3. **Policy-based selection** - only applies when no explicit or attached ORF is available
4. **Frame 0 fallback** - used when no ORFs can be determined

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.orf_store import attach_orf
from biotooler.core.translation import ensure_protein_record

dna = SeqRecord(Seq("NNNNATGAAAGCCCTGGTGTAA"), id="seq1")
dna.annotations["molecule_type"] = "DNA"

# Explicit ORF overrides longest_orf policy
result = ensure_protein_record(dna, orf=(0, 6), orf_policy="longest_orf")
print(result.annotations["translation_region_source"])  # "explicit_orf"

# Attached ORF overrides longest_orf policy
dna = attach_orf(dna, (0, 6))
result = ensure_protein_record(dna, orf_policy="longest_orf")
print(result.annotations["translation_region_source"])  # "attached_orf"
```

