# Reference Sequences

The `ReferenceSequenceSet` class provides a container for managing CDS (Coding DNA Sequence) and protein sequences with normalization, validation, and translation capabilities.

## Basic Usage

### Creating a ReferenceSequenceSet

```python
from biotooler.core.reference_sequences import ReferenceSequenceSet

# From dictionaries
cds = {"gene1": "ATGAAATAA", "gene2": "ATGGGGTGA"}
proteins = {"gene1": "MK*", "gene2": "MG*"}
ref_set = ReferenceSequenceSet(cds, proteins)

# From FASTA files
ref_set = ReferenceSequenceSet.from_fasta(
    "cds_sequences.fasta",
    protein_fasta_path="protein_sequences.fasta"
)

# From CSV file
ref_set = ReferenceSequenceSet.from_csv(
    "sequences.csv",
    id_column="gene_id",
    cds_column="cds_sequence",
    protein_column="protein_sequence"
)
```

## Validation

The `validate()` method checks sequences for common issues without modifying the instance.

```python
# Validate both CDS and proteins (default)
ref_set.validate()

# Validate only CDS sequences
ref_set.validate(kind="cds")

# Validate only protein sequences
ref_set.validate(kind="protein")
```

### Validation Guarantees

- **CDS validation** (`kind="cds"` or `kind="both"`):
  - No empty sequences (empty strings or whitespace-only)
  - All CDS lengths are multiples of 3 (required for translation)
  
- **Protein validation** (`kind="protein"` or `kind="both"`):
  - No empty sequences (empty strings or whitespace-only)
  - No internal stop codons (stop codons before the last position)

All validation errors include specific information about which sequence failed and why.

## Helper Methods

### with_proteins()

Add or update protein sequences, returning a new instance.

```python
# Add proteins to a CDS-only set
ref_set = ReferenceSequenceSet({"gene1": "ATGAAATAA"})
updated = ref_set.with_proteins({"gene1": "MK*"})

# Merge with existing proteins
ref_set = ReferenceSequenceSet(
    {"gene1": "ATGAAATAA", "gene2": "ATGGGGTGA"},
    {"gene1": "MK*"}
)
updated = ref_set.with_proteins({"gene2": "MG*"})
# Result has proteins for both gene1 and gene2
```

**Guarantees:**
- Returns a new instance; original is unchanged
- No silent changes: raises `ValueError` if a protein ID already exists with a different sequence
- Preserves the original CDS and genetic code table
- Empty sequences (empty strings or whitespace-only) are rejected

### extend_cds()

Add new CDS sequences, returning a new instance.

```python
# Add new CDS sequences
ref_set = ReferenceSequenceSet({"gene1": "ATGAAATAA"})
extended = ref_set.extend_cds({"gene2": "ATGGGGTGA", "gene3": "ATGCCCTAG"})

# Attempting to add a different sequence for the same ID raises an error
try:
    ref_set.extend_cds({"gene1": "ATGGGGTGA"})  # Different sequence!
except ValueError as e:
    print(f"Error: {e}")
```

**Guarantees:**
- Returns a new instance; original is unchanged
- No silent changes: raises `ValueError` if a CDS ID already exists with a different sequence
- Sequences are compared after normalization (uppercase, U→T conversion)
- If the same sequence (after normalization) is provided, the existing sequence format is preserved
- Preserves protein sequences and genetic code table
- Cache is invalidated in the new instance
- Empty sequences (empty strings or whitespace-only) are rejected

### merge()

Merge two `ReferenceSequenceSet` instances.

```python
# Merge disjoint sets
set1 = ReferenceSequenceSet({"gene1": "ATGAAATAA"})
set2 = ReferenceSequenceSet({"gene2": "ATGGGGTGA"})
merged = set1.merge(set2)
# merged contains both gene1 and gene2

# Merge with proteins
set1 = ReferenceSequenceSet({"gene1": "ATGAAATAA"}, {"gene1": "MK*"})
set2 = ReferenceSequenceSet({"gene2": "ATGGGGTGA"}, {"gene2": "MG*"})
merged = set1.merge(set2)
# merged has both CDS and protein sequences
```

**Guarantees:**
- Returns a new instance; originals are unchanged
- Sequences from the first set (self) appear first in the merged result
- No silent changes: raises `ValueError` if the same ID exists in both sets with different sequences
- CDS sequences are compared after normalization (uppercase, U→T conversion)
- Protein sequences are compared exactly (no normalization)
- If the same sequence appears in both sets, the version from the first set is kept
- Genetic code table from the first set is preserved
- If either set has proteins, the merged set will have proteins
- Empty sequences (empty strings or whitespace-only) are rejected

## Sequence Normalization

When accessing sequences through `cds_strings()` and `protein_strings()`:

- **CDS normalization:**
  - Converted to uppercase
  - U (uracil) replaced with T (thymine) for RNA→DNA conversion
  
- **Protein normalization:**
  - Used as-is from provided proteins, or translated from CDS if not provided

Results are cached for efficiency, so repeated calls return the same object reference.

## Error Handling

All helper methods provide explicit, specific error messages:

- Empty sequence errors identify the specific sequence ID
- Duplicate sequence errors show both the existing and provided sequences (first 20 characters)
- Validation errors include the sequence ID, index, and specific issue

**No silent changes:** If there's a conflict or issue, an exception is raised rather than silently modifying or dropping data.
